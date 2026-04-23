function stats = invertElastixTransformCP(transformDir, outputDir)
%INVERTELASTIXTRANSFORMCP Invert an elastix transform using elastix.log (LightSuite fork).
%
%   Same role as matlab_elastix invertElastixTransformCP / invertElastixTransform, with
%   parsing fixes for Elastix 5.1 logs: -f0 … -m0 with or without quotes (Windows logs),
%   a single -p flag (quoted or not), deduplication of repeated -p lines, and no false
%   matches on -fp / -mp. If the log omits the CLI, falls back to *_dual_f0.mhd in the
%   same folder (LightSuite dual-channel temp layout). Parameter file path can also be
%   taken from elastix lines "end of ParameterFile: <path>".
%
%   Forwards fixedscale/movingscale to elastix from ElementSpacing in the fixed MHD (same
%   idea as performMultObjBsplineRegistration) so temporary MHDs are not written with zero spacing.
%   If the parameter file includes CorrespondingPointsEuclideanDistanceMetric, passes
%   fixed.txt / moving.txt from transformDir as fixedpoints/movingpoints (-fp/-mp).
%
%   See also: elastix, elastixDualFixedSameMovingBspline, clearElastixWorkspaceForNewRun

% Based on Rob Campbell's invertElastixTransform (matlab_elastix), modified for LightSuite.

if nargin < 2
    outputDir = [];
end

addElastixRepoPaths;

logFile = fullfile(transformDir, 'elastix.log');
if exist(logFile, 'file') ~= 2
    error('LightSuite:invertElastixTransformCP:NoLog', ...
        'Can not find elastix.log in %s', transformDir);
end

txt = fileread(logFile);

% Fixed image: elastix logs on Windows often omit quotes; parse -f0 … -m0 / -f … -m.
fixedFile = localParseFixedImagePath(txt, transformDir);
if isempty(fixedFile)
    error('LightSuite:invertElastixTransformCP:NoFixed', ...
        ['Could not find -f0/-f fixed image path in %s. ', ...
        'If the log omits the CLI, ensure *_dual_f0.mhd exists in that folder.'], logFile);
end
if exist(fixedFile, 'file') ~= 2
    error('LightSuite:invertElastixTransformCP:FixedMissing', ...
        'Can not find fixed file at %s', fixedFile);
end

% Parameter files: -p path but not -fp / -mp (quoted or unquoted)
paths = localParseParameterFilePaths(txt);
if isempty(paths)
    error('LightSuite:invertElastixTransformCP:NoParams', ...
        ['Could not find parameter file path (-p or "end of ParameterFile:") in %s'], logFile);
end
params = unique(paths, 'stable');

for ii = 1:numel(params)
    if exist(params{ii}, 'file') ~= 2
        error('LightSuite:invertElastixTransformCP:ParamMissing', ...
            'Can not find parameter file at %s', params{ii});
    end
end

files = dir(fullfile(transformDir, 'TransformParameters.*.txt'));
files = files(~[files.isdir]);
if numel(files) ~= numel(params)
    error('LightSuite:invertElastixTransformCP:CountMismatch', ...
        ['Found %d TransformParameters*.txt but %d parameter file(s) from log. ', ...
        'Remove stale files in %s (see clearElastixWorkspaceForNewRun).'], ...
        numel(files), numel(params), transformDir);
end

coefFiles = fliplr({files.name});
for ii = 1:numel(coefFiles)
    coefFiles{ii} = fullfile(transformDir, coefFiles{ii});
    if exist(coefFiles{ii}, 'file') ~= 2
        error('LightSuite:invertElastixTransformCP:CoefMissing', ...
            'Can not find coef file at %s', coefFiles{ii});
    end
end

fprintf('Using fixed file: %s\n', fixedFile);
fprintf('Using parameter files:');
for ii = 1:numel(params)
    fprintf(' %s', params{ii});
    if ii < numel(params)
        fprintf(',');
    end
end
fprintf('\n');
fprintf('Using coef files:');
for ii = 1:numel(coefFiles)
    fprintf(' %s', coefFiles{ii});
    if ii < numel(coefFiles)
        fprintf(',');
    end
end
fprintf('\n');

fixedImage = mhd_read(fixedFile);
sp = localReadMhdElementSpacing(fixedFile);
nd = ndims(fixedImage);
if isempty(sp) || any(sp <= 0) || any(~isfinite(sp))
    error('LightSuite:invertElastixTransformCP:BadSpacing', ...
        ['Could not read positive ElementSpacing from %s. ', ...
        'Inversion needs the same voxel spacing (mm) as the forward B-spline run.'], fixedFile);
end
sp = sp(:).';
if numel(sp) < nd
    if isscalar(sp)
        sp = repmat(sp, 1, nd);
    else
        error('LightSuite:invertElastixTransformCP:SpacingDims', ...
            'ElementSpacing in %s has %d values but image has %d dimensions.', ...
            fixedFile, numel(sp), nd);
    end
elseif numel(sp) > nd
    sp = sp(1:nd);
end

% Multi-metric B-spline params (e.g. dual MI + landmarks) still include
% CorrespondingPointsEuclideanDistanceMetric; elastix needs -fp/-mp for that slot.
fixpath = fullfile(transformDir, 'fixed.txt');
movpath = fullfile(transformDir, 'moving.txt');
needsLandmarks = localParamFileNeedsLandmarks(params);
if needsLandmarks && (exist(fixpath, 'file') ~= 2 || exist(movpath, 'file') ~= 2)
    error('LightSuite:invertElastixTransformCP:LandmarkFilesMissing', ...
        ['Parameter file uses CorrespondingPointsEuclideanDistanceMetric but ', ...
        'fixed.txt / moving.txt were not found in %s. Keep elastix_temp until inversion finishes.'], ...
        transformDir);
end
elastixArgs = {'t0', coefFiles, 'fixedscale', sp, 'movingscale', sp};
if exist(fixpath, 'file') == 2 && exist(movpath, 'file') == 2
    fprintf('Using landmark files: %s , %s\n', fixpath, movpath);
    elastixArgs = [elastixArgs, {'fixedpoints', fixpath, 'movingpoints', movpath}];
end

[~, stats] = elastix(fixedImage, fixedImage, outputDir, params, elastixArgs{:});

if ~isa(stats, 'struct') || ~isfield(stats, 'TransformParameters') || ...
        ~iscell(stats.TransformParameters) || isempty(stats.TransformParameters)
    error('LightSuite:invertElastixTransformCP:ElastixFailed', ...
        ['elastix inversion did not return TransformParameters. ', ...
        'See elastix.log in %s and the elastix console output above.'], outputDir);
end

stats.TransformParameters{1}.InitialTransformParametersFileName = 'NoInitialTransform';
end

%--------------------------------------------------------------------------
function fixedFile = localParseFixedImagePath(txt, transformDir)
% Prefer log delimiters (-f0 … -m0) so quoted and unquoted Windows paths work.
fixedFile = '';
toks = regexp(txt, '-f0\s+(.+?)\s+-m0\b', 'tokens');
if isempty(toks)
    toks = regexp(txt, '-f\s+(.+?)(?=\s+(?<![fm])-m\b)', 'tokens');
end
if ~isempty(toks)
    fixedFile = localStripOuterQuotes(strtrim(toks{end}{1}));
end
if isempty(fixedFile) || exist(fixedFile, 'file') ~= 2
    fixedFile = localGuessDualFixedMhd(transformDir);
end
end

%--------------------------------------------------------------------------
function fixedFile = localGuessDualFixedMhd(transformDir)
fixedFile = '';
d = [dir(fullfile(transformDir, '*_dual_f0.mhd')); dir(fullfile(transformDir, '*_dual_f0.MHD'))];
d = d(~[d.isdir]);
if isempty(d)
    return
end
[~, ui] = unique({d.name}, 'stable');
d = d(ui);
[~, ix] = max([d.datenum]);
fixedFile = fullfile(transformDir, d(ix).name);
end

%--------------------------------------------------------------------------
function s = localStripOuterQuotes(s)
if numel(s) < 2
    return
end
if s(1) == '"' && s(end) == '"'
    s = s(2:end-1);
elseif s(1) == '''' && s(end) == ''''
    s = s(2:end-1);
end
end

%--------------------------------------------------------------------------
function paths = localParseParameterFilePaths(txt)
% Elastix logs the parameter path after "end of ParameterFile:" even when the CLI is
% not echoed. Also parse -p (quoted or unquoted); use separate single-capture regexes
% because alternation token layouts differ across MATLAB releases.
paths = {};
pfTok = [regexp(txt, 'end of ParameterFile:\s*([^\r\n]+)', 'tokens'), ...
    regexp(txt, 'begin of ParameterFile:\s*([^\r\n]+)', 'tokens')];
for k = 1:numel(pfTok)
    t = pfTok{k};
    p = localFirstNonemptyToken(t);
    if isempty(p)
        continue
    end
    p = strtrim(p);
    p = regexprep(p, '\s*=+\s*$', '');
    p = localStripOuterQuotes(p);
    if ~isempty(p)
        paths{end+1} = p; %#ok<AGROW>
    end
end
% Quoted -p (not -fp / -mp)
toksQ = regexp(txt, '(?<![fm])-p\s+"([^"]+)"', 'tokens');
for k = 1:numel(toksQ)
    p = localFirstNonemptyToken(toksQ{k});
    if ~isempty(p)
        paths{end+1} = p; %#ok<AGROW>
    end
end
% Unquoted -p: next token must not start with " (those are handled above)
toksU = regexp(txt, '(?<![fm])-p\s+(?!")(\S+)', 'tokens');
for k = 1:numel(toksU)
    p = localFirstNonemptyToken(toksU{k});
    if isempty(p)
        continue
    end
    p = localStripOuterQuotes(p);
    if ~isempty(p)
        paths{end+1} = p; %#ok<AGROW>
    end
end
end

%--------------------------------------------------------------------------
function p = localFirstNonemptyToken(t)
p = '';
if isempty(t) || ~iscell(t)
    return
end
for ii = 1:numel(t)
    if ~isempty(t{ii})
        p = t{ii};
        return
    end
end
end

%--------------------------------------------------------------------------
function sp = localReadMhdElementSpacing(mhdPath)
% Read ElementSpacing (or ElementSize) from an MHD header; values are in mm as written by mhd_write.
sp = [];
if exist(mhdPath, 'file') ~= 2
    return
end
txt = fileread(mhdPath);
tok = regexp(txt, '(?i)ElementSpacing\s*=\s*([^\r\n]+)', 'tokens', 'once');
if isempty(tok)
    tok = regexp(txt, '(?i)ElementSize\s*=\s*([^\r\n]+)', 'tokens', 'once');
end
if isempty(tok)
    return
end
line = strtrim(tok{1});
if contains(line, '#')
    line = strtrim(strtok(line, '#'));
end
if contains(line, '%')
    line = strtrim(strtok(line, '%'));
end
line = strrep(line, ',', ' ');
nums = sscanf(line, '%f');
if isempty(nums) || any(~isfinite(nums)) || any(nums <= 0)
    return
end
sp = nums(:).';
end

%--------------------------------------------------------------------------
function tf = localParamFileNeedsLandmarks(paramPaths)
tf = false;
for ii = 1:numel(paramPaths)
    p = paramPaths{ii};
    if exist(p, 'file') ~= 2
        continue
    end
    t = fileread(p);
    if contains(t, 'CorrespondingPointsEuclideanDistanceMetric')
        tf = true;
        return
    end
end
end
