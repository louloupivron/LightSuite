function stats = invertElastixTransformCP(transformDir, outputDir)
%INVERTELASTIXTRANSFORMCP Invert an elastix transform using elastix.log (LightSuite fork).
%
%   Same role as matlab_elastix invertElastixTransformCP / invertElastixTransform, with
%   parsing fixes for Elastix 5.1 logs: -f0 … -m0 with or without quotes (Windows logs),
%   a single -p flag (quoted or not), deduplication of repeated -p lines, and no false
%   matches on -fp / -mp. If the log omits the CLI, falls back to *_dual_f0.mhd in the
%   same folder (LightSuite dual-channel temp layout).
%
%   See also: elastix, elastixDualFixedSameMovingBspline, clearElastixWorkspaceForNewRun

% Based on Rob Campbell's invertElastixTransform (matlab_elastix), modified for LightSuite.

if nargin < 2
    outputDir = [];
end

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
        'Could not find -p parameter file path in %s', logFile);
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
[~, stats] = elastix(fixedImage, fixedImage, outputDir, params, 't0', coefFiles);

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
% Quoted "-p ""...""", or unquoted -p F:\path\file.txt (not -fp / -mp).
toks = regexp(txt, '(?<![fm])-p\s+(?:"([^"]+)"|(\S+))', 'tokens');
paths = {};
for k = 1:numel(toks)
    t = toks{k};
    if ~isempty(t{1})
        paths{end+1} = t{1}; %#ok<AGROW>
    elseif numel(t) >= 2 && ~isempty(t{2})
        paths{end+1} = t{2}; %#ok<AGROW>
    end
end
end
