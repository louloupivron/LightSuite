function stats = invertElastixTransformCP(transformDir, outputDir)
%INVERTELASTIXTRANSFORMCP Invert an elastix transform using elastix.log (LightSuite fork).
%
%   Same role as matlab_elastix invertElastixTransformCP / invertElastixTransform, with
%   parsing fixes for Elastix 5.1 logs that use -f0/-f1/-f2 and a single -p flag (one
%   space), and deduplication of repeated -p lines in the log. Also avoids matching -fp
%   / -mp as parameter flags.
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

% Fixed image: prefer -f0 (multi-image CLI), else legacy -f
tok = regexp(txt, '-f0\s+"([^"]+)"', 'tokens', 'once');
if isempty(tok)
    tok = regexp(txt, '-f\s+"([^"]+)"', 'tokens', 'once');
end
if isempty(tok)
    error('LightSuite:invertElastixTransformCP:NoFixed', ...
        'Could not find -f0 or -f fixed image path in %s', logFile);
end
fixedFile = tok{1};
if exist(fixedFile, 'file') ~= 2
    error('LightSuite:invertElastixTransformCP:FixedMissing', ...
        'Can not find fixed file at %s', fixedFile);
end

% Parameter files: -p "path" but not -fp / -mp (negative lookbehind before "p")
allTok = regexp(txt, '(?<![fm])-p\s+"([^"]+)"', 'tokens');
if isempty(allTok)
    error('LightSuite:invertElastixTransformCP:NoParams', ...
        'Could not find -p parameter file path in %s', logFile);
end
paths = cellfun(@(c)c{1}, allTok, 'UniformOutput', false);
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
