function [tfilesOut, bad] = filterSpinalCordSlices(tfiles, skipCorrupt)
%FILTERSPINALCORDSLICES Keep readable slice TIFFs from a plane-per-file series.
%   [tfilesOut, bad] = filterSpinalCordSlices(tfiles)
%   [tfilesOut, bad] = filterSpinalCordSlices(tfiles, skipCorrupt)
%
%   tfiles is a dir() struct array. Files are natsort-ed. Each file is tested
%   with readPlaneTiff. If skipCorrupt is false (default) and any file fails,
%   throws with a list of all corrupt paths. If skipCorrupt is true, bad files
%   are removed from tfilesOut and listed in table bad.

if nargin < 2
    skipCorrupt = false;
end

[~, sortIdx] = natsortfiles({tfiles.name});
tfiles = tfiles(sortIdx);
folder = tfiles(1).folder;
Nz = numel(tfiles);

tinfo = imfinfo(fullfile(folder, tfiles(1).name));
Ny0 = tinfo.Height;
Nx0 = tinfo.Width;

bad = table([], [], [], [], ...
    'VariableNames', {'sliceIndex', 'fileName', 'filePath', 'message'});
goodMask = true(Nz, 1);

fprintf('  validating %d slice TIFFs...\n', Nz);
for iz = 1:Nz
    path = fullfile(folder, tfiles(iz).name);
    try
        currim = readPlaneTiff(path);
        if ~isequal(size(currim, 1:2), [Ny0 Nx0])
            msg = sprintf('size %dx%d px, expected %dx%d px', ...
                size(currim, 1), size(currim, 2), Ny0, Nx0);
            error('filterSpinalCordSlices:SizeMismatch', '%s', msg);
        end
    catch ME
        goodMask(iz) = false;
        bad = [bad; {iz, tfiles(iz).name, path, ME.message}]; %#ok<AGROW>
    end
    if mod(iz, 200) == 0 || iz == Nz
        fprintf('    checked %d / %d (%d failed)\n', iz, Nz, height(bad));
    end
end

if ~isempty(bad) && ~skipCorrupt
    raiseCorruptSpinalCordSliceError(bad, Nz);
end

if ~isempty(bad) && skipCorrupt
    fprintf(['  warning: skipping %d corrupt slice(s); loading %d of %d\n'], ...
        height(bad), sum(goodMask), Nz);
end

tfilesOut = tfiles(goodMask);
if isempty(tfilesOut)
    error('readSpinalCordSample:NoReadableSlices', ...
        'No readable slice TIFFs remain in %s after filtering.', folder);
end
end
