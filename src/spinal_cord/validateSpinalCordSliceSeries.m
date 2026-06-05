function bad = validateSpinalCordSliceSeries(dp)
%VALIDATESPINALCORDSLICESERIES Check all slice TIFFs in a folder can be read.
%   bad = validateSpinalCordSliceSeries(dp)
%   Returns a table with fields sliceIndex, fileName, filePath, message.
%   Empty table means all slices read OK.
%
%   Example:
%     bad = validateSpinalCordSliceSeries(dpspinesample);
%     if ~isempty(bad), disp(bad); end

tiffiles = dir(fullfile(dp, '*.tiff'));
tifiles  = dir(fullfile(dp, '*.tif'));
tfiles   = cat(1, tifiles, tiffiles);
assert(~isempty(tfiles), 'No .tif files in %s', dp);

[~, sortIdx] = natsortfiles({tfiles.name});
tfiles = tfiles(sortIdx);
folder = tfiles(1).folder;
Nz = numel(tfiles);

try
    tinfo = imfinfo(fullfile(folder, tfiles(1).name));
    Ny0 = tinfo.Height;
    Nx0 = tinfo.Width;
catch ME
    error('validateSpinalCordSliceSeries:FirstSlice', ...
        'Could not read first slice %s: %s', tfiles(1).name, ME.message);
end

bad = table([], [], [], [], ...
    'VariableNames', {'sliceIndex', 'fileName', 'filePath', 'message'});
fprintf('Validating %d slice TIFFs in %s\n', Nz, dp);
for iz = 1:Nz
    path = fullfile(folder, tfiles(iz).name);
    try
        currim = readPlaneTiff(path);
        if ~isequal(size(currim, 1:2), [Ny0 Nx0])
            msg = sprintf('size %dx%d, expected %dx%d', ...
                size(currim, 1), size(currim, 2), Ny0, Nx0);
            bad = [bad; {iz, tfiles(iz).name, path, msg}]; %#ok<AGROW>
        end
    catch ME
        bad = [bad; {iz, tfiles(iz).name, path, ME.message}]; %#ok<AGROW>
    end
    if mod(iz, 200) == 0 || iz == Nz
        fprintf('  checked %d / %d (%d failed so far)\n', iz, Nz, height(bad));
    end
end
if isempty(bad)
    fprintf('All %d slices OK (%dx%d px).\n', Nz, Ny0, Nx0);
else
    fprintf('Found %d unreadable slice(s). See table output.\n', height(bad));
end
end
