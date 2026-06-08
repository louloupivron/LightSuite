function bad = validateSpinalCordSliceSeries(dp)
%VALIDATESPINALCORDSLICESERIES Check all slice TIFFs in a folder can be read.
%   bad = validateSpinalCordSliceSeries(dp)
%   Returns a table (sliceIndex, fileName, filePath, message). Empty if all OK.

tiffiles = dir(fullfile(dp, '*.tiff'));
tifiles  = dir(fullfile(dp, '*.tif'));
tfiles   = cat(1, tifiles, tiffiles);
assert(~isempty(tfiles), 'No .tif files in %s', dp);

fprintf('Scanning slice TIFFs in %s\n', dp);
[~, bad] = filterSpinalCordSlices(tfiles, true);
if isempty(bad)
    [~, sortIdx] = natsortfiles({tfiles.name});
    tfiles = tfiles(sortIdx);
    tinfo = imfinfo(fullfile(tfiles(1).folder, tfiles(1).name));
    fprintf('All %d slices OK (%dx%d px).\n', numel(tfiles), tinfo.Height, tinfo.Width);
else
    fprintf('Found %d unreadable slice(s).\n', height(bad));
end
end
