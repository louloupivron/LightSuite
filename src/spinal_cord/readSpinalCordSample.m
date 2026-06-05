function [finvol, opts] = readSpinalCordSample(dp, sampleres, tifftype, registrationres)
%READSPINALCORDSAMPLE Load a spinal cord lightsheet TIFF stack.
%   [finvol, opts] = readSpinalCordSample(dp, sampleres)
%   [finvol, opts] = readSpinalCordSample(dp, sampleres, tifftype, registrationres)
%
%   dp is a folder containing TIFF files. Supported layouts:
%     'planeperfile'   — one 2D TIFF per slice (Terastitcher-style series)
%     'channelperfile' — one multi-page TIFF, or one stack file per channel
%     'auto' (default) — infer from file count and page count
%
%   sampleres must be the native voxel size in micrometers (not the registration
%   grid). For large slice series, data are downsampled while loading to
%   registrationres (default 20 x 20 x 20 um) so the full stack is never held
%   in memory at native resolution.
%
%   finvol is [height x width x depth x channels] in uint16.

if nargin < 3 || isempty(tifftype)
    tifftype = 'auto';
end
if nargin < 4 || isempty(registrationres)
    registrationres = [20 20 20];
end
tifftype = lower(strtrim(char(tifftype)));

tiffiles = dir(fullfile(dp, '*.tiff'));
tifiles  = dir(fullfile(dp, '*.tif'));
tfiles   = cat(1, tifiles, tiffiles);
if isempty(tfiles)
    error('readSpinalCordSample:NoTiffs', ...
        'No .tif or .tiff files found in:\n  %s', dp);
end

layout = localResolveTiffLayout(tfiles, tifftype);
nativeSize = [];
tic;
switch layout
    case 'planeperfile'
        fprintf('Using planeperfile loading (%d slice TIFFs)\n', numel(tfiles));
        [finvol, nativeSize] = localLoadPlanePerFileStack(tfiles, sampleres, registrationres);
    case 'channelperfile'
        if numel(tfiles) == 1
            finvol = localLoadSingleTiffStack(tfiles(1));
        else
            finvol = localLoadChannelPerFileStack(tfiles);
        end
        nativeSize = size(finvol, 1:3);
    otherwise
        error('readSpinalCordSample:BadLayout', ...
            'Unknown tifftype ''%s''. Use ''auto'', ''planeperfile'', or ''channelperfile''.', tifftype);
end

[Nslices, Ny, Nx, Nchan] = size(finvol);
assert(all([Nslices, Ny, Nx, Nchan] > 0), ...
    'readSpinalCordSample:EmptyVolume', ...
    'Loaded volume is empty (file: %s). Check that Bio-Formats can read this TIFF.', ...
    fullfile(tfiles(1).folder, tfiles(1).name));
assert(Nslices > 1 && Ny > 1 && Nx > 1, ...
    'readSpinalCordSample:DegenerateStack', ...
    ['Loaded stack is %d x %d x %d. Expected a 3D volume with depth > 1 along ' ...
    'the cord. If depth is 1, Bio-Formats may not have parsed Z/T correctly.'], ...
    Nslices, Ny, Nx);
fprintf('Parsed spinal cord sample in %2.1f s. Size %d x %d x %d with %d channels\n', ...
    toc, Nslices, Ny, Nx, Nchan)
if ~isempty(nativeSize) && ~isequal(nativeSize, [Nslices Ny Nx])
    fprintf('  (native size was %d x %d x %d px; downsampled to registration grid)\n', ...
        nativeSize(1), nativeSize(2), nativeSize(3));
end
%--------------------------------------------------------------------------
opts.datafolder     = tfiles(1).folder;
opts.lsfolder       = fullfile(opts.datafolder, 'lightsuite');
makeNewDir(opts.lsfolder)
opts.orisize         = nativeSize;
if isempty(opts.orisize)
    opts.orisize = [Nslices, Ny, Nx];
end
opts.Nchan           = Nchan;
opts.sampleres       = sampleres;
opts.registrationres = registrationres;
opts.tifftype        = layout;
%--------------------------------------------------------------------------

end

%--------------------------------------------------------------------------
function layout = localResolveTiffLayout(tfiles, tifftype)
valid = {'auto', 'planeperfile', 'channelperfile'};
assert(any(strcmp(tifftype, valid)), ...
    'readSpinalCordSample:BadTiffType', ...
    'tifftype must be ''auto'', ''planeperfile'', or ''channelperfile''.');
if ~strcmp(tifftype, 'auto')
    layout = tifftype;
    return
end
if numel(tfiles) == 1
    layout = 'channelperfile';
    return
end
if localLooksLikePlanePerFile(tfiles)
    layout = 'planeperfile';
else
    layout = 'channelperfile';
end
end

%--------------------------------------------------------------------------
function tf = localLooksLikePlanePerFile(tfiles)
nCheck = min(numel(tfiles), 8);
idxCheck = unique(round(linspace(1, numel(tfiles), nCheck)));
for k = idxCheck
    try
        nPages = numel(imfinfo(fullfile(tfiles(k).folder, tfiles(k).name)));
    catch
        tf = false;
        return
    end
    if nPages ~= 1
        tf = false;
        return
    end
end
tf = true;
end

%--------------------------------------------------------------------------
function [finvol, nativeSize] = localLoadPlanePerFileStack(tfiles, sampleres, registrationres)
[~, sortIdx] = natsortfiles({tfiles.name});
tfiles = tfiles(sortIdx);
folder = tfiles(1).folder;
Nz = numel(tfiles);
tinfo = imfinfo(fullfile(folder, tfiles(1).name));
Ny0 = tinfo.Height;
Nx0 = tinfo.Width;
nativeSize = [Ny0, Nx0, Nz];

[resfac, targetSize] = cordVolumeDownsampleSpec(nativeSize, sampleres, registrationres);
nativeBytes = double(Ny0) * double(Nx0) * double(Nz) * 2;
localAssertPlanePerFileLoadable(nativeSize, sampleres, registrationres, nativeBytes, resfac);

if isequal(targetSize, nativeSize)
    finvol = zeros(Ny0, Nx0, Nz, 1, 'uint16');
    for iz = 1:Nz
        finvol(:, :, iz, 1) = imread(fullfile(folder, tfiles(iz).name));
        localPrintSliceProgress(iz, Nz);
    end
    return
end

fprintf(['  downsampling while loading: %d x %d x %d -> %d x %d x %d px ' ...
    '(sampleres %s um -> registration %s um)\n'], ...
    Ny0, Nx0, Nz, targetSize(1), targetSize(2), targetSize(3), ...
    mat2str(sampleres), mat2str(registrationres));

targetNy = targetSize(1);
targetNx = targetSize(2);
targetNz = targetSize(3);
backvol = zeros(targetNy, targetNx, Nz, 'uint16');
for iz = 1:Nz
    currim = imread(fullfile(folder, tfiles(iz).name));
    if ~isequal(size(currim, 1:2), [Ny0 Nx0])
        error('readSpinalCordSample:SliceSizeMismatch', ...
            'Slice %d (%s) is %dx%d px, expected %dx%d px.', ...
            iz, tfiles(iz).name, size(currim, 1), size(currim, 2), Ny0, Nx0);
    end
    if resfac(1) == 1 && resfac(2) == 1
        backvol(:, :, iz) = currim;
    else
        backvol(:, :, iz) = imresize(currim, [targetNy targetNx]);
    end
    localPrintSliceProgress(iz, Nz);
end

finvol = zeros(targetNy, targetNx, targetNz, 1, 'uint16');
if resfac(3) == 1
    finvol(:, :, :, 1) = backvol;
else
    finvol(:, :, :, 1) = imresize3(backvol, targetSize);
end
end

%--------------------------------------------------------------------------
function localAssertPlanePerFileLoadable(nativeSize, sampleres, registrationres, nativeBytes, resfac)
maxBytes = localMaxSafeArrayBytes();
if nativeBytes <= maxBytes && all(abs(resfac - 1) < 1e-6)
    return
end
if all(abs(resfac - 1) < 1e-6)
    error('readSpinalCordSample:VolumeTooLarge', ...
        ['Native stack is %.1f GB (%d x %d x %d px) but sampleres matches registrationres ', ...
        '(no downsampling).\nSet sampleres to your native voxel size in micrometers, e.g.\n', ...
        '  sampleres = [2.0, 2.0, 5.0];  %% xy and z step from your microscope\n', ...
        'Registration still runs at %s um.'], ...
        nativeBytes / 1e9, nativeSize(1), nativeSize(2), nativeSize(3), mat2str(registrationres));
end
[~, targetSize] = cordVolumeDownsampleSpec(nativeSize, sampleres, registrationres);
targetBytes = prod(double(targetSize)) * 2;
if targetBytes > maxBytes
    error('readSpinalCordSample:VolumeTooLarge', ...
        ['Downsampled stack would still be %.1f GB (%s px). Check sampleres (%s um) ', ...
        'and registrationres (%s um).'], ...
        targetBytes / 1e9, mat2str(targetSize), mat2str(sampleres), mat2str(registrationres));
end
end

%--------------------------------------------------------------------------
function maxBytes = localMaxSafeArrayBytes()
maxBytes = 200 * 1024^3;
try
    s = memory;
    if isfield(s, 'MemAvailableAllArrays') && s.MemAvailableAllArrays > 0
        maxBytes = min(maxBytes, s.MemAvailableAllArrays * 0.5);
    end
catch
    % memory() unavailable on some platforms
end
end

%--------------------------------------------------------------------------
function localPrintSliceProgress(iz, Nz)
if mod(iz, 100) == 0 || iz == Nz
    fprintf('  read slice %d / %d\n', iz, Nz);
end
end

%--------------------------------------------------------------------------
function finvol = localLoadSingleTiffStack(tfile)
tiffpath = fullfile(tfile.folder, tfile.name);
dataim   = BioformatsImage(tiffpath);
Nchans   = dataim.sizeC;
[ny, nx, Nz, planesInTime, preferImread] = tiffChannelStackDims(tiffpath);
if Nz == 1 && endsWith(lower(tiffpath), {'.tif', '.tiff'})
    ninfo = numel(imfinfo(tiffpath));
    if ninfo > 1
        Nz = ninfo;
        planesInTime = true;
        preferImread = true;
        fprintf(['Bio-Formats reported 1 Z-plane; imfinfo found %d ' ...
            'TIFF directories — using that as depth.\n'], ninfo);
    end
elseif planesInTime
    fprintf(['Bio-Formats reported depth along T (%d planes); ' ...
        'reading with planes_in_time.\n'], Nz);
end
finvol = zeros(ny, nx, Nz, Nchans, 'uint16');
if preferImread && Nchans == 1
    tinfo = imfinfo(tiffpath);
    for iz = 1:Nz
        finvol(:, :, iz, 1) = imread(tiffpath, 'Index', iz, 'Info', tinfo);
    end
else
    for ichan = 1:Nchans
        for iz = 1:Nz
            if planesInTime
                try
                    plane = dataim.getPlane(1, ichan, iz);
                catch
                    plane = dataim.getPlane(iz);
                end
            else
                plane = dataim.getPlane(iz, ichan, 1, 1);
            end
            finvol(:, :, iz, ichan) = plane;
        end
    end
end
end

%--------------------------------------------------------------------------
function finvol = localLoadChannelPerFileStack(tfiles)
Nchans = numel(tfiles);
finvol = cell(Nchans, 1);
for ichan = 1:Nchans
    finvol{ichan} = localLoadSingleTiffStack(tfiles(ichan));
    if size(finvol{ichan}, 4) ~= 1
        error('readSpinalCordSample:MultiChannelFile', ...
            ['Each TIFF in channelperfile mode must be a single-channel stack. ', ...
            'File %s has %d channels. Use one multi-page TIFF or planeperfile instead.'], ...
            tfiles(ichan).name, size(finvol{ichan}, 4));
    end
    finvol{ichan} = finvol{ichan}(:, :, :, 1);
end
finvol = cat(4, finvol{:});
end
