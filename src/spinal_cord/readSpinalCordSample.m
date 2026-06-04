function [finvol, opts] = readSpinalCordSample(dp, sampleres)
%READSPINALCORDSAMPLE Load a spinal cord lightsheet TIFF stack.
%   [finvol, opts] = readSpinalCordSample(dp, sampleres)
%   dp is a folder with one multi-page TIFF or one TIFF per channel.
%   finvol is [height x width x depth x channels] in uint16.

tiffiles = dir(fullfile(dp, '*.tiff'));
tifiles  = dir(fullfile(dp, '*.tif'));
tfiles   = cat(1, tifiles, tiffiles);
if isempty(tfiles)
    error('readSpinalCordSample:NoTiffs', ...
        'No .tif or .tiff files found in:\n  %s', dp);
end
tic;
if numel(tfiles) == 1
    tiffpath = fullfile(tfiles(1).folder, tfiles(1).name);
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
else
    Nchans = numel(tfiles);
    finvol = cell(Nchans, 1);
    for ichan = 1 : Nchans
        tiffpath = fullfile(tfiles(ichan).folder, tfiles(ichan).name);
        [ny, nx, Nz, planesInTime, preferImread] = tiffChannelStackDims(tiffpath);
        dataim = BioformatsImage(tiffpath);
        currvol = zeros(ny, nx, Nz, 'uint16');
        if preferImread
            tinfo = imfinfo(tiffpath);
            for iz = 1:Nz
                currvol(:, :, iz) = imread(tiffpath, 'Index', iz, 'Info', tinfo);
            end
        else
            for iz = 1:Nz
                if planesInTime
                    try
                        plane = dataim.getPlane(1, 1, iz);
                    catch
                        plane = dataim.getPlane(iz);
                    end
                else
                    plane = dataim.getPlane(iz, 1, 1, 1);
                end
                currvol(:, :, iz) = plane;
            end
        end
        finvol{ichan} = currvol;
    end
    finvol = cat(4, finvol{:});
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
%--------------------------------------------------------------------------
targetres = [20 20 20];
%--------------------------------------------------------------------------
opts.datafolder     = tfiles(1).folder;
opts.lsfolder       = fullfile(opts.datafolder, 'lightsuite');
makeNewDir(opts.lsfolder)
opts.orisize         = [Nslices, Ny, Nx];
opts.Nchan           = Nchan;
opts.sampleres       = sampleres;
opts.registrationres = targetres;
%--------------------------------------------------------------------------

end
