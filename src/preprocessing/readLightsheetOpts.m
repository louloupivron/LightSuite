function opts = readLightsheetOpts(opts)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%--------------------------------------------------------------------------
fprintf('Looking for data in %s\n', opts.datafolder)
%--------------------------------------------------------------------------
tifftype = getOr(opts, 'tifftype', 'planeperfile');
tiffiles = dir(fullfile(opts.datafolder, '*.tiff'));
tifiles  = dir(fullfile(opts.datafolder, '*.tif'));
tfiles   = cat(1, tifiles, tiffiles);
opts.multitiffs = false;
%--------------------------------------------------------------------------
switch tifftype
    %======================================================================
    case 'planeperfile' % classic for Terastitcher
        fprintf('Using planeperfile loading\n')
        %--------------------------------------------------------------------------
        opts.Nz  = numel(tifiles);
        fprintf('Found %d tiff files ', opts.Nz)
        %--------------------------------------------------------------------------
        tinfo    = imfinfo(fullfile(tifiles(1).folder, tifiles(1).name));
        Ny       = tinfo.Height;
        Nx       = tinfo.Width;
        %--------------------------------------------------------------------------
        allnyx = nan(opts.Nz, 2);
        for ii = 1:5:opts.Nz
            tinfo         = imfinfo(fullfile(tifiles(ii).folder, tifiles(ii).name));
            allnyx(ii, :) = [tinfo.Height tinfo.Width];
        end
        icheck = ~isnan(sum(allnyx,2));
        assert(all(allnyx(icheck, :) == [Ny Nx], "all"), ...
            "some slices do not have matching size, LightSuite cannot proceed")
        %--------------------------------------------------------------------------
        opts.Nx     = Nx; 
        opts.Ny     = Ny; 
        opts.Nchans = 1;
        fprintf('of size %d x %d px\n', opts.Ny, opts.Nx)
        %======================================================================
    case 'channelperfile' % classic for BigStitcher
        fprintf('Using channelperfile loading, every channel should be in a different tiff\n')
        Nfiles = numel(tfiles);
        opts.planes_in_time = false;
        opts.use_imread_channelperfile = false;
        tiffpath = fullfile(tfiles(1).folder, tfiles(1).name);
        if Nfiles == 1
             datainfo = BioformatsImage(tiffpath);
             opts.Nchans = datainfo.sizeC;
             fprintf('Found a single tiff with %d channels \n', opts.Nchans);
             sizeZ = datainfo.sizeZ;
             if isprop(datainfo, 'sizeT')
                 sizeT = datainfo.sizeT;
             else
                 sizeT = 1;
             end
             % Multi-page / OME stacks often report depth as T (time) not Z
             if sizeZ == 1 && sizeT > 1
                 opts.Nz = sizeT;
                 opts.planes_in_time = true;
             else
                 opts.Nz = sizeZ;
             end
             % When Bio-Formats still reports a single plane, count TIFF IFDs
             if opts.Nz == 1 && endsWith(lower(tiffpath), {'.tif', '.tiff'})
                 ninfo = numel(imfinfo(tiffpath));
                 if ninfo > 1
                     opts.Nz = ninfo;
                     opts.planes_in_time = true;
                     fprintf(['Bio-Formats reported 1 Z-plane; imfinfo found %d ' ...
                         'TIFF directories — using that as depth.\n'], ninfo);
                 end
             end
             allnyxz = [datainfo.height datainfo.width opts.Nz];
        else
            opts.Nchans  = numel(tfiles);
            fprintf('Assuming each channel is a separate tiff. Found %d channels \n', opts.Nchans)
            allnyxz = nan(opts.Nchans, 3);
            allTiffChannels = true;
            allPreferImread = true;
            for ichan = 1 : opts.Nchans
                currpath = fullfile(tfiles(ichan).folder, tfiles(ichan).name);
                if endsWith(lower(currpath), {'.tif', '.tiff'})
                    [nyC, nxC, nzC, planeT, preferImread] = tiffChannelStackDims(currpath);
                    allnyxz(ichan, :) = [nyC nxC nzC];
                    allPreferImread = allPreferImread && preferImread;
                    if planeT
                        opts.planes_in_time = true;
                    end
                else
                    allTiffChannels = false;
                    allPreferImread = false;
                    datainfo          = BioformatsImage(currpath);
                    allnyxz(ichan, :) = [datainfo.height datainfo.width datainfo.sizeZ];
                end
            end
            opts.multitiffs = true;
            assert(all(allnyxz == allnyxz(1,:), "all"), ...
                "some volumes do not have matching size, LightSuite cannot proceed")
            opts.Nz = allnyxz(1, 3);
            opts.use_imread_channelperfile = allTiffChannels && allPreferImread;
            if opts.use_imread_channelperfile
                fprintf('Using native TIFF page reading for channel-per-file stacks.\n');
            elseif allTiffChannels && ~opts.use_imread_channelperfile && opts.Nz > 1
                fprintf(['Reading channel TIFFs with Bio-Formats (imfinfo did not list ' ...
                    'multiple IFDs; planes_in_time=%d).\n'], opts.planes_in_time);
            end
        end
        %--------------------------------------------------------------------------
        opts.Nx  = allnyxz(1, 2); 
        opts.Ny  = allnyxz(1, 1); 
        fprintf('Volume size is %d x %d x %d px\n', opts.Ny, opts.Nx,  opts.Nz)
        %======================================================================
end
%--------------------------------------------------------------------------
opts.tfiles = fullfile({tfiles(:).folder},{tfiles(:).name})';
makeNewDir(opts.savepath)
%--------------------------------------------------------------------------
% Brain atlas (Allen CCF default; Perens LSFM optional — see resolveBrainAtlasConfig)
if ~isfield(opts, 'brain_atlas') || isempty(opts.brain_atlas)
    opts.brain_atlas = 'allen';
else
    opts.brain_atlas = lower(strtrim(char(opts.brain_atlas)));
end
assert(any(strcmp(opts.brain_atlas, {'allen', 'perens'})), ...
    'opts.brain_atlas must be ''allen'' or ''perens''.');
if ~isfield(opts, 'atlas_dir')
    opts.atlas_dir = [];
end
%--------------------------------------------------------------------------

end

function [ny, nx, nz, planesInTime, preferImread] = tiffChannelStackDims(path)
%TIFFCHANNELSTACKDIMS Height/width/depth for one channel file.
%   preferImread is true when numel(imfinfo)>1 (classic multi-IFD stack).
%   When imfinfo reports a single directory (OME-TIFF, SubIFDs, etc.),
%   use Bio-Formats sizeZ/sizeT instead and set preferImread false.
tinfo  = imfinfo(path);
nPages = numel(tinfo);
ny     = tinfo(1).Height;
nx     = tinfo(1).Width;
if nPages > 1
    nz             = nPages;
    planesInTime   = false;
    preferImread   = true;
    return;
end
preferImread = false;
datainfo = BioformatsImage(path);
ny       = datainfo.height;
nx       = datainfo.width;
sizeZ    = datainfo.sizeZ;
if isprop(datainfo, 'sizeT')
    sizeT = datainfo.sizeT;
else
    sizeT = 1;
end
if sizeZ == 1 && sizeT > 1
    nz           = sizeT;
    planesInTime = true;
else
    nz           = sizeZ;
    planesInTime = false;
end
end