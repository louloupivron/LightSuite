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
        opts.use_tiffreadvolume_channelperfile = false;
        tiffpath = fullfile(tfiles(1).folder, tfiles(1).name);
        if Nfiles == 1
             datainfo = BioformatsImage(tiffpath);
             opts.Nchans = datainfo.sizeC;
             fprintf('Found a single tiff with %d channels \n', opts.Nchans);
             allnyxz = [datainfo.height datainfo.width datainfo.sizeZ];
             opts.Nz = datainfo.sizeZ;
        else
            opts.Nchans  = numel(tfiles);
            fprintf('Assuming each channel is a separate tiff. Found %d channels \n', opts.Nchans)
            allnyxz = nan(opts.Nchans, 3);
            for ichan = 1 : opts.Nchans
                currpath = fullfile(tfiles(ichan).folder, tfiles(ichan).name);
                if endsWith(lower(currpath), {'.tif', '.tiff'})
                    tinfo = imfinfo(currpath);
                    Nzcurr = numel(tinfo);
                    if Nzcurr == 1 && isfield(tinfo(1), 'ImageDescription')
                        desc = tinfo(1).ImageDescription;
                        tok = regexp(desc, '(images|slices)\s*=\s*(\d+)', 'tokens', 'once');
                        if ~isempty(tok)
                            Nzcurr = str2double(tok{2});
                            if Nzcurr > 1
                                opts.use_tiffreadvolume_channelperfile = true;
                            end
                        end
                    end
                    allnyxz(ichan, :) = [tinfo(1).Height tinfo(1).Width Nzcurr];
                    if Nzcurr > 1 && ~opts.use_tiffreadvolume_channelperfile
                        opts.use_imread_channelperfile = true;
                    end
                else
                    datainfo          = BioformatsImage(currpath);
                    allnyxz(ichan, :) = [datainfo.height datainfo.width datainfo.sizeZ];
                end
            end
            opts.multitiffs = true;
            assert(all(allnyxz == allnyxz(1,:), "all"), ...
                "some volumes do not have matching size, LightSuite cannot proceed")
            opts.Nz = allnyxz(1, 3);
            if opts.use_imread_channelperfile
                fprintf('Using native TIFF page reading for channel-per-file stacks.\n');
            elseif opts.use_tiffreadvolume_channelperfile
                fprintf('Using tiffreadVolume for channel-per-file TIFF stacks.\n');
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

end