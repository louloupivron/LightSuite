function allmedians = generateRegisteredBrainVolumes(savepath, varargin)
%GENERATEREGISTEREDBRAINVOLUMES Registers brain volumes and calculates intensities.
%
%   ALLMEDIANS = GENERATEREGISTEREDBRAINVOLUMES(SAVEPATH) loads transformation 
%   parameters and options from SAVEPATH, applies the registration to the 
%   volumes, and calculates the background fluorescence in atlas coordinates.
%
%   ALLMEDIANS = GENERATEREGISTEREDBRAINVOLUMES(..., 'writetocsv', VAL) 
%   specifies whether to write the parcellation intensities to a CSV file.
%
%   ALLMEDIANS = GENERATEREGISTEREDBRAINVOLUMES(..., 'saveregisteredvolume', VAL) 
%   specifies whether to save the registered volumes to disk.
%
%   Inputs:
%       savepath             - (char/string) Directory containing the volume 
%                              files, 'transform_params.mat', and 'regopts.mat'.
%
%   Optional Name-Value Parameters:
%       writetocsv           - (logical) If true, writes intensity data to CSV. 
%                              [Defaults to opts.writetocsv or false]
%       saveregisteredvolume - (logical) If true, saves the registered volume 
%                              stack. [Defaults to opts.saveregisteredvol or false]
%
%   Parcellation (loaded from regopts.mat): Allen uses parcellation membership CSV.
%   Perens requires ARA2_annotation_info_avail_regions.csv beside the NIfTIs in
%   opts.atlas_dir; hemisphere split uses opts.perens_ml_axis (default 2) at COM.
%
%   Outputs:
%       allmedians           - (single) Ngroups x 2 x Nchans array containing 
%                              median intensities over brain areas and
%                              hemispheres.
%--------------------------------------------------------------------------
% 1. Load configuration and transforms
%--------------------------------------------------------------------------
trstruct   = load(fullfile(savepath, 'transform_params.mat'));
loadedOpts = load(fullfile(savepath, 'regopts.mat'));
opts       = loadedOpts.opts;
atlas_opts = struct('brain_atlas', getOr(opts, 'brain_atlas', 'allen'), ...
    'atlas_dir', getOr(opts, 'atlas_dir', []));
atlas_cfg  = resolveBrainAtlasConfig(atlas_opts);
%--------------------------------------------------------------------------
% 2. Parse Inputs and establish defaults
%--------------------------------------------------------------------------
% Determine defaults from the loaded opts structure, falling back to false
if isfield(opts, 'writetocsv')
    defaultWriteCsv = opts.writetocsv;
else
    defaultWriteCsv = false;
end

if isfield(opts, 'saveregisteredvol')
    defaultSaveVol  = opts.saveregisteredvol;
else
    defaultSaveVol  = false;
end

% Set up the input parser
p = inputParser;
addRequired(p, 'savepath', @(x) ischar(x) || isstring(x));
addParameter(p, 'writetocsv', defaultWriteCsv, @(x) islogical(x) || isscalar(x));
addParameter(p, 'saveregisteredvolume', defaultSaveVol, @(x) islogical(x) || isscalar(x));

% Parse the arguments passed to the function
parse(p, savepath, varargin{:});

% Assign the final parsed variables for the rest of the script
writetocsv = p.Results.writetocsv;
saveregvol = p.Results.saveregisteredvolume;
%--------------------------------------------------------------------------
% 3. Initialize paths and channel names
%--------------------------------------------------------------------------
registerpath = fullfile(savepath, 'volume_registered');
makeNewDir(registerpath);

channames = cell(opts.Nchans, 1);
if isfield(opts, 'channames')
    channames = opts.channames;
end
%==========================================================================
fprintf('Applying transforms... \n'); savetic = tic;

straightvol = zeros([trstruct.atlassize opts.Nchans], 'uint16');

for ichan = 1:opts.Nchans
    %--------------------------------------------------------------------------
    % load and transform background volume
    volpath    = dir(fullfile(savepath, sprintf('chan_%d_*register*.tif', ichan)));
    currfname  = fullfile(volpath.folder, volpath.name);
    fprintf('Registering %s\n', currfname)
    %--------------------------------------------------------------------------
    backvolume = readDownStack(currfname);
    volume     = permuteBrainVolume(backvolume, trstruct.how_to_perm);

    volumereg  = transformix(volume,trstruct.tform_bspline_samp20um_to_atlas_20um_px,...
        'movingscale', opts.registres*1e-3*[1 1 1]);
    volumereg                   = uint16(abs(volumereg));
    Rmoving                     = imref3d(size(volumereg));
    Rfixed                      = imref3d(trstruct.atlassize);
    straightvol(:, :, :, ichan) = imwarp(volumereg, Rmoving, trstruct.tform_affine_samp20um_to_atlas_10um_px,...
        'OutputView',Rfixed);
    %--------------------------------------------------------------------------
    fprintf('Channel %d/%d done. Time %2.2f s. \n', ichan, opts.Nchans, toc(savetic));
end
%==========================================================================
if saveregvol
    fprintf('Saving registered volumes... '); savetic = tic;
    saveLargeSliceVolume(permute(straightvol, [1 2 4 3]), channames, registerpath);
    fprintf('Done! Took %2.2f s. \n', toc(savetic));
end
%==========================================================================
fprintf('Calculating background fluoresence in atlas coords...\n'); proctic = tic;

if strcmp(atlas_cfg.brain_atlas, 'allen') && atlas_cfg.supports_parcellation
    av                      = niftiread(atlas_cfg.annotation_path);
    allen_atlas_parcel_path = fileparts(which('parcellation_to_parcellation_term_membership.csv'));

    parcelinfo       = readtable(fullfile(allen_atlas_parcel_path, 'parcellation_to_parcellation_term_membership.csv'));
    substridx        = strcmp(parcelinfo.parcellation_term_set_name, 'substructure');
    [areaidx, ib]    = unique(parcelinfo.parcellation_index(substridx));
    namessub         = parcelinfo.parcellation_term_name(substridx);
    stridx           = strcmp(parcelinfo.parcellation_term_set_name, 'structure');
    [~, ibstr]       = unique(parcelinfo.parcellation_index(stridx));
    namesstruct      = parcelinfo.parcellation_term_name(stridx);
    dividx           = strcmp(parcelinfo.parcellation_term_set_name, 'division');
    [~, ibdiv]       = unique(parcelinfo.parcellation_index(dividx));
    namesdiv         = parcelinfo.parcellation_term_name(dividx);

    Ngroups          = numel(areaidx);
    Nforaccum        = max(av, [], 'all') + 1;
    Npxlr            = size(av,3)/2;

    allmedians = nan(Ngroups, 2, opts.Nchans, 'single');
    for ichan = 1:opts.Nchans

        medianoverareas  = nan(Ngroups, 2, 'single');
        volumeoverareas  = nan(Ngroups, 2, 'single');
        for iside = 1:2
            istart = (iside - 1) * Npxlr + 1;
            iend   = istart + Npxlr - 1;
        
            sideav    = reshape(av(:, :, istart:iend), [], 1);
            sidevals  = reshape(straightvol(:, :, istart:iend, ichan), [], 1);
            ikeep     = sideav>0;
            medareas  = single(accumarray(sideav(ikeep)+1, sidevals(ikeep), [Nforaccum 1], @median));
            medareas  = medareas(areaidx+1);

            volareas  = accumarray(sideav+1, 1, [Nforaccum 1], @sum);
            volareas  = volareas(areaidx+1);
            volumeoverareas(:, iside) = volareas * (opts.atlasres*1e-3)^3;

            % get index 0 level for background
            backlevel                 = single(median(sidevals(~ikeep)));
            medareas(areaidx == 0)    = backlevel;
            medianoverareas(:, iside) = medareas;
        end

        allmedians(:, :, ichan) = medianoverareas;

        % save as mat file for later processing
        fmatname  = fullfile(registerpath, sprintf('chan%02d_intensities.mat', ichan));
        save(fmatname, 'medianoverareas', 'areaidx', 'volumeoverareas')

        % save as csv if asked
        if writetocsv
            currtable = array2table([areaidx medianoverareas volumeoverareas], ...
                'VariableNames',...
                {'parcellation_index', ...
                'RightSideIntensity', 'LeftSideIntensity',...
                 'RightSideVolume[mm3]', 'LeftSideVolume[mm3]'});
            currtable = addvars(currtable, namessub(ib), namesstruct(ibstr),  namesdiv(ibdiv),...
                'NewVariableNames',{'name', 'structure', 'division'},'Before','parcellation_index');
            fsavename      = fullfile(registerpath, sprintf('chan%02d_intensities.csv', ichan));
            writetable(currtable, fsavename)
        end
        %--------------------------------------------------------------------------
        fprintf('Channel %d/%d done. Time %2.2f s. \n', ichan, opts.Nchans, toc(proctic));
        %--------------------------------------------------------------------------
    end

elseif strcmp(atlas_cfg.brain_atlas, 'perens') && atlas_cfg.supports_parcellation
    % Perens LSFM: ARA-style label IDs (sparse large integers); hemisphere split by ML plane
    % through COM (opts.perens_ml_axis), not Allen-style Z-slabs (see resolveBrainAtlasConfig).
    av = niftiread(atlas_cfg.annotation_path);
    assert(isequal(size(av), trstruct.atlassize), 'LightSuite:PerensParcel', ...
        'Annotation volume size %s must match transform_params.atlassize %s.', ...
        mat2str(size(av)), mat2str(trstruct.atlassize));

    stTab = readtable(atlas_cfg.structures_csv_path);
    areaidx = stTab.id(:);
    Ngroups = numel(areaidx);

    ml_axis = double(getOr(opts, 'perens_ml_axis', 2));
    assert(ml_axis >= 1 && ml_axis <= 3 && ml_axis == floor(ml_axis), ...
        'LightSuite:PerensParcel', 'opts.perens_ml_axis must be 1, 2, or 3 (got %g).', ml_axis);

    sz = size(av);
    brain = av > 0;
    switch ml_axis
        case 1
            coord = repmat(single((1:sz(1))'), [1, sz(2), sz(3)]);
        case 2
            coord = repmat(single(1:sz(2)), [sz(1), 1, sz(3)]);
        otherwise
            coord = repmat(reshape(single(1:sz(3)), [1, 1, sz(3)]), [sz(1), sz(2), 1]);
    end
    split_plane = round(mean(coord(brain)));
    sideLower = coord <= split_plane & brain;
    sideUpper = coord > split_plane & brain;

    fprintf(['Perens parcellation: ML split along axis %d at plane %d (lower idx = Right slot, ' ...
        'higher = Left slot; compare Allen orientation before Brainglobe export).\n'], ...
        ml_axis, split_plane);

    parent_acronym = cell(Ngroups, 1);
    id_col = double(stTab.id(:));
    pid_col = stTab.parent_id(:);
    acr_raw = stTab.acronym;

    for g = 1:Ngroups
        pid = pid_col(g);
        if isempty(pid) || (isnumeric(pid) && isnan(double(pid)))
            parent_acronym{g} = '';
            continue
        end
        ix = find(id_col == double(pid), 1);
        if isempty(ix)
            parent_acronym{g} = '';
        else
            parent_acronym{g} = char(acr_raw(ix));
        end
    end

    acronyms_cell = cellstr(stTab.acronym);
    structure_names = cellstr(stTab.name);

    voxel_mm3 = (opts.atlasres * 1e-3)^3;
    allmedians = nan(Ngroups, 2, opts.Nchans, 'single');

    for ichan = 1:opts.Nchans
        medianoverareas = nan(Ngroups, 2, 'single');
        volumeoverareas = nan(Ngroups, 2, 'single');
        sv = straightvol(:, :, :, ichan);

        for iside = 1:2
            if iside == 1
                sm = sideLower;
            else
                sm = sideUpper;
            end
            lab = av(sm);
            vals = double(sv(sm));
            mask = lab > 0;
            lab = lab(mask);
            vals = vals(mask);
            if isempty(lab)
                continue
            end
            [G, uids] = findgroups(lab);
            meds = splitapply(@median, vals, G);
            nvox_per = splitapply(@numel, vals, G);
            for k = 1:numel(uids)
                row = find(areaidx == uids(k), 1);
                if isempty(row)
                    continue
                end
                medianoverareas(row, iside) = single(meds(k));
                volumeoverareas(row, iside) = nvox_per(k) * voxel_mm3;
            end
        end

        allmedians(:, :, ichan) = medianoverareas;

        fmatname = fullfile(registerpath, sprintf('chan%02d_intensities.mat', ichan));
        save(fmatname, 'medianoverareas', 'areaidx', 'volumeoverareas')

        if writetocsv
            currtable = array2table([double(areaidx), medianoverareas, volumeoverareas], ...
                'VariableNames', ...
                {'parcellation_index', ...
                'RightSideIntensity', 'LeftSideIntensity', ...
                'RightSideVolume[mm3]', 'LeftSideVolume[mm3]'});
            currtable = addvars(currtable, acronyms_cell, structure_names, ...
                parent_acronym, ...
                'NewVariableNames', {'name', 'structure', 'division'}, 'Before', 'parcellation_index');
            fsavename = fullfile(registerpath, sprintf('chan%02d_intensities.csv', ichan));
            writetable(currtable, fsavename)
        end
        fprintf('Channel %d/%d done (Perens). Time %2.2f s. \n', ichan, opts.Nchans, toc(proctic));
    end

elseif strcmp(atlas_cfg.brain_atlas, 'perens') && isempty(atlas_cfg.structures_csv_path)
    fprintf(['Skipping Perens parcellation: ARA2_annotation_info_avail_regions.csv not found ' ...
        'in atlas dir %s\n'], atlas_cfg.atlas_dir);
    allmedians = [];
else
    fprintf(['Skipping intensity parcellation (brain_atlas=%s). ' ...
        'Registered volumes were still saved if requested.\n'], atlas_cfg.brain_atlas);
    allmedians = [];
end
%==========================================================================
end