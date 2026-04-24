% Make sure which invertElastixTransformCP -all returns 
% C:\Users\ALICE_lightsheet\Documents\LightSuite-feature-multi-channel-registration\src\helpers\invertElastixTransformCP.m
% otherwise add to path 
opts = struct();
%=========================================================================
% options to change
%--------------------------------------------------------------------------
% for naming
opts.mousename  = 'Giulia';
% change for the folder that contains stitched tiff files
opts.datafolder = 'F:\Louis\mesoSPIM\Fanny\2515\LLS';
opts.fproc      = fullfile('D:\DATA_LightSuite_Temp'); % where the processed volume is saved as a binary (fast SSD, at
% least 500 GB), will be deleted
parallel.gpu.enableCUDAForwardCompatibility(true)
% path to save results
opts.savepath   = fullfile(opts.datafolder, 'test');
%--------------------------------------------------------------------------
% some processing options
opts.tifftype           = 'channelperfile'; % can be planeperfile or channelperfile
opts.pxsize             = [5.26 5.26 5]; % voxel size, xy and z, in um
opts.atlasres           = 10; % atlas isotropic voxel size in um (20 for perens, 25 for perens2023)
opts.registres          = 20; % resolution to do the nonrigid registration, keep fixed, in um
% opts.brain_atlas      = 'perens';    % optional: Perens 2020 Gubra LSFM (gubra_template_olf.nii.gz + gubra_ano_olf.nii.gz)
% opts.brain_atlas      = 'perens2023'; % optional: Perens 2023 multimodal LSFM only (lsfm_temp.nii.gz + lsfm_ano.nii.gz from LSFM_space_oriented/)
% opts.atlas_dir        = '';      % optional: explicit atlas folder; if empty, which() finds the template
% cell detection parameters
opts.debug              = true; % toggle plotting (takes longer) for cell detections
opts.savecellimages     = false; % toggle saving of individual cell images
opts.celldiam           = 14; % approximate cell size in um
opts.thres_cell_detect  = [0.5 0.4]; % thresholds for detecting cells relative to background, first should be bigger than second
opts.channelforcells    = []; % channel to use for cell detection, leave empty ([]) for none
opts.writetocsv         = false; % write cells to csv file
%  registration
opts.channelforregister = 1; % channel to use for registration (autofluorescence / structural)
opts.channelforregister_secondary = 2; % Optional: second channel for dual-AMI B-spline (same Allen template, two sample contrasts)
%--------------------------------------------------------------------------
opts                   = readLightsheetOpts(opts);
%=========================================================================
%% (auto) main processing pipeline, preprocess and detect cells
opts = preprocessLightSheetVolume(opts);

%% (auto) initialize registration
opts = initializeRegistration(opts.savepath);

%% (manual) Refine control points
optsfile = dir(fullfile(opts.savepath, 'regopts.mat'));
if ~isempty(optsfile)
    opts = load(fullfile(optsfile.folder, optsfile.name));
    opts = opts.opts;
end
matchControlPoints_unified(opts);

%% (auto) perform registration
optsfile = dir(fullfile(opts.savepath, 'regopts.mat'));
if ~isempty(optsfile)
    opts = load(fullfile(optsfile.folder, optsfile.name));
    opts = opts.opts;
end
% options for registration
opts.augmentpoints         = false; % whether to augment user-defined points with automatic ones
opts.weight_usr_pts        = 0.2;   % weight of user-defined points for atlas fitting, set to zero for image-only information
opts.dual_channel_mi_weight_autofluor = 0.4;
opts.dual_channel_mi_weight_signal    = 0.4;
opts.bspline_spatial_scale = 0.64;  % in mm, how much you allow the bspline to bend (smaller is higher, but more prone to noise)
transform_params           = multiobjRegistration(opts, opts.weight_usr_pts, true);

%% (auto) apply registration to volume
% change function arguments if you want to skip saving
generateRegisteredBrainVolumes(opts.savepath, ...
    'writetocsv', true, 'saveregisteredvolume', false);

%% (auto) apply registration to cell detections
% you can also run this function with a cell matrix of detections
% calculated outside LightSuite, as long as the transform parameters are
% provided as an argument
transformPointsToAtlas(opts.savepath, 'writetocsv', true);
