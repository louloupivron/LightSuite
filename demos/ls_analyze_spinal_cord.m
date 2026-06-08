% set main data path - where the tiffs are
dpspinesample      = 'F:\Louis\spinal-cord\All_Channels';
bcpdpath           = which('bcpd.exe'); % path for point-cloud registration


%% load sample and atlas - set resolution and channel for registration
% sampleres: native voxel size in micrometers — NOT the registration grid (20 um).
sampleres          = [1.8, 1.8, 4]; % in micrometers across sides
% List corrupt slice TIFFs before loading (recommended for large series):
% bad = validateSpinalCordSliceSeries(dpspinesample);
% tifftype: 'auto' | 'planeperfile' (one 2D TIFF per slice) | 'channelperfile'
[cordvol, opts]    = readSpinalCordSample(dpspinesample, sampleres, 'planeperfile');
% If a few slices are corrupt and cannot be re-exported, use:
% [cordvol, opts] = readSpinalCordSample(dpspinesample, sampleres, 'planeperfile', [], 'SkipCorruptSlices', true);
%%
opts.regchan       = 2; % choose registration channel
opts.bcpdpath      = bcpdpath;
regopts            = prepareCordSampleForRegistration(cordvol, opts);

%% (manual) trace spinal cord to untwist and unbend
regopts = loadRegOpts(dpspinesample);
spinal_cord_aligner(regopts);

%% (auto) initialize registration
regopts = loadRegOpts(dpspinesample);
initializeCordRegistration(regopts);

%% (manual) add control points
regopts = loadRegOpts(dpspinesample);
matchControlPointsSpine(regopts);

%% (auto) perform nonlinear registration (b-spline)
regopts = loadRegOpts(dpspinesample);
control_point_weight = 0.2;
tparams = multiobjCordRegistration(regopts, control_point_weight);

%% (auto) register all volume channels to the atlas
transform_params = load(fullfile(regopts.lsfolder, 'transform_params.mat'));
regvol = generateRegisteredCordVolume(regopts, transform_params);