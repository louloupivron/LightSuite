function [regimg,tform_bspline, tformpath, pathtemp] = performMultObjBsplineRegistration(movingvol,fixedvol,...
    volscale, movingpts, fixedpts, savepath, optsreg, varargin)
%UNTITLED3 B-spline registration with elastix (single or dual sample-channel AMI).
%
%   Optional name-value pair:
%     'FixedSampleSecondary' — second fixed image (e.g. fluorescent signal at the same
%       grid as fixedvol). When provided, elastix uses two AdvancedMattesMutualInformation
%       terms (shared moving atlas volume) plus the landmark metric, with weights
%       dual_channel_mi_weight_autofluor, dual_channel_mi_weight_signal, and cpwt.
%
%   See also: multiobjRegistration, elastixDualFixedSameMovingBspline

%==========================================================================
ip = inputParser;
ip.addParameter('FixedSampleSecondary', [], @(x) isempty(x) || isnumeric(x));
ip.parse(varargin{:});
fixedSecondary = ip.Results.FixedSampleSecondary;
%==========================================================================
cpwt         = getOr(optsreg, 'cpwt', 0.3);
nhistbins    = getOr(optsreg, 'n_histogram_bins', 88);

if isempty(movingpts) || isempty(fixedpts)
    cpwt = 0;
end

dualMi = ~isempty(fixedSecondary);
if dualMi
    assert(isequal(size(fixedSecondary), size(fixedvol)), ...
        'LightSuite:performMultObjBsplineRegistration:sizeMismatch', ...
        'FixedSampleSecondary must match fixedvol size [%s], got [%s].', ...
        mat2str(size(fixedvol)), mat2str(size(fixedSecondary)));
end

params = localBuildBsplineParams(optsreg, dualMi, cpwt, nhistbins, fixedvol, volscale);

%==========================================================================
pathtemp = fullfile(savepath, 'elastix_temp');
makeNewDir(pathtemp);
clearElastixWorkspaceForNewRun(pathtemp);

movpath  = fullfile(pathtemp, 'moving.txt');
fixpath  = fullfile(pathtemp, 'fixed.txt');

% -1 because offset is always zero. this way, the first pixel is at 0
writePointsFile(movpath, (movingpts-1)*volscale)
writePointsFile(fixpath, (fixedpts-1)*volscale)

addElastixRepoPaths;
fprintf('Performing B-spline registration with elastix...\n'); tic;
if dualMi
    fprintf(['Dual-channel MI: atlas moving vs fixed AF + fixed signal ', ...
        '(weights AF=%.3g signal=%.3g landmarks=%.3g).\n'], ...
        getOr(optsreg, 'dual_channel_mi_weight_autofluor', 1), ...
        getOr(optsreg, 'dual_channel_mi_weight_signal', 0.5), ...
        cpwt);
    [regimg,tform_bspline] = elastixDualFixedSameMovingBspline(movingvol, fixedvol, fixedSecondary, ...
        pathtemp, volscale, params, movpath, fixpath);
else
    [regimg,tform_bspline] = elastix(movingvol, fixedvol, pathtemp,'elastix_default.yml','paramstruct',params,...
        'movingpoints', movpath,  'fixedpoints', fixpath,...
        'movingscale', volscale*[1 1 1],  'fixedscale', volscale*[1 1 1]);
end
fprintf('Done! Took %2.2f s.\n', toc)

% copy file and delete temporary folder
tformpath = fullfile(savepath, 'bspline_atlas_to_samp_20um.txt');
copyfile(tform_bspline.TransformParametersFname{1}, tformpath)
% rmdir(pathtemp, 's');
%==========================================================================
end

%--------------------------------------------------------------------------
function params = localBuildBsplineParams(optsreg, dualMi, cpwt, nhistbins, fixedvol, volscale)

usemultistep = getOr(optsreg, 'usemultistep', true);
bspscale     = getOr(optsreg, 'bspline_spatial_scale', 0.64);

params = struct();
params.Registration                  = 'MultiMetricMultiResolutionRegistration';
if dualMi
    w_af = getOr(optsreg, 'dual_channel_mi_weight_autofluor', 0.5);
    w_sg = getOr(optsreg, 'dual_channel_mi_weight_signal', 1.0);
    params.Metric                      = {'AdvancedMattesMutualInformation', ...
        'AdvancedMattesMutualInformation', 'CorrespondingPointsEuclideanDistanceMetric'};
    % Elastix 5.1 ITK: counts of pyramids / samplers / interpolators must be 1 or equal
    % NumberOfMetrics (3). Third image pair on disk duplicates AF+atlas (see
    % elastixDualFixedSameMovingBspline) so pyramid 2 has valid inputs for the
    % landmark metric slot.
    params.FixedImagePyramid           = {'FixedRecursiveImagePyramid', 'FixedRecursiveImagePyramid', ...
        'FixedRecursiveImagePyramid'};
    params.MovingImagePyramid          = {'MovingRecursiveImagePyramid', 'MovingRecursiveImagePyramid', ...
        'MovingRecursiveImagePyramid'};
    params.ImageSampler                = {'RandomCoordinate', 'RandomCoordinate', 'RandomCoordinate'};
    params.Interpolator                = {'BSplineInterpolator', 'BSplineInterpolator', 'BSplineInterpolator'};
    params.Metric0Weight               = w_af;
    params.Metric1Weight               = w_sg;
    params.Metric2Weight               = cpwt;
else
    params.Metric                      = {'AdvancedMattesMutualInformation', ...
        'CorrespondingPointsEuclideanDistanceMetric'};
    params.FixedImagePyramid           = 'FixedRecursiveImagePyramid';
    params.MovingImagePyramid          = 'MovingRecursiveImagePyramid';
    params.ImageSampler                = 'RandomCoordinate';
    params.Metric1Weight                = cpwt;
    params.Metric0Weight               = 1.0;
end

params.Transform                       = 'RecursiveBSplineTransform';
params.Optimizer                       = 'AdaptiveStochasticGradientDescent';
params.AutomaticParameterEstimation    = true;
params.AutomaticScalesEstimation       = false;
params.BSplineInterpolationOrder       = 3;
params.FinalBSplineInterpolationOrder  = 3;
params.FixedImageDimension             = 3;
params.MovingImageDimension            = 3;
params.UseRandomSampleRegion           = true;
params.NewSamplesEveryIteration        = true;
params.NumberOfResolutions             = 4;
params.NumberOfHistogramBins           = nhistbins;
params.SP_A                            = 20;
params.MaximumNumberOfIterations       = [500  1000 1500 2000];
params.NumberOfSpatialSamples          = 5000;
params.ImagePyramidSchedule            = [8*ones(1,3) 4*ones(1,3) 2*ones(1,3) 1*ones(1,3)];
params.FinalGridSpacingInPhysicalUnits   = bspscale*ones(1,3);

fixed_size_mm     = size(fixedvol) .* volscale;
max_sample_region = min(fixed_size_mm) / 3;
if usemultistep
    base_sizes = [4.5 4 3 2];
else
    base_sizes = [2];
end
base_sizes = min(base_sizes, max_sample_region);
if usemultistep
    params.SampleRegionSize = [base_sizes(1)*ones(1,3) base_sizes(2)*ones(1,3) base_sizes(3)*ones(1,3) base_sizes(4)*ones(1,3)];
else
    params.SampleRegionSize = base_sizes(1)*ones(1,3);
end
end
