function [regimg, tform_bspline] = elastixDualFixedSameMovingBspline(movingVol, fixedVol1, fixedVol2, pathtemp, volscale, params, movpath, fixpath)
%   VOLSCALE is voxel spacing in mm (scalar or [sx sy sz]) passed to mhd_write as
%   ElementSpacing for all four volumes (must match single-channel registration).
%ELASTIXDUALFIXEDSAMEMOVINGBSPLINE Run elastix with two fixed images and one moving image.
%
%   Elastix is invoked as:
%     elastix -f0 fixedAF.mhd -m0 movingAtlas.mhd -f1 fixedSignal.mhd -m1 movingAtlas.mhd ...
%
%   The same moving volume is written twice (m0, m1) so each
%   AdvancedMattesMutualInformation term pairs the warped Allen template with
%   the autofluorescence and signal downsamples respectively. One shared
%   B-spline transform is optimised (MultiMetricMultiResolutionRegistration).
%
%   Requires matlab_elastix on the path (addElastixRepoPaths), elastix binary
%   on the system PATH, and mhd_write / elastix_parameter_write from that package.
%
%   See also: performMultObjBsplineRegistration, multiobjRegistration

addElastixRepoPaths;

[~, dirName] = fileparts(pathtemp);
root = pathtemp;

baseF0 = fullfile(root, [dirName '_dual_f0']);
baseM0 = fullfile(root, [dirName '_dual_m0']);
baseF1 = fullfile(root, [dirName '_dual_f1']);
baseM1 = fullfile(root, [dirName '_dual_m1']);

sp = volscale(:).';
if isscalar(sp)
    sp = sp * [1, 1, 1];
elseif numel(sp) ~= 3
    error('elastixDualFixedSameMovingBspline:BadScale', ...
        'volscale must be scalar or length-3, got %d elements.', numel(sp));
end

mhd_write(fixedVol1, baseF0, sp);
mhd_write(movingVol, baseM0, sp);
mhd_write(fixedVol2, baseF1, sp);
mhd_write(movingVol, baseM1, sp);

paramFname = fullfile(root, sprintf('%s_parameters_dual.txt', dirName));
elastix_parameter_write(paramFname, 'elastix_default.yml', params);

ptsCmd = '';
if nargin >= 8 && ~isempty(movpath) && ~isempty(fixpath) && ...
        exist(movpath, 'file') == 2 && exist(fixpath, 'file') == 2
    ptsCmd = sprintf(' -fp "%s" -mp "%s" ', fixpath, movpath);
end

cmd = sprintf(['elastix -f0 "%s.mhd" -m0 "%s.mhd" -f1 "%s.mhd" -m1 "%s.mhd" -out "%s" ' ...
    '%s-p "%s" '], ...
    baseF0, baseM0, baseF1, baseM1, root, ptsCmd, paramFname);

cmdFid = fopen(fullfile(root, 'CMD_dual'), 'w');
fprintf(cmdFid, '%s\n', cmd);
fclose(cmdFid);

fprintf('Running dual-fixed elastix:\n%s\n', cmd);
[status, result] = system(cmd);
if status ~= 0
    error('elastixDualFixedSameMovingBspline:ElastixFailed', ...
        'Dual-channel elastix failed (exit %d).\n%s', status, result);
end

d = dir(fullfile(root, 'result*.*'));
d(cellfun(@(x) endsWith(x, '.raw'), {d.name})) = [];
if isempty(d)
    error('elastixDualFixedSameMovingBspline:NoResult', ...
        'No result image found in %s.\n%s', root, result);
end
fullPath = fullfile(root, d(end).name);
regimg = localLoadElastixResult(fullPath);

dT = dir(fullfile(root, 'TransformParameters.*.txt'));
if isempty(dT)
    error('elastixDualFixedSameMovingBspline:NoTform', 'No TransformParameters in %s', root);
end
tform_bspline = struct();
tform_bspline.TransformParametersFname = {fullfile(root, dT(1).name)};
tform_bspline.TransformParameters = {elastix_parameter_read(tform_bspline.TransformParametersFname{1})};

scaleStr = sprintf('%.12g ', volscale(:));
fprintf('Done dual-channel B-spline. Voxel spacing (mm) used in MHD: %s\n', strtrim(scaleStr));
end

function im = localLoadElastixResult(fname)
[~, ~, ext] = fileparts(fname);
if strcmpi(ext, '.mhd')
    im = mhd_read(fname);
elseif strcmpi(ext, '.tif') || strcmpi(ext, '.tiff')
    im = load3Dtiff(fname);
else
    error('elastixDualFixedSameMovingBspline:UnknownFormat', 'Unknown result extension: %s', ext);
end
end
