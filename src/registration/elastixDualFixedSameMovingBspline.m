function [regimg, tform_bspline] = elastixDualFixedSameMovingBspline(movingVol, fixedVol1, fixedVol2, pathtemp, volscale, params, movpath, fixpath)
%   VOLSCALE is voxel spacing in mm (scalar or [sx sy sz]) passed to mhd_write as
%   ElementSpacing for all four volumes (must match single-channel registration).
%ELASTIXDUALFIXEDSAMEMOVINGBSPLINE Run elastix with dual sample contrasts + landmarks.
%
%   Elastix is invoked with three fixed/moving pairs (-f0/-m0 .. -f2/-m2) because
%   MultiMetricMultiResolutionRegistration (Elastix 5.1) requires the number of
%   fixed/moving pyramids to be 1 or equal to the number of metrics (3 here).
%   Pairs 0 and 1 are AF+atlas and signal+atlas; pair 2 duplicates pair 0 so the
%   landmark metric's pyramid chain receives valid images (data identical to f0/m0).
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
baseF2 = fullfile(root, [dirName '_dual_f2']);
baseM2 = fullfile(root, [dirName '_dual_m2']);

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
% Third pair duplicates f0/m0 (same arrays) so ITK has three image inputs for three metrics.
mhd_write(fixedVol1, baseF2, sp);
mhd_write(movingVol, baseM2, sp);

paramFname = fullfile(root, sprintf('%s_parameters_dual.txt', dirName));
elastix_parameter_write(paramFname, 'elastix_default.yml', params);
% matlab_elastix YAML validation omits Metric2Weight (not in default.yml); append if missing.
if isfield(params, 'Metric2Weight')
    localAppendElastixParamIfMissing(paramFname, 'Metric2Weight', params.Metric2Weight);
end

ptsCmd = '';
if nargin >= 8 && ~isempty(movpath) && ~isempty(fixpath) && ...
        exist(movpath, 'file') == 2 && exist(fixpath, 'file') == 2
    ptsCmd = sprintf(' -fp "%s" -mp "%s" ', fixpath, movpath);
end

cmd = sprintf(['elastix -f0 "%s.mhd" -m0 "%s.mhd" -f1 "%s.mhd" -m1 "%s.mhd" ' ...
    '-f2 "%s.mhd" -m2 "%s.mhd" -out "%s" %s-p "%s" '], ...
    baseF0, baseM0, baseF1, baseM1, baseF2, baseM2, root, ptsCmd, paramFname);

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

function localAppendElastixParamIfMissing(fname, key, val)
% Append (Key value) line if the parameter file does not already define Key.
if exist(fname, 'file') ~= 2
    return
end
txt = fileread(fname);
token = sprintf('(%s ', key);
if contains(txt, token)
    return
end
fid = fopen(fname, 'a');
if fid < 0
    return
end
fprintf(fid, '\n(%s %g)\n', key, val);
fclose(fid);
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
