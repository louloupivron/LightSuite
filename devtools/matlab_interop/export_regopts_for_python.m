function export_regopts_for_python(savepath)
%EXPORT_REGOPTS_FOR_PYTHON  Save numeric transforms for Python parity imports.
%
%   export_regopts_for_python(savepath)
%
%   Writes regopts_python_export.mat with 4x4 double matrices that scipy.io.loadmat
%   can read (MATLAB affinetform3d objects are not readable in Python).
%
%   Run once after init-registration (and again after register for the final affine):
%     export_regopts_for_python('F:\...\register_Yosi_test');

if nargin < 1 || isempty(savepath)
    savepath = pwd;
end
savepath = char(savepath);

S = load(fullfile(savepath, 'regopts.mat'));
opts = S.opts;

if ~isfield(opts, 'original_trans')
    error('export_regopts_for_python:MissingTransform', ...
        'regopts.mat has no original_trans — run initializeRegistration first.');
end

original_trans_matrix = affinetform_to_matrix(opts.original_trans);

tform_affine_samp20um_to_atlas_10um_px_matrix = [];
tfile = fullfile(savepath, 'transform_params.mat');
if isfile(tfile)
    tp = load(tfile);
    if isfield(tp, 'tform_affine_samp20um_to_atlas_10um_px')
        tform_affine_samp20um_to_atlas_10um_px_matrix = ...
            affinetform_to_matrix(tp.tform_affine_samp20um_to_atlas_10um_px);
    end
end

out = fullfile(savepath, 'regopts_python_export.mat');
if isempty(tform_affine_samp20um_to_atlas_10um_px_matrix)
    save(out, 'original_trans_matrix', '-v7');
else
    save(out, 'original_trans_matrix', 'tform_affine_samp20um_to_atlas_10um_px_matrix', '-v7');
end
fprintf('Wrote %s\n', out);
end

%--------------------------------------------------------------------------
function M = affinetform_to_matrix(tform)
% Row-vector convention used by LightSuite Python (same as affinetform3d.T).
%
% Python ``transform_points`` applies ``[x y z 1] @ M`` via ``hom @ M.T``,
% which matches MATLAB ``transformPointsForward`` with row-vector points.
if isnumeric(tform) && isequal(size(tform), [4 4])
    M = double(tform);
    return;
end
if ~isa(tform, 'affinetform3d')
    error('Expected affinetform3d or 4x4 numeric, got %s', class(tform));
end
if isprop(tform, 'T')
    % MATLAB row-vector layout: [x y z 1] * T  (translation in T(4,1:3)).
    % Python stores the premultiply transpose: hom @ M.T with translation in M(1:3,4).
    M = double(tform.T.');
    return;
end
% Legacy fallback (older Image Processing Toolbox builds).
A = double(tform.A);
if isprop(tform, 'Translation')
    t = double(tform.Translation(:))';
elseif isprop(tform, 't')
    t = double(tform.t(:))';
else
    error('affinetform3d has no T, Translation, or t property.');
end
M = [A, t'; 0 0 0 1];
end
