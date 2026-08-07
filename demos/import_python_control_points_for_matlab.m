function import_python_control_points_for_matlab(savepath)
%IMPORT_PYTHON_CONTROL_POINTS_FOR_MATLAB  Verify Python-exported control points.
%
%   import_python_control_points_for_matlab(savepath)
%
%   After running on the Python side:
%     lightsuite brain export-matlab-control-points --save-path <savepath>
%
%   Copy the generated atlas2histology_tform.mat into ``savepath`` (or run
%   export with -o pointing here), then call this function before register.
%
%   It prints the same matched-pair counts and coordinate ranges that
%   multiobjRegistration.m will use (including the [2 1 3] axis swap and
%   original_trans inverse on sample points).

if nargin < 1 || isempty(savepath)
    savepath = pwd;
end
savepath = char(savepath);

cpfile = fullfile(savepath, 'atlas2histology_tform.mat');
if ~isfile(cpfile)
    error('import_python_control_points_for_matlab:MissingFile', ...
        'Expected %s — run lightsuite brain export-matlab-control-points first.', cpfile);
end

regopts = load(fullfile(savepath, 'regopts.mat'));
regopts = regopts.opts;

cpdata = load(cpfile);
np1 = cellfun(@(x) size(x, 1), cpdata.atlas_control_points);
np2 = cellfun(@(x) size(x, 1), cpdata.histology_control_points);
ikeep = np1 == np2 & np1 > 0;
cptsatlas = cat(1, cpdata.atlas_control_points{ikeep});
cptshistology = cat(1, cpdata.histology_control_points{ikeep});
cptsatlas = cptsatlas(:, 1:3);
cptshistology = cptshistology(:, 1:3);

if isfield(cpdata, 'ori_trans')
    original_trans = affinetform3d(cpdata.ori_trans);
else
    original_trans = regopts.original_trans;
end

cptsatlas = cptsatlas(:, [2 1 3]);
cptshistology = cptshistology(:, [2 1 3]);
cptshistology = original_trans.transformPointsInverse(cptshistology);
cptsatlas = cptsatlas / regopts.downfac_reg;

fprintf('Python-exported control points in %s\n', cpfile);
fprintf('  matched slices: %d\n', nnz(ikeep));
fprintf('  manual pairs:     %d\n', size(cptshistology, 1));
fprintf('  atlas range:      [%g %g] x [%g %g] x [%g %g]\n', ...
    min(cptsatlas(:,1)), max(cptsatlas(:,1)), ...
    min(cptsatlas(:,2)), max(cptsatlas(:,2)), ...
    min(cptsatlas(:,3)), max(cptsatlas(:,3)));
fprintf('  sample range:     [%g %g] x [%g %g] x [%g %g]\n', ...
    min(cptshistology(:,1)), max(cptshistology(:,1)), ...
    min(cptshistology(:,2)), max(cptshistology(:,2)), ...
    min(cptshistology(:,3)), max(cptshistology(:,3)));
fprintf('Ready for multiobjRegistration — ensure regopts.mat and registration volumes are in %s\n', savepath);
end
