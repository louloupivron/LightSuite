function out = import_arivis_csv_cells(csvpath, lightsuite_savepath, varargin)
%IMPORT_ARIVIS_CSV_CELLS Import Arivis segmentation CSV as LightSuite cells and optionally register.
%
% Usage:
%   out = import_arivis_csv_cells(csvpath, lightsuite_savepath)
%   out = import_arivis_csv_cells(..., 'prefix', 'arivis_')
%   out = import_arivis_csv_cells(..., 'do_transform', true)
%
% Inputs
%   csvpath            - Arivis "features.csv" export containing X/Y/Z in pixels.
%   lightsuite_savepath- Folder containing regopts.mat and (optionally) transform_params.mat.
%
% Name-value
%   prefix       - filename prefix for saved outputs (default: 'arivis_')
%   do_transform - if true and transform_params.mat exists, run transformPointsToAtlas
%                  and save atlas-space points into volume_registered/ (default: true)
%
% Output (struct)
%   .cell_locations_sample_mat   - path written
%   .cell_locations_sample_csv   - path written (simple numeric)
%   .was_0_based                 - heuristic decision
%   .n_points                    - number of segments imported
%   .atlas_csv                   - written if do_transform=true
%   .atlas_mat                   - written if do_transform=true

ip = inputParser;
ip.addParameter('prefix', 'arivis_', @(x) ischar(x) || isstring(x));
ip.addParameter('do_transform', true, @(x) islogical(x) && isscalar(x));
ip.parse(varargin{:});
prefix = char(ip.Results.prefix);
do_transform = ip.Results.do_transform;

if ~exist(csvpath, 'file')
    error('CSV not found: %s', csvpath);
end
if ~exist(lightsuite_savepath, 'dir')
    error('LightSuite savepath not found: %s', lightsuite_savepath);
end

optfile = fullfile(lightsuite_savepath, 'regopts.mat');
if ~exist(optfile, 'file')
    error('regopts.mat not found in %s', lightsuite_savepath);
end
opts = load(optfile);
opts = opts.opts;

T = readtable(csvpath, 'TextType', 'string');

colX = find(contains(T.Properties.VariableNames, "X (px)"), 1, 'first');
colY = find(contains(T.Properties.VariableNames, "Y (px)"), 1, 'first');
colZ = find(contains(T.Properties.VariableNames, "Z (px)"), 1, 'first');
if isempty(colX) || isempty(colY) || isempty(colZ)
    error('Could not find X/Y/Z (px) center-of-mass columns in %s', csvpath);
end

xyz = [T{:, colX}, T{:, colY}, T{:, colZ}];
xyz = double(xyz);
xyz = xyz(all(isfinite(xyz), 2), :);

% LightSuite expects sample voxel coordinates in the original sample grid.
% Heuristic: detect whether Arivis exported 0-based pixel coordinates.
Nx = opts.Nx; Ny = opts.Ny; Nz = opts.Nz;
mins = min(xyz, [], 1);
maxs = max(xyz, [], 1);

fits_1_based = all(mins >= 1) && (maxs(1) <= Nx) && (maxs(2) <= Ny) && (maxs(3) <= Nz);
fits_0_based = all(mins >= 0) && (maxs(1) <= (Nx-1)) && (maxs(2) <= (Ny-1)) && (maxs(3) <= (Nz-1));

was_0_based = false;
if ~fits_1_based && fits_0_based
    xyz = xyz + 1;
    was_0_based = true;
elseif fits_1_based
    was_0_based = false;
elseif fits_0_based
    % Ambiguous (both could fit if coords are safely away from borders).
    % Prefer 1-based because LightSuite internal detections are 1-based.
    was_0_based = false;
else
    warning(['Arivis XYZ ranges do not fit regopts size. ' ...
        'X range [%g,%g] vs Nx=%d, Y range [%g,%g] vs Ny=%d, Z range [%g,%g] vs Nz=%d. ' ...
        'Proceeding without offset; registration may be incorrect.'], ...
        mins(1), maxs(1), Nx, mins(2), maxs(2), Ny, mins(3), maxs(3), Nz);
end

% Build LightSuite-style cell_locations: [x y z intensity diameter ellipticity].
% Arivis export may have volume/voxel count etc, but those are not used by mapping.
n = size(xyz, 1);
cell_locations = single([xyz, nan(n, 3)]);

% Save alongside other LightSuite artifacts so transformPointsToAtlas folder-mode can find it.
cell_locations_sample_mat = fullfile(lightsuite_savepath, sprintf('%scell_locations_sample.mat', prefix));
save(cell_locations_sample_mat, 'cell_locations');

cell_locations_sample_csv = fullfile(lightsuite_savepath, sprintf('%scell_locations_sample.csv', prefix));
writematrix(cell_locations, cell_locations_sample_csv);

out = struct();
out.cell_locations_sample_mat = cell_locations_sample_mat;
out.cell_locations_sample_csv = cell_locations_sample_csv;
out.was_0_based = was_0_based;
out.n_points = n;
out.atlas_csv = '';
out.atlas_mat = '';

trfile = fullfile(lightsuite_savepath, 'transform_params.mat');
if do_transform && exist(trfile, 'file')
    trstruct = load(trfile);
    atlaspts = transformPointsToAtlas(cell_locations, 'transform_params', trstruct);

    registerpath = fullfile(lightsuite_savepath, 'volume_registered');
    if ~exist(registerpath, 'dir'); mkdir(registerpath); end

    atlas_mat = fullfile(registerpath, sprintf('%scell_locations_atlas.mat', prefix));
    save(atlas_mat, 'atlaspts');

    atlas_csv = fullfile(registerpath, sprintf('%scell_locations_atlas.csv', prefix));
    writematrix(atlaspts, atlas_csv);

    out.atlas_csv = atlas_csv;
    out.atlas_mat = atlas_mat;
elseif do_transform
    warning('transform_params.mat not found in %s; skipping atlas transform.', lightsuite_savepath);
end

end

