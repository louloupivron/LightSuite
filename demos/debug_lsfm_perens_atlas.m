%DEBUG_LSFM_PERENS_ATLAS Quick visualization of Gubra LSFM atlas files (MATLAB).
%
%   Set atlas_root to your LSFM-mouse-brain-atlas-master folder, run this demo,
%   then inspect the figure handles. Uses FINDLSFMATLASNIFTIS to resolve the split
%   layout where the template may live under LSFM_atlas_files/perens/.
%
%   Requires niftiread (common MATLAB releases). Prefer Python inspect_lsfm_atlas.py
%   for CSV-vs-volume statistics and PNG exports.

atlas_root = ''; % set to full path of LSFM-mouse-brain-atlas-master, or leave empty to pick folder

if isempty(atlas_root)
    atlas_root = uigetdir('', 'Select LSFM-mouse-brain-atlas-master folder');
    if atlas_root == 0
        error('LightSuite:debug_lsfm', 'Cancelled.');
    end
end

paths = findLsfmAtlasNiftis(atlas_root);
if isfield(paths, 'note')
    fprintf('%s\n', paths.note);
end

tv = niftiread(paths.template_path);
av = niftiread(paths.annotation_path);

assert(isequal(size(tv), size(av)), 'LightSuite:debug_lsfm', ...
    'Template and annotation must share the same grid (%s vs %s).', ...
    mat2str(size(tv)), mat2str(size(av)));

fprintf('Shape %s — template %s, annotation %s\n', mat2str(size(tv)), ...
    class(tv), class(av));

% Mid-slice montage: rows = template / annotation; cols = cuts along dim 1..3
figure('Color', 'w', 'Name', 'LSFM atlas mid-slices');
tiledlayout(2, 3, 'Padding', 'compact', 'TileSpacing', 'compact');

sz = size(tv);
mid = ceil(sz ./ 2);

for col = 1:3
    nexttile(col);
    if col == 1
        sl = squeeze(tv(mid(1), :, :));
        ttl = 'axis 1 (dim 1 fixed)';
    elseif col == 2
        sl = squeeze(tv(:, mid(2), :));
        ttl = 'axis 2 (dim 2 fixed)';
    else
        sl = squeeze(tv(:, :, mid(3)));
        ttl = 'axis 3 (dim 3 fixed)';
    end
    imagesc(sl); axis image off; colormap(gca, gray);
    title(ttl);
end

for col = 1:3
    nexttile(col + 3);
    if col == 1
        sl = squeeze(av(mid(1), :, :));
    elseif col == 2
        sl = squeeze(av(:, mid(2), :));
    else
        sl = squeeze(av(:, :, mid(3)));
    end
    imagesc(sl); axis image off; colormap(gca, parula);
    title('annotation');
end

sgtitle(sprintf('LSFM Perens atlas (%s)', paths.template_path));
