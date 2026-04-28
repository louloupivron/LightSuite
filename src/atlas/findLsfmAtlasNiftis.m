function paths = findLsfmAtlasNiftis(atlas_root)
%FINDLSFMATLASNIFTIS Locate Perens LSFM template + annotation under Gubra repo layout.
%
%   Some downloads place gubra_template_olf.nii.gz only under LSFM_atlas_files/perens/
%   while gubra_ano_olf.nii.gz sits in LSFM_atlas_files/. LightSuite expects both files
%   in the same opts.atlas_dir — copy/symlink one next to the other, or pass atlas_root
%   here for diagnostics that resolve both paths.
%
%   PATHS = FINDLSFMATLASNIFTIS(ATLAS_ROOT)
%   ATLAS_ROOT is the folder containing LSFM_atlas_files (e.g. LSFM-mouse-brain-atlas-master).
%
%   Output struct fields:
%     template_path, annotation_path, structures_csv (may be empty),
%     suggested_atlas_dir — folder to use as opts.atlas_dir after copying both NIfTIs there.

atlas_root = char(atlas_root);
lf = fullfile(atlas_root, 'LSFM_atlas_files');
assert(isfolder(lf), 'LightSuite:FindLsfmAtlas', ...
    'Expected LSFM_atlas_files under %s', atlas_root);

ann = fullfile(lf, 'gubra_ano_olf.nii.gz');
tpl = fullfile(lf, 'gubra_template_olf.nii.gz');
if ~isfile(tpl)
    tpl_alt = fullfile(lf, 'perens', 'gubra_template_olf.nii.gz');
    if isfile(tpl_alt)
        tpl = tpl_alt;
    end
end
assert(isfile(ann), 'LightSuite:FindLsfmAtlas', 'Missing %s', ann);
assert(isfile(tpl), 'LightSuite:FindLsfmAtlas', ...
    'Missing gubra_template_olf.nii.gz in LSFM_atlas_files or LSFM_atlas_files/perens');

csv = fullfile(lf, 'ARA2_annotation_info_avail_regions.csv');
if ~isfile(csv)
    csv = '';
end

paths = struct();
paths.template_path = tpl;
paths.annotation_path = ann;
paths.structures_csv = csv;

ann_dir = fileparts(ann);
tpl_dir = fileparts(tpl);
if strcmp(ann_dir, tpl_dir)
    paths.suggested_atlas_dir = ann_dir;
else
    paths.suggested_atlas_dir = ann_dir;
    paths.note = sprintf([ ...
        'Template and annotation are in different folders. Copy or symlink both into one folder ' ...
        'and set opts.atlas_dir to that folder:\n  annotation: %s\n  template: %s\n'], ...
        ann, tpl);
end
end
