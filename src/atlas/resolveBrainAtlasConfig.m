function cfg = resolveBrainAtlasConfig(opts)
%RESOLVEBRAINATLASCONFIG Resolve brain atlas folder and NIfTI paths from opts.
%
%   opts.brain_atlas - Atlas id (default 'allen'):
%       'allen'        - Allen CCF v3 (10 um): average_template_10.nii.gz,
%                        annotation_10.nii.gz
%       'perens'       - Perens et al. ~2020 / Gubra LSFM mouse brain atlas
%                        (20 um, orientation ial in BrainGlobe packaging):
%                        gubra_template_olf.nii.gz, gubra_ano_olf.nii.gz
%
%   opts.atlas_dir   - Optional char/string. Folder that contains the template
%                      and annotation files for the selected atlas. If empty,
%                      the folder is found via which() on the template filename
%                      (add that atlas directory to the MATLAB path).
%                      For split Gubra layouts (template under LSFM_atlas_files/perens/),
%                      copy both NIfTIs into one folder or use FINDLSFMATLASNIFTIS.
%
%   Output cfg fields: brain_atlas, atlas_dir, template_path, annotation_path,
%   template_file, annotation_file, structures_csv_path (Perens: ARA2 CSV if present),
%   supports_parcellation (Allen always; Perens when ARA2 CSV is next to the NIfTIs).

if nargin < 1 || isempty(opts)
    opts = struct();
end

brain_atlas = lower(strtrim(char(getOr(opts, 'brain_atlas', 'allen'))));
valid = {'allen', 'perens'};
assert(any(strcmp(brain_atlas, valid)), ...
    'LightSuite:BrainAtlas', ...
    'opts.brain_atlas must be one of: %s (got "%s").', ...
    strjoin(valid, ', '), brain_atlas);

switch brain_atlas
    case 'allen'
        tpl = 'average_template_10.nii.gz';
        ann = 'annotation_10.nii.gz';
    case 'perens'
        tpl = 'gubra_template_olf.nii.gz';
        ann = 'gubra_ano_olf.nii.gz';
end

if isfield(opts, 'atlas_dir') && ~isempty(opts.atlas_dir)
    atlas_dir = char(opts.atlas_dir);
    assert(isfolder(atlas_dir), ...
        'LightSuite:BrainAtlas', ...
        'opts.atlas_dir is not a folder: %s', atlas_dir);
else
    atlas_dir = fileparts(which(tpl));
end

if isempty(atlas_dir)
    error('LightSuite:AtlasNotFound', ...
        ['Brain atlas "%s" not found. Add the folder containing %s to the ' ...
        'MATLAB path, or set opts.atlas_dir to that folder.'], ...
        brain_atlas, tpl);
end

cfg.brain_atlas = brain_atlas;
cfg.atlas_dir = atlas_dir;
cfg.template_file = tpl;
cfg.annotation_file = ann;
cfg.template_path = fullfile(atlas_dir, tpl);
cfg.annotation_path = fullfile(atlas_dir, ann);

cfg.structures_csv_path = '';
if strcmp(brain_atlas, 'perens')
    cf = fullfile(atlas_dir, 'ARA2_annotation_info_avail_regions.csv');
    if isfile(cf)
        cfg.structures_csv_path = cf;
    end
end

cfg.supports_parcellation = strcmp(brain_atlas, 'allen') ...
    || (strcmp(brain_atlas, 'perens') && ~isempty(cfg.structures_csv_path));
end
