"""Tests for division label map generation."""

from __future__ import annotations

from pathlib import Path

import nibabel as nib
import numpy as np
import pandas as pd
import pytest

from lightsuite.analysis.division_map import build_division_labels, ensure_division_map
from lightsuite.atlas.registry import AtlasPaths


def _write_allen_membership(path: Path) -> None:
    rows = []
    spec = {
        0: ("unassigned", "unassigned"),
        1: ("Isocortex", "Isocortex"),
        2: ("Thalamus", "TH"),
    }
    for pidx, (name, acr) in spec.items():
        rows.append(
            {
                "parcellation_index": pidx,
                "parcellation_term_set_name": "division",
                "parcellation_term_name": name,
                "parcellation_term_acronym": acr,
            }
        )
    pd.DataFrame(rows).to_csv(path, index=False)


def test_build_division_labels_allen_brainglobe_uses_ccf_ids(tmp_path: Path) -> None:
    """BrainGlobe Allen annotations store CCF ids, not ABC parcellation indices."""
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (2, 2, 2)
    ccf_isocortex = 320
    ccf_thalamus = 718
    annotation = np.array(
        [
            [[0, ccf_isocortex], [ccf_isocortex, 0]],
            [[ccf_thalamus, ccf_thalamus], [0, ccf_isocortex]],
        ],
        dtype=np.uint32,
    )
    nib.save(nib.Nifti1Image(annotation, np.eye(4)), str(atlas_dir / "annotation_10.nii.gz"))
    nib.save(nib.Nifti1Image(np.zeros(shape), np.eye(4)), str(atlas_dir / "average_template_10.nii.gz"))
    membership = atlas_dir / "parcellation_to_parcellation_term_membership.csv"
    pd.DataFrame(
        [
            {
                "parcellation_index": 1,
                "parcellation_term_set_name": "division",
                "parcellation_term_name": "Isocortex",
                "parcellation_term_acronym": "Isocortex",
                "parcellation_term_label": "AllenCCF-Ontology-2017-315",
            },
            {
                "parcellation_index": 1,
                "parcellation_term_set_name": "substructure",
                "parcellation_term_name": "Primary motor area, Layer 1",
                "parcellation_term_acronym": "MOp1",
                "parcellation_term_label": f"AllenCCF-Ontology-2017-{ccf_isocortex}",
            },
            {
                "parcellation_index": 2,
                "parcellation_term_set_name": "division",
                "parcellation_term_name": "Thalamus",
                "parcellation_term_acronym": "TH",
                "parcellation_term_label": "AllenCCF-Ontology-2017-549",
            },
            {
                "parcellation_index": 2,
                "parcellation_term_set_name": "substructure",
                "parcellation_term_name": "Ventral posterolateral nucleus",
                "parcellation_term_acronym": "VPL",
                "parcellation_term_label": f"AllenCCF-Ontology-2017-{ccf_thalamus}",
            },
        ]
    ).to_csv(membership, index=False)
    pd.DataFrame(
        {
            "id": [ccf_isocortex, ccf_thalamus],
            "parent_id": [985, 637],
            "acronym": ["MOp1", "VPL"],
            "name": ["Primary motor area, Layer 1", "Ventral posterolateral nucleus"],
        }
    ).to_csv(atlas_dir / "structures.csv", index=False)

    atlas = AtlasPaths(
        brain_atlas="allen",
        atlas_dir=atlas_dir,
        template_path=atlas_dir / "average_template_10.nii.gz",
        annotation_path=atlas_dir / "annotation_10.nii.gz",
        boundary_path=None,
        structures_csv_path=atlas_dir / "structures.csv",
        supports_parcellation=True,
        atlas_source="brainglobe",
        brainglobe_name="allen_mouse_10um",
    )
    labels, legend = build_division_labels(atlas)
    assert labels.shape == shape
    assert set(np.unique(labels)) <= {0, 1, 2}
    iso_id = int(legend.loc[legend["division_name"] == "Isocortex", "division_id"].iloc[0])
    th_id = int(legend.loc[legend["division_name"] == "Thalamus", "division_id"].iloc[0])
    assert np.any(labels == iso_id)
    assert np.any(labels == th_id)


def test_build_division_labels_allen(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (2, 2, 2)
    # voxels: 0=bg, 1=isocortex, 2=thalamus
    annotation = np.array(
        [
            [[0, 1], [1, 0]],
            [[2, 2], [0, 1]],
        ],
        dtype=np.uint16,
    )
    nib.save(nib.Nifti1Image(annotation, np.eye(4)), str(atlas_dir / "annotation_10.nii.gz"))
    nib.save(nib.Nifti1Image(np.zeros(shape), np.eye(4)), str(atlas_dir / "average_template_10.nii.gz"))
    _write_allen_membership(atlas_dir / "parcellation_to_parcellation_term_membership.csv")

    atlas = AtlasPaths(
        brain_atlas="allen",
        atlas_dir=atlas_dir,
        template_path=atlas_dir / "average_template_10.nii.gz",
        annotation_path=atlas_dir / "annotation_10.nii.gz",
        boundary_path=None,
        structures_csv_path=None,
        supports_parcellation=True,
    )
    labels, legend = build_division_labels(atlas)
    assert labels.shape == shape
    assert set(np.unique(labels)) <= {0, 1, 2}
    iso_id = int(legend.loc[legend["division_name"] == "Isocortex", "division_id"].iloc[0])
    assert np.any(labels == iso_id)


def test_ensure_division_map_caches(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (2, 2, 2)
    annotation = np.ones(shape, dtype=np.uint16)
    nib.save(nib.Nifti1Image(annotation, np.eye(4)), str(atlas_dir / "annotation_10.nii.gz"))
    nib.save(nib.Nifti1Image(np.zeros(shape), np.eye(4)), str(atlas_dir / "average_template_10.nii.gz"))
    _write_allen_membership(atlas_dir / "parcellation_to_parcellation_term_membership.csv")

    atlas = AtlasPaths(
        brain_atlas="allen",
        atlas_dir=atlas_dir,
        template_path=atlas_dir / "average_template_10.nii.gz",
        annotation_path=atlas_dir / "annotation_10.nii.gz",
        boundary_path=None,
        structures_csv_path=None,
        supports_parcellation=True,
    )
    first = ensure_division_map(atlas)
    assert first.paths.labels_tiff.is_file()
    assert first.paths.legend_csv.is_file()

    second = ensure_division_map(atlas)
    assert second.labels.shape == first.labels.shape


def test_build_division_labels_bad_membership_columns(tmp_path: Path) -> None:
    atlas_dir = tmp_path / "atlas"
    atlas_dir.mkdir()
    shape = (2, 2, 2)
    nib.save(nib.Nifti1Image(np.zeros(shape), np.eye(4)), str(atlas_dir / "annotation_10.nii.gz"))
    nib.save(nib.Nifti1Image(np.zeros(shape), np.eye(4)), str(atlas_dir / "average_template_10.nii.gz"))
    (atlas_dir / "parcellation_to_parcellation_term_membership.csv").write_text(
        "parcellation_index\n1\n", encoding="utf-8"
    )
    atlas = AtlasPaths(
        brain_atlas="allen",
        atlas_dir=atlas_dir,
        template_path=atlas_dir / "average_template_10.nii.gz",
        annotation_path=atlas_dir / "annotation_10.nii.gz",
        boundary_path=None,
        structures_csv_path=None,
        supports_parcellation=True,
    )
    with pytest.raises(ValueError, match="missing columns"):
        build_division_labels(atlas)
