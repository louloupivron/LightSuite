"""Atlas content trim with scratch caching (Option A — working grid only)."""

from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import tifffile

from lightsuite.atlas.io import load_atlas_volume
from lightsuite.atlas.registry import AtlasPaths
from lightsuite.config.models import AtlasConfig, ContentTrimMode
from lightsuite.registration.content_bbox import (
    ContentBox,
    bbox_from_annotation,
    bbox_from_template_threshold,
    crop_volume_yxz,
)


ATLAS_TRIM_MANIFEST = "atlas_content_manifest.json"


@dataclass(frozen=True)
class AtlasContentManifest:
    """Trim metadata for a cached atlas working copy."""

    brain_atlas: str
    resolution_um: float
    native_shape: list[int]
    crop_start: list[int]
    crop_size: list[int]
    source_paths: dict[str, str]
    source_fingerprint: str

    @property
    def crop_start_yxz(self) -> tuple[int, int, int]:
        return tuple(int(v) for v in self.crop_start)

    @property
    def is_trimmed(self) -> bool:
        return any(int(v) > 0 for v in self.crop_start) or list(self.crop_size) != list(
            self.native_shape
        )

    def to_dict(self) -> dict:
        return {
            "brain_atlas": self.brain_atlas,
            "resolution_um": self.resolution_um,
            "native_shape": self.native_shape,
            "crop_start": self.crop_start,
            "crop_size": self.crop_size,
            "source_paths": self.source_paths,
            "source_fingerprint": self.source_fingerprint,
        }

    @classmethod
    def from_dict(cls, raw: dict) -> AtlasContentManifest:
        return cls(
            brain_atlas=str(raw["brain_atlas"]),
            resolution_um=float(raw["resolution_um"]),
            native_shape=[int(v) for v in raw["native_shape"]],
            crop_start=[int(v) for v in raw["crop_start"]],
            crop_size=[int(v) for v in raw["crop_size"]],
            source_paths={str(k): str(v) for k, v in raw["source_paths"].items()},
            source_fingerprint=str(raw["source_fingerprint"]),
        )

    def save(self, path: Path) -> None:
        path = path.expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(json.dumps(self.to_dict(), indent=2), encoding="utf-8")

    @classmethod
    def load(cls, path: Path) -> AtlasContentManifest:
        return cls.from_dict(json.loads(path.expanduser().read_text(encoding="utf-8")))


@dataclass(frozen=True)
class ResolvedAtlasContent:
    """Atlas paths for registration/GUI plus native-space trim offsets (Option A)."""

    paths: AtlasPaths
    manifest: AtlasContentManifest | None
    native_shape: tuple[int, int, int]
    crop_start_yxz: tuple[int, int, int]

    @property
    def is_trimmed(self) -> bool:
        return self.manifest is not None and self.manifest.is_trimmed


def _source_fingerprint(atlas: AtlasPaths) -> str:
    parts: list[str] = []
    for key in ("template", "annotation", "boundary"):
        path = {
            "template": atlas.template_path,
            "annotation": atlas.annotation_path,
            "boundary": atlas.boundary_path,
        }[key]
        if path is None or not path.is_file():
            parts.append(f"{key}:missing")
            continue
        stat = path.stat()
        parts.append(f"{key}:{path}:{stat.st_mtime_ns}:{stat.st_size}")
    digest = hashlib.sha256("|".join(parts).encode()).hexdigest()
    return digest[:16]


def _trim_cache_root(scratch: Path, atlas: AtlasPaths, resolution_um: float) -> Path:
    tag = atlas.brainglobe_name or atlas.brain_atlas
    safe = tag.replace("/", "_").replace(" ", "_")
    return scratch.expanduser() / "atlas_trim" / f"{safe}_{int(resolution_um)}um"


def resolve_atlas_content_box(
    annotation: np.ndarray,
    template: np.ndarray,
    cfg: AtlasConfig,
) -> ContentBox:
    mode = cfg.content_trim
    if mode == ContentTrimMode.MANUAL:
        if cfg.content_box is None:
            msg = "atlas.content_box is required when content_trim is manual"
            raise ValueError(msg)
        return ContentBox.from_manual_box(cfg.content_box)

    box = bbox_from_annotation(annotation, margin_vox=0)
    if box is None:
        box = bbox_from_template_threshold(template, margin_vox=0)
    if box is None:
        msg = "Could not detect atlas foreground for content_trim=auto"
        raise RuntimeError(msg)
    return ContentBox(
        y0=max(0, box.y0 - cfg.content_margin_vox),
        y1=min(annotation.shape[0] - 1, box.y1 + cfg.content_margin_vox),
        x0=max(0, box.x0 - cfg.content_margin_vox),
        x1=min(annotation.shape[1] - 1, box.x1 + cfg.content_margin_vox),
        z0=max(0, box.z0 - cfg.content_margin_vox),
        z1=min(annotation.shape[2] - 1, box.z1 + cfg.content_margin_vox),
    )


def _write_trimmed_volume(path: Path, volume: np.ndarray) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        path.unlink()
    name = path.name.lower()
    if name.endswith((".tif", ".tiff")):
        tifffile.imwrite(path, np.ascontiguousarray(volume))
        return
    import nibabel as nib

    img = nib.Nifti1Image(np.ascontiguousarray(volume), np.eye(4))
    nib.save(img, str(path))


def trim_atlas_to_cache(
    atlas: AtlasPaths,
    cfg: AtlasConfig,
    *,
    scratch: Path,
) -> ResolvedAtlasContent:
    """Load, optionally trim, and cache atlas volumes for registration."""
    template = load_atlas_volume(atlas.template_path)
    annotation = load_atlas_volume(atlas.annotation_path)
    native_shape = tuple(int(v) for v in template.shape)
    fingerprint = _source_fingerprint(atlas)

    if cfg.content_trim == ContentTrimMode.OFF:
        return ResolvedAtlasContent(
            paths=atlas,
            manifest=None,
            native_shape=native_shape,
            crop_start_yxz=(0, 0, 0),
        )

    box = resolve_atlas_content_box(annotation, template, cfg)
    if box.size_yxz == native_shape and box.start_yxz == (0, 0, 0):
        return ResolvedAtlasContent(
            paths=atlas,
            manifest=None,
            native_shape=native_shape,
            crop_start_yxz=(0, 0, 0),
        )

    cache_dir = _trim_cache_root(scratch, atlas, cfg.resolution_um)
    manifest_path = cache_dir / ATLAS_TRIM_MANIFEST
    if manifest_path.is_file():
        manifest = AtlasContentManifest.load(manifest_path)
        if manifest.source_fingerprint == fingerprint:
            trimmed_paths = AtlasPaths(
                brain_atlas=atlas.brain_atlas,
                atlas_dir=cache_dir,
                template_path=cache_dir / "template.tif",
                annotation_path=cache_dir / "annotation.tif",
                boundary_path=cache_dir / "boundary.tif"
                if (cache_dir / "boundary.tif").is_file()
                else None,
                structures_csv_path=atlas.structures_csv_path,
                supports_parcellation=atlas.supports_parcellation,
                atlas_source=atlas.atlas_source,
                brainglobe_name=atlas.brainglobe_name,
            )
            return ResolvedAtlasContent(
                paths=trimmed_paths,
                manifest=manifest,
                native_shape=native_shape,
                crop_start_yxz=manifest.crop_start_yxz,
            )

    cache_dir.mkdir(parents=True, exist_ok=True)
    tpl_crop = crop_volume_yxz(template, box)
    ann_crop = crop_volume_yxz(annotation, box)
    tpl_out = cache_dir / "template.tif"
    ann_out = cache_dir / "annotation.tif"
    _write_trimmed_volume(tpl_out, tpl_crop)
    _write_trimmed_volume(ann_out, ann_crop)

    source_paths = {
        "template": str(atlas.template_path),
        "annotation": str(atlas.annotation_path),
    }
    boundary_out: Path | None = None
    if atlas.boundary_path is not None and atlas.boundary_path.is_file():
        boundary = load_atlas_volume(atlas.boundary_path)
        boundary_out = cache_dir / "boundary.tif"
        _write_trimmed_volume(boundary_out, crop_volume_yxz(boundary, box))
        source_paths["boundary"] = str(atlas.boundary_path)

    manifest = AtlasContentManifest(
        brain_atlas=atlas.brain_atlas,
        resolution_um=float(cfg.resolution_um),
        native_shape=list(native_shape),
        crop_start=list(box.start_yxz),
        crop_size=list(box.size_yxz),
        source_paths=source_paths,
        source_fingerprint=fingerprint,
    )
    manifest.save(manifest_path)

    trimmed_paths = AtlasPaths(
        brain_atlas=atlas.brain_atlas,
        atlas_dir=cache_dir,
        template_path=tpl_out,
        annotation_path=ann_out,
        boundary_path=boundary_out,
        structures_csv_path=atlas.structures_csv_path,
        supports_parcellation=atlas.supports_parcellation,
        atlas_source=atlas.atlas_source,
        brainglobe_name=atlas.brainglobe_name,
    )
    return ResolvedAtlasContent(
        paths=trimmed_paths,
        manifest=manifest,
        native_shape=native_shape,
        crop_start_yxz=box.start_yxz,
    )


def save_atlas_manifest_copy(manifest: AtlasContentManifest, save_path: Path) -> Path:
    """Write ``atlas_content_manifest.json`` beside sample checkpoints."""
    out = save_path.expanduser() / ATLAS_TRIM_MANIFEST
    manifest.save(out)
    return out
