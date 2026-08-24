"""Manual rostrocaudal orientation for spinal cord (replaces auto tofliprc detection).

The registration longitudinal axis is volume dimension 3 (Z, put last during
preprocess). ``cord_orientation.txt`` records the anatomical direction of the sample
along +Z so the atlas can be flipped to match the sample's storage order:

- ``rostrocaudal`` — +Z runs rostral → caudal (first slice is rostral/brain side).
  Atlas is left in native order (``tofliprc = False``).
- ``caudorostral`` — +Z runs caudal → rostral (last slice is rostral/brain side).
  Atlas rostrocaudal axis is flipped (``tofliprc = True``).

The orientation is set once with ``lightsuite spinal check-orientation`` (a small
Napari window) or by hand in ``cord_orientation.txt`` / the YAML config, and consumed
by ``preprocess``.
"""

from __future__ import annotations

from pathlib import Path

ROSTROCAUDAL = "rostrocaudal"
CAUDOROSTRAL = "caudorostral"
CORD_DIRECTIONS = (ROSTROCAUDAL, CAUDOROSTRAL)
CORD_ORIENTATION_FILENAME = "cord_orientation.txt"


class CordOrientationRequiredError(FileNotFoundError):
    """Raised when cord orientation must be confirmed before preprocess can finish."""


def cord_orientation_path(save_path: Path) -> Path:
    return save_path.expanduser() / CORD_ORIENTATION_FILENAME


def normalize_direction(direction: str) -> str:
    """Validate and canonicalize a longitudinal direction string."""
    value = str(direction).strip().lower()
    if value not in CORD_DIRECTIONS:
        msg = (
            f"Longitudinal direction must be one of {CORD_DIRECTIONS}, got {direction!r}."
        )
        raise ValueError(msg)
    return value


def tofliprc_from_direction(direction: str) -> bool:
    """Caudorostral samples need the atlas flipped along its rostrocaudal axis."""
    return normalize_direction(direction) == CAUDOROSTRAL


def direction_from_tofliprc(tofliprc: bool) -> str:
    return CAUDOROSTRAL if tofliprc else ROSTROCAUDAL


def save_cord_orientation(
    save_path: Path,
    direction: str,
    *,
    source: str = "manual",
) -> Path:
    """Write ``cord_orientation.txt`` with a human-readable, editable layout."""
    direction = normalize_direction(direction)
    tofliprc = tofliprc_from_direction(direction)
    path = cord_orientation_path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    text = (
        "# LightSuite spinal cord longitudinal orientation\n"
        "# direction: anatomical direction of the sample along the registration +Z axis\n"
        "#   rostrocaudal = first slice is rostral (brain side) -> tofliprc = false\n"
        "#   caudorostral = last slice is rostral  (brain side) -> tofliprc = true\n"
        "# Edit 'direction' and re-run 'lightsuite spinal preprocess' to apply.\n"
        f"direction: {direction}\n"
        f"tofliprc: {str(tofliprc).lower()}\n"
        f"source: {source}\n"
    )
    path.write_text(text, encoding="utf-8")
    return path


def load_cord_orientation(save_path: Path) -> str | None:
    """Return the stored direction from ``cord_orientation.txt`` (None when absent)."""
    path = cord_orientation_path(save_path)
    if not path.is_file():
        return None
    fields: dict[str, str] = {}
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        key, _, value = line.partition(":")
        fields[key.strip().lower()] = value.strip()

    direction = fields.get("direction")
    if direction is None and "tofliprc" in fields:
        direction = direction_from_tofliprc(
            fields["tofliprc"].lower() in ("true", "1", "yes")
        )
    if direction is None:
        return None
    try:
        return normalize_direction(direction)
    except ValueError:
        return None


def resolve_cord_orientation(
    save_path: Path,
    *,
    config_direction: str | None = None,
    require: bool = False,
) -> str:
    """Resolve the longitudinal direction from YAML, ``cord_orientation.txt``, or fail.

    Precedence: explicit ``config_direction`` (YAML) > ``cord_orientation.txt``.
    When neither is available and ``require`` is True, raise ``FileNotFoundError``.
    """
    if config_direction is not None:
        return normalize_direction(config_direction)

    direction = load_cord_orientation(save_path)
    if direction is not None:
        return direction

    if require:
        path = cord_orientation_path(save_path)
        msg = (
            f"Missing {path}. Run 'lightsuite spinal preprocess' (opens the orientation "
            "GUI when needed), 'lightsuite spinal check-orientation', or set "
            "registration.longitudinal_direction in the YAML config."
        )
        raise CordOrientationRequiredError(msg)
    return ROSTROCAUDAL
