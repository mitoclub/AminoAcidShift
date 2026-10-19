"""Repository paths.

Analysis scripts and notebooks should not depend on the process working
directory. Shared inputs live under ``data/``. Sibling projects stay outside
this repository and can be overridden with environment variables.
"""

import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
DATA = REPO_ROOT / "data"
EXTERNAL = DATA / "external"
REFERENCE = DATA / "reference"

# Parent of this repository: ../192, ../human_mtDNA_megatree,
# ../mtdnaMutSpecOfCancers, ../nemu-pipeline
_PARENT = REPO_ROOT.parent


def _env_path(name: str, default: Path) -> Path:
    override = os.environ.get(name)
    return Path(override) if override else default


NEMU_192 = _env_path("AAS_NEMU_192", _PARENT / "192")
MEGATREE = _env_path("AAS_MEGATREE", _PARENT / "human_mtDNA_megatree")
CANCER_MUTSPEC = _env_path("AAS_CANCER_MUTSPEC", _PARENT / "mtdnaMutSpecOfCancers")
NEMU_PIPELINE = _env_path("AAS_NEMU_PIPELINE", _PARENT / "nemu-pipeline")

VERTEBRATES_OBS = REPO_ROOT / "analyses" / "02_vertebrates" / "dataset" / "obs_muts.csv"
