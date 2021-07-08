"""
mmic_openmm
Tactic MMIC for OpenMM/MMSchema translation.
"""

# Add imports here
from . import components, models

# Handle versioneer
from ._version import get_versions

versions = get_versions()
__version__ = versions["version"]
__git_revision__ = versions["full-revisionid"]
del get_versions, versions

# Need to update these lists
molread_ext_maps = {
    ".gro": "gro",
    ".pdb": "pdb",
    ".inpcrd": "inpcrd",
}

molwrite_ext_maps = {".pdb": "pdb"}

ffread_ext_maps = {".psf": "psf", ".top": "top", ".prmtop": "prmtop"}
ffwrite_ext_maps = {".psf": "psf", ".top": "top", ".prmtop": "prmtop"}

_classes_map = {
    "Molecule": models.OpenMMMol,
    "ForceField": models.OpenMMFF,
}
