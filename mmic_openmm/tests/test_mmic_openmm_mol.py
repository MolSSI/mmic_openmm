"""
Unit and regression test for the mmic_openmm package.
"""

# Import package, test suite, and other packages as needed
import mmic_openmm
import pytest
import sys
import os
from simtk.openmm.app import (
    PDBFile,
    GromacsGroFile,
    GromacsTopFile,
    AmberPrmtopFile,
    AmberInpcrdFile,
)
import mmelemental as mm
import mm_data

pdbfile = mm_data.mols["dialanine.pdb"]
grofile = mm_data.mols["dialanine.gro"]
inpcrdfile = mm_data.mols["dialanine.inpcrd"]
topfile = mm_data.ffs["dialanine.top"]


def test_mmic_openmm_imported():
    """Sample test, will always pass so long as import statement worked"""
    assert "mmic_openmm" in sys.modules


def test_mmic_to_mol_from_pdb(**kwargs):
    struct = PDBFile(pdbfile)
    kwargs.setdefault("positions", struct.positions)
    kwargs.setdefault("positions_units", struct.positions.unit.get_name())
    inputs = {"data_object": struct.topology, "keywords": kwargs}

    return mmic_openmm.components.OpenMMToMolComponent.compute(inputs)


def test_mmic_to_mol_from_gro(**kwargs):
    struct = GromacsGroFile(grofile)
    gro = GromacsTopFile(topfile)
    kwargs.setdefault("positions", struct.positions)
    kwargs.setdefault("positions_units", struct.positions.unit.get_name())
    inputs = {"data_object": gro.topology, "keywords": kwargs}
    return mmic_openmm.components.OpenMMToMolComponent.compute(inputs)


def test_mol_to_openmm(**kwargs):
    mmol = mm.models.molecule.mm_mol.Molecule.from_file(grofile, topfile)
    inputs = {"schema_object": mmol, "keywords": kwargs}
    return mmic_openmm.components.MolToOpenMMComponent.compute(inputs)


def test_io_methods(**kwargs):
    omol = mmic_openmm.models.OpenMMMol.from_file(pdbfile, topfile)
    assert isinstance(omol.data, omol.dtype)

    omol.to_file("tmp.pdb")
    os.remove("tmp.pdb")
    # mmol = omol.to_schema()
    # assert isinstance(mmol, mm.models.molecule.Molecule)
