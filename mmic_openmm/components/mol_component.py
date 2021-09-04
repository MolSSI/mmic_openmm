from mmelemental.models import Molecule
from typing import List, Tuple, Optional
from mmelemental.util.units import convert
from mmic_translator import (
    TransComponent,
    TransInput,
    TransOutput,
    __version__,
)
from simtk.openmm import app
import numpy

provenance_stamp = {
    "creator": "mmic_openmm",
    "version": __version__,
    "routine": __name__,
}

__all__ = ["MolToOpenMMComponent", "OpenMMToMolComponent"]


class MolToOpenMMComponent(TransComponent):
    """A component for converting MMSchema to OpenMM Molecule."""

    def execute(
        self,
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        mmol = inputs.schema_object
        ndim = mmol.ndim

        if ndim != 3:
            raise NotImplementedError("mmic_openmm supports only 3D molecules.")

        top = app.Topology()
        natoms = len(mmol.symbols)

        if mmol.masses is not None:
            masses = convert(
                mmol.masses, mmol.masses_units, "amu"  # check unit in OpenMM
            )
        else:
            raise TypeError("mmic_openmm required atomic masses to be defined.")

        if mmol.atom_labels is not None:
            atom_labels = mmol.atom_labels
        else:
            atom_labels = ["" for i in range(natoms)]

        if mmol.atomic_numbers is not None:
            atomic_numbers = mmol.atomic_numbers
        else:
            raise NotImplementedError(
                "mmic_openmm is supported only for atomic/molecular systems. Molecule.atomic_numbers must be defined."
            )

        chains = [top.addChain()]

        if isinstance(mmol.extras, dict):
            if mmol.extras.get("chains"):
                chains = [top.addChain() for i in range(len(mmol.chains))]

        if mmol.substructs is not None:
            residues = [
                top.addResidue(resname, chain=chains[0], id=resid)
                for resname, resid in mmol.get_substructs()
            ]
        atoms = [None] * natoms

        for index, symb in enumerate(mmol.symbols):

            label = atom_labels[index] if atom_labels is not None else symb

            atomic_number = (
                atomic_numbers[index] if atomic_numbers is not None else None
            )
            mass = masses[index] if masses is not None else None

            element = app.element.Element.getByAtomicNumber(atomic_numbers[index])

            atoms.append(
                top.addAtom(label, element, residues[mmol.substructs[index][1]])
            )

        if mmol.geometry is not None:
            coordinates = mmol.geometry.reshape(natoms, ndim)

        if mmol.connectivity is not None:
            for (
                i,
                j,
                order,
            ) in mmol.connectivity:
                top.addBond(atoms[i], atoms[j])

        omol = (
            {
                "data": top,
                "positions": coordinates,
                "positions_units": mmol.geometry_units,
            },
        )

        return True, TransOutput(
            proc_input=inputs,
            data_object=omol,
            success=True,
            provenance=provenance_stamp,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )


class OpenMMToMolComponent(TransComponent):
    """A component for converting OpenMM to MMSchema Molecule object."""

    def execute(
        self,
        inputs: TransInput,
        extra_outfiles: Optional[List[str]] = None,
        extra_commands: Optional[List[str]] = None,
        scratch_name: Optional[str] = None,
        timeout: Optional[int] = None,
    ) -> Tuple[bool, TransOutput]:

        if isinstance(inputs, dict):
            inputs = self.input()(**inputs)

        top = inputs.data_object
        geo = inputs.keywords.get("positions", None)
        if geo is not None:
            geo_units = inputs.keywords.get("positions_units", geo.unit.get_name())

            geo = numpy.array([(pos.x, pos.y, pos.z) for pos in geo]).T.flatten()

        atomic_data = [
            (atom.name, atom.element.atomic_number, atom.element.symbol)
            for atom in top.atoms()
        ]
        names, atomic_nums, element_names = zip(*atomic_data)

        try:
            masses = [atom.element.mass._value for atom in top.atoms()]
            masses_units = next(top.atoms()).element.mass.unit.get_name()
        except Exception:
            masses = None
            masses_units = "dalton"

        # If bond order is none, set it to 1.
        connectivity = [
            (bond.atom1.index, bond.atom2.index, bond.order or 1)
            for bond in top.bonds()
        ]

        residues = [(atom.residue.name, atom.residue.index) for atom in top.atoms()]

        input_dict = {
            "atomic_numbers": atomic_nums,
            "symbols": element_names,
            "atom_labels": names,
            "geometry": geo,
            "geometry_units": geo_units,
            "substructs": residues,
            "connectivity": connectivity,
            "masses": masses,
            "masses_units": masses_units,
        }

        return True, TransOutput(
            proc_input=inputs,
            schema_object=Molecule(**input_dict),
            success=True,
            provenance=provenance_stamp,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )
