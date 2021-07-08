from mmelemental.models.molecule import Molecule
from typing import List, Tuple, Optional
from mmelemental.util.units import convert
from mmic_translator import (
    TransComponent,
    TransInput,
    TransOutput,
)
from simtk.openmm import app

__all__ = ["MolToOpenMMComponent", "OpenMMToMolComponent"]


class MolToOpenMMComponent(TransComponent):
    """A component for converting Molecule to ParmEd molecule object."""

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

        if mmol.extras.get("chains"):
            chains = [top.addChain() for i in range(len(mmol.chains))]
        else:
            chains = [top.addChain()]

        if mmol.substructs:
            residues = list(_fast_set(mmol.substructs))
            nres = len(residues)
            resnames, _ = zip(*residues)
            _, resids = zip(*mmol.substructs)
            resids = [i for i in resids]
        else:
            nres = 1
            resnames, resids = "UNK", [0]

        residues = [
            top.addResidue(resname, chain=chains[0], id=resid)
            for resname, resid in residues
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

        if mmol.connectivity:
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
        )


class OpenMMToMolComponent(TransComponent):
    """A component for converting ParmEd molecule to Molecule object."""

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

        # I think parmed.Structure does not store forces
        pmol = inputs.data_object
        geo_units, vel_units = None, None

        geo = TransComponent.get(pmol, "coordinates")
        if geo is not None:  # General enough? hackish?
            geo = geo.flatten()
            geo_units = pmol.positions.unit.get_name()

        vel = TransComponent.get(pmol, "velocities")
        if vel is not None:
            vel = vel.flatten()
            vel_units = "angstrom/picosecond"  # hard-coded in ParmEd

        atomic_nums = [atom.atomic_number for atom in pmol.atoms]
        names = [atom.name for atom in pmol.atoms]
        element_names = [atom.element_name for atom in pmol.atoms]

        masses = [atom.mass for atom in pmol.atoms]
        masses_units = pmol.atoms[0].umass.unit.get_name()

        # If bond order is none, set it to 1.
        if hasattr(pmol, "bonds"):
            connectivity = [
                (bond.atom1.idx, bond.atom2.idx, bond.order or 1) for bond in pmol.bonds
            ]
        else:
            connectivity = None

        if hasattr(pmol, "residues"):
            residues = [(atom.residue.name, atom.residue.idx) for atom in pmol.atoms]

        input_dict = {
            "atomic_numbers": atomic_nums,
            "symbols": element_names,
            "atom_labels": names,
            "geometry": geo,
            "geometry_units": geo_units,
            "velocities": vel,
            "velocities_units": vel_units,
            "substructs": residues,
            "connectivity": connectivity,
            "masses": masses,
            "masses_units": masses_units,
        }

        return True, TransOutput(
            proc_input=inputs, schema_object=Molecule(**input_dict)
        )


def _fast_set(seq: List) -> List:
    """Removes duplicate entries in a list while preserving the order."""
    seen = set()
    seen_add = seen.add
    return [x for x in seq if not (x in seen or seen_add(x))]
