from mmelemental.models import forcefield
from mmic_translator.models.io import (
    TransInput,
    TransOutput,
)
from mmic_translator.components import TransComponent
from typing import List, Tuple, Optional
from collections.abc import Iterable
from mmelemental.util.units import convert
from simtk.openmm.app import forcefield as openmm_ff
import simtk

__all__ = ["FFToOpenMMComponent", "OpenMMToFFComponent"]


class FFToOpenMMComponent(TransComponent):
    """A component for converting MMSchema to OpenMM ForceField object."""

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

        empty_atom = parmed.topologyobjects.Atom()
        mmff = inputs.schema_object
        pff = parmed.structure.Structure()

        masses = convert(
            mmff.masses, mmff.masses_units, empty_atom.umass.unit.get_symbol()
        )

        charges = TransComponent.get(mmff, "charges")
        charges = convert(
            charges, mmff.charges_units, empty_atom.ucharge.unit.get_symbol()
        )

        atomic_numbers = TransComponent.get(mmff, "atomic_numbers")
        atom_types = TransComponent.get(mmff, "defs")

        rmin, epsilon = self._get_nonbonded(mmff, empty_atom)

        for index, symb in enumerate(mmff.symbols):

            # Will likely lose FF-related info ... but then Molecule is not supposed to store any params specific to FFs
            if atomic_numbers is not None:
                atomic_number = atomic_numbers[index]
            else:
                atomic_number = None

            if atom_types is not None:
                atom_type = atom_types[index]
            else:
                atom_type = None

            if masses is not None:
                mass = masses[index]
            else:
                mass = None

            if charges is not None:
                charge = charges[index]
            else:
                charge = None

            atom = parmed.topologyobjects.Atom(
                list=None,
                atomic_number=atomic_number,
                name=symb,
                type=atom_type,
                mass=mass,
                charge=charge,
                nb_idx=0,
                solvent_radius=0.0,
                screen=0.0,
                tree="BLA",
                join=0.0,
                irotat=0.0,
                occupancy=1.0,
                bfactor=0.0,
                altloc="",
                number=-1,
                rmin=rmin[index],
                epsilon=epsilon[index],
                rmin14=None,
                epsilon14=None,
                # bonds=..., faster than connecting atoms one by one as done below?
                # angles=...,
                # dihedrals=...,
                # impropers=...,
                # polarizable=...,
            )

            residues = TransComponent.get(mmff, "substructs")

            if residues:
                resname, resnum = residues[index]
            else:
                raise NotImplementedError(
                    "Residues must be supplied for forcefields based on atom typing."
                )

            pff.add_atom(atom, resname, resnum, chain="", inscode="", segid="")

        # Bonds
        bonds = TransComponent.get(mmff, "bonds")
        if bonds is not None:
            assert (
                mmff.bonds.form == "Harmonic"
            ), "Only Harmonic potential supported for now"

            spring = convert(
                bonds.params.spring, bonds.params.spring_units, "kcal/mol/angstroms**2"
            )
            req = convert(bonds.lengths, bonds.lengths_units, "angstroms")

            for (
                bi,
                (
                    i,
                    j,
                    order,
                ),
            ) in enumerate(mmff.bonds.indices):
                btype = parmed.topologyobjects.BondType(
                    k=spring[bi], req=req[bi], list=pff.bond_types
                )
                pff.bonds.append(
                    parmed.topologyobjects.Bond(
                        pff.atoms[i], pff.atoms[j], order=order, type=btype
                    )
                )
                pff.bond_types.append(btype)
                # both implementations seem to perform almost the same:
                # pff.atoms[i].bond_to(pff.atoms[j])

        # Angles
        angles = TransComponent.get(mmff, "angles")
        if angles is not None:
            assert (
                mmff.angles.form == "Harmonic"
            ), "Only Harmonic potential supported for now"

            spring = convert(
                angles.params.spring, angles.params.spring_units, "kcal/mol/radians^2"
            )
            angles_eq = convert(angles.angles, angles.angles_units, "degrees")

            for ai, (i, j, k) in enumerate(mmff.angles.indices):
                atype = parmed.topologyobjects.AngleType(
                    k=spring[ai], theteq=angles_eq[ai], list=pff.angle_types
                )
                pff.angles.append(
                    parmed.topologyobjects.Angle(
                        pff.atoms[i], pff.atoms[j], pff.atoms[k], type=atype
                    )
                )
                pff.angle_types.append(atype)

        # Dihedrals
        dihedrals = TransComponent.get(mmff, "dihedrals")
        if dihedrals is not None:
            dihedrals = (
                dihedrals.pop() if isinstance(dihedrals, list) else dihedrals
            )  # For now, keep track of only a single type
            # Need to change this ASAP! Must take multiple types into account!
            assert (
                dihedrals.form == "Charmm" or dihedrals.form == "CharmmMulti"
            ), "Only Charmm-style potentials supported for now"

            energy = convert(
                dihedrals.params.energy, dihedrals.params.energy_units, "kcal/mol"
            )
            phase = convert(
                dihedrals.params.phase, dihedrals.params.phase_units, "degrees"
            )
            periodicity = dihedrals.params.periodicity

            for di, (i, j, k, l) in enumerate(dihedrals.indices):
                if isinstance(energy[di], Iterable):
                    dtype = [
                        parmed.topologyobjects.DihedralType(
                            phi_k=energy[di][dj],
                            per=periodicity[di][dj],
                            phase=phase[di][dj],
                            # scee,
                            # scnb,
                            list=pff.dihedral_types,
                        )
                        for dj in range(len(energy[di]))
                    ]
                else:
                    dtype = parmed.topologyobjects.DihedralType(
                        phi_k=energy[di],
                        per=periodicity[di],
                        phase=phase[di],
                        # scee
                        # scnb
                        list=pff.dihedral_types,
                    )
                # assert:
                # dtype.funct = (
                #    9  # hackish: assume all dihedrals are proper and charmm-style
                # )
                pff.dihedrals.append(
                    parmed.topologyobjects.Dihedral(
                        pff.atoms[i],
                        pff.atoms[j],
                        pff.atoms[k],
                        pff.atoms[l],
                        improper=False,
                        type=dtype,
                    )
                )
                pff.dihedral_types.append(dtype)

        return True, TransOutput(proc_input=inputs, data_object=pff)

    def _get_nonbonded(
        self,
        mmff: forcefield.ForceField,
        empty_atom: "parmed.topologyobjects.Atom",
    ) -> Tuple["numpy.ndarray", "numpy.ndarray"]:

        assert (
            mmff.nonbonded.form == "LennardJones"
        ), "Only LJ potential supported for now"

        lj_units = forcefield.nonbonded.potentials.lenjones.LennardJones.get_units()
        scaling_factor = 2 ** (1.0 / 6.0)  # rmin = 2^(1/6) sigma

        rmin = mmff.nonbonded.params.sigma * scaling_factor
        rmin = convert(rmin, lj_units["sigma_units"], empty_atom.urmin.unit.get_name())
        # atom.rmin_14 * rmin_14_factor * scaling_factor,
        epsilon = convert(
            mmff.nonbonded.params.epsilon,
            lj_units["epsilon_units"],
            empty_atom.uepsilon.unit.get_name(),
        )
        # atom.epsilon_14 * epsilon_14_factor,
        return rmin, epsilon


class OpenMMToFFComponent(TransComponent):
    """A component for converting OpenMM ForceField to MMSchema object."""

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

        ff = inputs.data_object
        mm_units = forcefield.ForceField.get_units()

        atoms = ff._atomTypes
        templates = ff._templates
        templates_ = [
            {resname: [atom.name for atom in residue.atoms]}
            for resname, residue in templates.items()
        ]

        # names-defs = [(atom.name, symbols[atom.type]) for atom in residue.atoms for residue in templates.items()]

        for gen in ff.getGenerators():
            if isinstance(gen, simtk.openmm.app.forcefield.NonbondedGenerator):
                nonbond_ff = gen
            elif isinstance(gen, simtk.openmm.app.forcefield.HarmonicBondGenerator):
                bond_ff = gen
            elif isinstance(gen, simtk.openmm.app.forcefield.HarmonicAngleGenerator):
                angle_ff = gen
            elif isinstance(gen, simtk.openmm.app.forcefield.PeriodicTorsionGenerator):
                dihedral_ff = gen
            else:
                raise NotImplementedError

        # Need to map these to potential types
        lj_units = forcefield.nonbonded.potentials.lenjones.LennardJones.get_units()

        data = [
            (
                atom.atomClass,
                atom.element.symbol,
                atom.element.atomic_number,
                atom.mass,
            )
            for _, atom in ff._atomTypes.items()
        ]
        types, symbols, atomic_numbers, masses = zip(*data)
        masses_units = atoms["0"].element.mass.unit.get_name()
        masses = convert(masses, masses_units, mm_units["masses_units"])

        # Non-bonded params
        nb_params = self._get_nonbonded(nonbond_ff)

        lj = forcefield.nonbonded.potentials.LennardJones(
            sigma=nb_params["sigma"], epsilon=nb_params["epsilon"]
        )
        # Need to include sigma_14 and epsilon_14
        nonbonded = forcefield.nonbonded.NonBonded(params=lj, form="LennardJones")

        # charge_groups = None ... should support charge_groups?
        exclusions = None
        inclusions = None

        input_dict = {
            "masses": masses,
            "charges": nb_params["charge"],
            # "bonds": bonds,
            # "angles": angles,
            # "dihedrals": dihedrals,
            "nonbonded": nonbonded,
            "exclusions": exclusions,
            "inclusions": inclusions,
            "defs": types,  # or names?
            "symbols": symbols,
            # "substructs": residues,
            "atomic_numbers": atomic_numbers,
        }

        ff = forcefield.ForceField(**input_dict)
        success = True
        return success, TransOutput(
            proc_input=inputs,
            schema_object=ff,
            success=success,
            schema_version=inputs.schema_version,
            schema_name=inputs.schema_name,
        )

    def _get_nonbonded(self, ff):
        lj_scale = ff.lj14scale
        coul_scale = ff.coulomb14scale

        params = {}
        # How to deal with lj/coul 1-4 scale? relevant to sigma/epsilon scale in parmed?

        for paramName in ff.params.paramNames:
            params[paramName] = [
                ff.params.paramsForType[str(i)][paramName]
                for i in range(len(ff.params.paramsForType))
            ]

        # How do we access units?
        return params

    def _get_bonded(self):
        if bond:
            bonds_units = (
                forcefield.bonded.bonds.potentials.harmonic.Harmonic.get_units()
            )
            bonds_units.update(forcefield.bonded.Bonds.get_units())

            bond_req_factor = convert(
                1.0, bond.type.ureq.unit.get_name(), bonds_units["lengths_units"]
            )
            bond_k_factor = convert(
                1.0, bond.type.uk.unit.get_name(), bonds_units["spring_units"]
            )
            bonds_lengths = [bond.type.req * bond_req_factor for bond in ff.bonds]

            bonds_type = [bond.funct for bond in ff.bonds]
            bonds_k = [bond.type.k * bond_k_factor for bond in ff.bonds]
            connectivity = [
                (bond.atom1.idx, bond.atom2.idx, bond.order or 1) for bond in ff.bonds
            ]

            unique_bonds_type = set(bonds_type)

            if len(unique_bonds_type) > 1:
                raise NotImplementedError("Multiple bond types not yet supported.")
                # params = [
                #    bond_types.get(btype)(spring=bonds_k[bonds_type == btype])
                #    for btype in unique_bonds_type
                #    if bond_types.get(btype)
                # ]
            else:
                bond_funct = unique_bonds_type.pop()
                assert (
                    bond_funct == 1
                ), "Only Harmonic bond potentials supported in mmic_openmm."
                params = bond_types.get(bond_funct)(spring=bonds_k)

            bonds = forcefield.bonded.Bonds(
                params=params,
                lengths=bonds_lengths,
                indices=connectivity,
                form="Harmonic",
            )
        else:
            bonds = None

        if angle:
            angles_units = (
                forcefield.bonded.angles.potentials.harmonic.Harmonic.get_units()
            )
            angles_units.update(forcefield.bonded.Angles.get_units())

            angle_theta_factor = convert(
                1.0, angle.type.utheteq.unit.get_name(), angles_units["angles_units"]
            )
            angle_k_factor = convert(
                1.0, angle.type.uk.unit.get_name(), angles_units["spring_units"]
            )
            angles_ = [angle.type.theteq * angle_theta_factor for angle in ff.angles]
            angles_k = [angle.type.k * angle_k_factor for angle in ff.angles]
            angles_type = [angle.funct for angle in ff.angles]

            angles_indices = [
                (angle.atom1.idx, angle.atom2.idx, angle.atom3.idx)
                for angle in ff.angles
            ]

            unique_angles_type = set(angles_type)
            if len(unique_angles_type) > 1:
                raise NotImplementedError("Multiple angle types not yet supported.")
                # params = [
                #    angleTypes.get(btype)(spring=angles_k[angles_type == btype])
                #    for btype in unique_angles_type
                #    if angleTypes.get(btype)
                # ]
            else:
                angle_funct = unique_angles_type.pop()
                assert (
                    angle_funct == 1
                ), "Only Harmonic angle potentials supported in mmic_openmm."
                params = angle_types.get(angle_funct)(spring=angles_k)

            angles = forcefield.bonded.Angles(
                params=params, angles=angles_, indices=angles_indices, form="Harmonic"
            )
        else:
            angles = None

        if dihedral:
            dihedrals_funct = set([di.funct for di in ff.dihedrals if di.funct != 4])
            # funct == 4: periodic improper dihedrals see https://manual.gromacs.org/documentation/2019/reference-manual/topologies/topology-file-formats.html#tab-topfile2
            dihedrals = [self._get_dihedrals(funct, ff) for funct in dihedrals_funct]

            if len(dihedrals) == 1:  # no point in keeping a list for single-type
                dihedrals = dihedrals.pop()
        else:
            dihedrals = None

        if hasattr(ff, "residues"):
            residues = [(atom.residue.name, atom.residue.idx) for atom in ff.atoms]

    def _get_dihedrals(self, funct, ff):

        assert dihedral_types.get(
            funct
        ), f"Functional form {funct} not supported in mmic_openmm"

        dihedrals_units = dihedral_types[funct].get_units()  # model param units
        dihedrals_units.update(
            forcefield.bonded.Dihedrals.get_units()
        )  # global param units

        dihedral_phi_factor = convert(
            1.0,
            "kcal/mol",  # dihedrals_type.uphi_k.unit.get_name(),
            dihedrals_units["energy_units"],
        )
        dihedral_phase_factor = convert(
            1.0,
            "degrees",  # dihedrals_type.uphase.unit.get_name(),
            dihedrals_units["phase_units"],
        )
        dihedrals_indices = [
            (
                dihedral.atom1.idx,
                dihedral.atom2.idx,
                dihedral.atom3.idx,
                dihedral.atom4.idx,
            )
            for dihedral in ff.dihedrals
        ]

        # Need to look into per, scee, and scnb ... their meaning, units, etc.
        if funct == 9:  # multi-dihedrals ... special case
            phase, energy, per = [], [], []
            for dihedral in ff.dihedrals:
                if dihedral.funct == funct:
                    phase.append(
                        [item.phase * dihedral_phase_factor for item in dihedral.type]
                    )
                    energy.append(
                        [item.phi_k * dihedral_phi_factor for item in dihedral.type]
                    )
                    per.append([item.per for item in dihedral.type])
        else:
            dihedrals_params = [
                [
                    (
                        dihedral.type.phase * dihedral_phase_factor,
                        dihedral.type.phi_k * dihedral_phi_factor,
                        dihedral.type.per,
                        dihedral.type.scee,
                        dihedral.type.scnb,
                    )
                ]
                for dihedral in ff.dihedrals
                if dihedral.funct == funct
            ]
            phase, energy, per, _, _ = zip(*dihedrals_params)

        params = dihedral_types.get(funct)(
            phase=phase,
            energy=energy,
            periodicity=per,
        )

        return forcefield.bonded.Dihedrals(
            params=params,
            indices=dihedrals_indices,
            form=dihedral_types.get(funct).__name__,
        )
