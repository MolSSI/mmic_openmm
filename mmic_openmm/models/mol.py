from typing import Dict, Any, Optional
from mmic_translator.models.base import ToolkitModel
from mmelemental.models import Molecule
from mmelemental.types import Array
import numpy
from pydantic import Field
from pathlib import Path

# OpenMM converter components
from mmic_openmm.components.mol_component import OpenMMToMolComponent
from mmic_openmm.components.mol_component import MolToOpenMMComponent

# OpenMM library
from simtk.unit.quantity import Quantity
from simtk.unit.unit import Unit, BaseUnit, BaseDimension
from simtk.openmm.vec3 import Vec3
from simtk.openmm import app
from simtk.openmm import __version__ as openmm_version

__all__ = ["OpenMMMol"]


class OpenMMUnit(Unit):
    @classmethod
    def __get_validators__(cls):
        yield cls.validate

    @classmethod
    def validate(cls, v):
        try:
            v.get_name()
            v.get_symbol()
        except AttributeError as e:
            raise AttributeError(f"{v} is not a valid OpenMM Unit object: {e}")
        return v


class OpenMMMol(ToolkitModel):
    """A model for OpenMM molecule class storing a Topology object."""

    positions: Array[float] = Field(
        ...,
        description="Particle (e.g. atomic) positions numpy array of length (natoms*3,).",
    )
    positions_units: OpenMMUnit = Field(
        ...,
        description="Unit object for positions e.g. Unit(base_dim=BaseDimension('length'), name='nanometer', symbol='nm')",
    )

    @classmethod
    def engine(cls):
        return "openmm", openmm_version

    @classmethod
    def dtype(cls):
        """Returns the fundamental molecule object type."""
        return app.topology.Topology

    @classmethod
    def isvalid(cls, data):
        """Makes sure the Structure object stores atoms."""
        if data.getNumAtoms() > 0:
            return data
        raise ValueError("OpenMM Topology object does not contain any atoms!")

    @classmethod
    def from_file(
        cls, filename: str, top_filename: str = None, **kwargs
    ) -> "OpenMMMol":
        """
        Constructs an OpenMMMol object from file(s).

        Parameters
        ----------
        filename : str
            The atomic positions filename to read
        top_filename: str, optional
            The topology filename to read
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        OpenMMMol
            A constructed OpenMMMol class.
        """

        ext = Path(filename).suffix
        top = None

        # Add support for PSF files
        if ext == ".pdb":
            fileobj = app.PDBFile(filename, **kwargs)
            top = fileobj.topology
        elif ext == ".gro":
            fileobj = app.GromacsGroFile(filename, **kwargs)
            if top_filename:
                gro = app.GromacsTopFile(
                    top_filename,
                    periodicBoxVectors=fileobj.getPeriodicBoxVectors(),
                    includeDir="/usr/share/gromacs/top",
                )
                top = gro.topology
        elif ext == ".inpcrd":
            fileobj = app.AmberInpcrdFile(filename, **kwargs)
            if top_filename:
                top = app.AmberPrmtopFile(top_filename)
        else:
            raise NotImplementedError(
                "Only PDB (.pdb), GROMACS (.gro), and AMBER (.inpcrd) coords files supported by mmic_pdb."
            )

        if top is None:
            top = app.topology.Topology()
            top._numAtoms = (
                len(fileobj.atomNames)
                if hasattr(fileobj, "atomNames")
                else len(fileobj.positions)
            )

        positions = numpy.array(
            [(pos.x, pos.y, pos.z) for pos in fileobj.positions]
        ).flatten()

        return cls(
            data=top,
            positions=positions,
            positions_units=fileobj.positions.unit,
        )

    @classmethod
    def from_schema(
        cls, data: Molecule, version: Optional[int] = None, **kwargs: Dict[str, Any]
    ) -> "OpenMMMol":
        """
        Constructs an OpenMMMol object from an MMSchema Molecule object.
        Parameters
        ----------
        data: Molecule
            Data to construct Molecule from.
        version: int, optional
            Schema version e.g. 1. Overwrites data.schema_version.
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        OpenMMMol
            A constructed OpenMMMol class.
        """
        inputs = {
            "schema_object": data,
            "schema_version": version or data.schema_version,
        }
        out = MolToOpenMMComponent.compute(inputs)
        return cls(data=out.data_object, units=out.data_units)

    def to_file(self, filename: str, dtype: str = None, mode: str = "w", **kwargs):
        """Writes the molecule to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        **kwargs
            Additional kwargs to pass to the constructors. kwargs takes precedence over  data.
        """
        if dtype is None:
            ext = Path(filename).suffix
        else:
            ext = "." + dtype

        positions = Quantity(
            [
                Vec3(x=x, y=y, z=z)
                for x, y, z in self.positions.reshape(self.data._numAtoms, 3)
            ],
            unit=self.positions_units,
        )  # assume dim=3 always

        if ext == ".pdb":
            writeFile = app.PDBFile.writeFile
        else:
            raise NotImplementedError(
                "mmic_openmm ssupports writing only PDB files for now."
            )

        writeFile(self.data, positions, open(filename, mode=mode), **kwargs)

    def to_schema(self, version: Optional[int] = 0, **kwargs) -> Molecule:
        """Converts the molecule to MMSchema molecule.
        Parameters
        ----------
        version: str, optional
            Schema specification version to comply with e.g. 1.0.1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        kwargs["positions"] = kwargs.get("positions", self.positions)
        kwargs["positions_units"] = kwargs.get("positions_units", self.positions_units)

        inputs = {"data_object": self.data, "schema_version": version, "kwargs": kwargs}
        out = OpenMMToMolComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
