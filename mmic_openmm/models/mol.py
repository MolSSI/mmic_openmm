from typing import Dict, Any, Optional
from mmic_translator.models.base import ToolkitModel
from mmelemental.models import Molecule
from mmelemental.types import Array
from simtk.openmm import app
import numpy
from pydantic import Field
from pathlib import Path

# OpenMM converter components
from mmic_openmm.components.mol_component import OpenMMToMolComponent
from mmic_openmm.components.mol_component import MolToOpenMMComponent

__all__ = ["OpenMMMol"]


class OpenMMMol(ToolkitModel):
    """A model for OpenMM molecule class storing a Topology object."""

    positions: Array[float] = Field(
        ...,
        description="Particle (e.g. atomic) positions numpy array of length natomsx3.",
    )
    positions_units: str = Field(..., description="Unit for positions e.g. nanometer.")

    @property
    def dtype(self):
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

        if ext == ".pdb":
            fileobj = app.PDBFile(filename, **kwargs)
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
            top._numAtoms = len(fileobj.atomNames)

        positions = numpy.array(
            [(pos.x, pos.y, pos.z) for pos in fileobj.positions]
        ).T.flatten()

        return cls(
            data=top,
            positions=positions,
            positions_units=fileobj.positions.unit.get_name(),
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
        if dtype:
            kwargs["format"] = kwargs.get("format", dtype)
        if mode == "w":
            kwargs["overwrite"] = True
        elif mode == "a":
            kwargs["overwrite"] = False
        else:
            raise NotImplementedError(
                "File write mode can be either 'w' (write) or 'a' (append) for now."
            )

        self.data.save(filename, **kwargs)

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
