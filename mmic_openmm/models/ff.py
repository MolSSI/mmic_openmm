from typing import Dict, Any, Optional, Union, List
from mmic_translator.models.base import ToolkitModel
from mmelemental.models.forcefield.mm_ff import ForceField

# Import Components
from mmic_openmm.components.ff_component import FFToOpenMMComponent
from mmic_openmm.components.ff_component import OpenMMToFFComponent

# OpenMM library
from simtk.unit.quantity import Quantity
from simtk.unit.unit import Unit, BaseUnit, BaseDimension
from simtk.openmm.vec3 import Vec3
from simtk.openmm import app
from simtk.openmm import __version__ as openmm_version
from simtk.openmm.app import ForceField as OpenMMForceField

__all__ = ["OpenMMFF"]


class OpenMMFF(ToolkitModel):
    """A model for storing a ForceField object in OpenMM."""

    @classmethod
    def engine(cls):
        return "openmm", openmm_version

    @classmethod
    def dtype(cls):
        """Returns the fundamental molecule object type."""
        return OpenMMForceField

    @classmethod
    def isvalid(cls, data):
        """Makes sure the Structure object stores atoms."""
        if hasattr(data, "createSystem"):
            if callable(data.createSystem):
                return data
        raise ValueError("OpenMM Forcefield object is invalid. Cannot create a system!")

    @classmethod
    def from_file(cls, filename: Union[str, List[str]], **kwargs) -> "OpenMMFF":
        """
        Constructs an ParmedFF object from file(s).

        Parameters
        ----------
        filename : str or List[str]
            The forcefield filename(s) to read
        **kwargs
            Any additional keywords to pass to the constructor
        Returns
        -------
        OpenMMFF
            A constructed OpenMMFF object.
        """
        if isinstance(filename, str):
            ff = OpenMMForceField(filename, **kwargs)
        else:
            ff = OpenMMForceField(*filename, **kwargs)

        return cls(data=ff)

    @classmethod
    def from_schema(
        cls, data: ForceField, version: Optional[int] = None, **kwargs: Dict[str, Any]
    ) -> "ParmedFF":
        """
        Constructs an ParmedFF object from an MMSchema ForceField object.
        Parameters
        ----------
        data: ForceField
            Data to construct the forcefield object from.
        version: int, optional
            Schema version e.g. 1. Overwrites data.schema_version.
        **kwargs
            Additional kwargs to pass to the constructors.
        Returns
        -------
        OpenMMFF
            A constructed OpenMMFF object.
        """
        inputs = {
            "schema_object": data,
            "schema_version": version or data.schema_version,
        }
        out = FFToOpenMMComponent.compute(inputs)
        return cls(data=out.data_object, units=out.data_units)

    def to_file(self, filename: str, dtype: str = None, **kwargs):
        """Writes the forcefield to a file.
        Parameters
        ----------
        filename : str
            The filename to write to
        dtype : Optional[str], optional
            File format
        **kwargs
            Additional kwargs to pass to the constructors.
        """
        if dtype:
            kwargs["format"] = dtype
        self.data.save(filename, **kwargs)

    def to_schema(self, version: Optional[int] = 0, **kwargs) -> ForceField:
        """Converts the forcefield to MMSchema ForceField object.
        Parameters
        ----------
        version: int, optional
            Schema specification version to comply with e.g. 1.
        **kwargs
            Additional kwargs to pass to the constructor.
        """
        inputs = {"data_object": self.data, "schema_version": version, "kwargs": kwargs}
        out = ParmedToFFComponent.compute(inputs)
        if version:
            assert version == out.schema_version
        return out.schema_object
