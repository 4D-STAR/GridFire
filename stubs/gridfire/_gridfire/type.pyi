"""
GridFire type bindings
"""
from __future__ import annotations
import fourdst._phys.composition
import typing
__all__: list[str] = ['NetIn', 'NetOut']
class NetIn:
    composition: fourdst._phys.composition.Composition
    def __init__(self) -> None:
        ...
    def __repr__(self) -> str:
        ...
    @property
    def density(self) -> float:
        ...
    @density.setter
    def density(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def dt0(self) -> float:
        ...
    @dt0.setter
    def dt0(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def energy(self) -> float:
        ...
    @energy.setter
    def energy(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def tMax(self) -> float:
        ...
    @tMax.setter
    def tMax(self, arg0: typing.SupportsFloat) -> None:
        ...
    @property
    def temperature(self) -> float:
        ...
    @temperature.setter
    def temperature(self, arg0: typing.SupportsFloat) -> None:
        ...
class NetOut:
    def __repr__(self) -> str:
        ...
    @property
    def composition(self) -> fourdst._phys.composition.Composition:
        ...
    @property
    def dEps_dRho(self) -> float:
        ...
    @property
    def dEps_dT(self) -> float:
        ...
    @property
    def dlnEps_dlnRho(self) -> float:
        """
        Logarithmic derivative dln(eps_nuc)/dln(rho) at the final state (0 if eps_nuc == 0)
        """
    @property
    def dlnEps_dlnT(self) -> float:
        """
        Logarithmic derivative dln(eps_nuc)/dln(T) at the final state (0 if eps_nuc == 0)
        """
    @property
    def energy(self) -> float:
        ...
    @property
    def eps_nuc(self) -> float:
        """
        Instantaneous specific nuclear energy generation rate (erg/g/s) at the final state, net of reaction neutrino losses
        """
    @property
    def neutrino_energy_loss_rate(self) -> float:
        """
        Instantaneous reaction neutrino energy loss rate (erg/g/s) at the final state
        """
    @property
    def num_steps(self) -> int:
        ...
    @property
    def specific_neutrino_energy_loss(self) -> float:
        ...
    @property
    def specific_neutrino_flux(self) -> float:
        ...
