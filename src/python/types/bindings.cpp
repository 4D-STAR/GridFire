#include <pybind11/pybind11.h>
#include <pybind11/stl.h> // Needed for vectors, maps, sets, strings
#include <pybind11/stl_bind.h> // Needed for binding std::vector, std::map etc. if needed directly
#include <format>

#include "bindings.h"

namespace py = pybind11;

#include "gridfire/types/types.h"

void register_type_bindings(const pybind11::module &m) {
    py::class_<gridfire::NetIn>(m, "NetIn")
    .def(py::init<>())
    .def_readwrite("composition", &gridfire::NetIn::composition)
    .def_readwrite("tMax", &gridfire::NetIn::tMax)
    .def_readwrite("dt0", &gridfire::NetIn::dt0)
    .def_readwrite("temperature", &gridfire::NetIn::temperature)
    .def_readwrite("density", &gridfire::NetIn::density)
    .def_readwrite("energy", &gridfire::NetIn::energy)
    .def("__repr__", [](const gridfire::NetIn &netIn) {
        std::stringstream ss;
        ss << "NetIn(composition=" << netIn.composition
           << ", tMax=" << netIn.tMax
           << ", dt0=" << netIn.dt0
           << ", temperature=" << netIn.temperature
           << ", density=" << netIn.density
           << ", energy=" << netIn.energy << ")";
        return ss.str();
    });

    py::class_<gridfire::NetOut>(m, "NetOut")
        .def_readonly("composition", &gridfire::NetOut::composition)
        .def_readonly("num_steps", &gridfire::NetOut::num_steps)
        .def_readonly("energy", &gridfire::NetOut::energy)
        .def_readonly("dEps_dT", &gridfire::NetOut::dEps_dT)
        .def_readonly("dEps_dRho", &gridfire::NetOut::dEps_dRho)
        .def_readonly("specific_neutrino_flux", &gridfire::NetOut::specific_neutrino_flux)
        .def_readonly("specific_neutrino_energy_loss", &gridfire::NetOut::specific_neutrino_energy_loss)
        .def_readonly("eps_nuc", &gridfire::NetOut::eps_nuc, "Instantaneous specific nuclear energy generation rate (erg/g/s) at the final state, net of reaction neutrino losses")
        .def_readonly("dlnEps_dlnT", &gridfire::NetOut::dlnEps_dlnT, "Logarithmic derivative dln(eps_nuc)/dln(T) at the final state (0 if eps_nuc == 0)")
        .def_readonly("dlnEps_dlnRho", &gridfire::NetOut::dlnEps_dlnRho, "Logarithmic derivative dln(eps_nuc)/dln(rho) at the final state (0 if eps_nuc == 0)")
        .def_readonly("neutrino_energy_loss_rate", &gridfire::NetOut::neutrino_energy_loss_rate, "Instantaneous reaction neutrino energy loss rate (erg/g/s) at the final state")
        .def("__repr__", [](const gridfire::NetOut &netOut) {
            std::string repr = std::format(
                "NetOut(<μ> = {} steps = {}, ε = {}, dε/dT = {}, dε/dρ = {}, ε_nuc = {})",
                netOut.composition.getMeanParticleMass(),
                netOut.num_steps,
                netOut.energy,
                netOut.dEps_dT,
                netOut.dEps_dRho,
                netOut.eps_nuc
            );
            return repr;
        });

}
