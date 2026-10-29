#ifndef GF_GRIDFIRE_CONTEXT_H
#define GF_GRIDFIRE_CONTEXT_H

#include "gridfire/gridfire.h"
#include "fourdst/atomic/atomicSpecies.h"

#include <memory>
#include <vector>

/**
 * @brief Optional destinations for the instantaneous final-state quantities in NetOut.
 *
 * For a point context each pointer refers to a single double; for a grid context each refers to an array with one
 * entry per zone. A null member (or a null struct) means that quantity is not written.
 */
struct GFInstantaneousOutputs {
    double* eps_nuc = nullptr;
    double* dlnEps_dlnT = nullptr;
    double* dlnEps_dlnRho = nullptr;
    double* neutrino_energy_loss_rate = nullptr;

    void write(const gridfire::NetOut& netOut, size_t index) const;
};

enum class GFContextType {
    POINT,
    GRID
};

struct GFContext {
    virtual ~GFContext() = default;

    std::unique_ptr<gridfire::policy::NetworkPolicy> policy;
    const gridfire::engine::DynamicEngine* engine;
    std::unique_ptr<gridfire::engine::scratch::StateBlob> engine_ctx;
    std::vector<fourdst::atomic::Species> speciesList;

    virtual void init_species_map(const std::vector<std::string>& species_names);
    virtual void init_engine_from_policy(const std::string& policy_name, const double *abundances, size_t num_species);
    virtual void init_solver_from_engine() = 0;

    fourdst::composition::Composition init_composition_from_abundance_vector(const std::vector<double> &abundances, size_t num_species) const;
    std::string last_error;
};

struct GFPointContext final: GFContext{
    std::unique_ptr<gridfire::solver::SingleZoneDynamicNetworkSolver> solver;
    std::unique_ptr<gridfire::solver::SolverContextBase> solver_ctx;

    void init_solver_from_engine() override;

    int evolve(
        const double* Y_in,
        size_t num_species,
        double T,
        double rho,
        double tMax,
        double dt0,
        double* Y_out,
        double& energy_out,
        double& dEps_dT,
        double& dEps_dRho,
        double& specific_neutrino_energy_loss,
        double& specific_neutrino_flux,
        double& mass_lost,
        const GFInstantaneousOutputs* instantaneous = nullptr
    ) const;

};

struct GFGridContext final : GFContext {
    std::unique_ptr<gridfire::solver::SingleZoneDynamicNetworkSolver> local_solver;
    std::unique_ptr<gridfire::solver::MultiZoneDynamicNetworkSolver> solver;
    std::unique_ptr<gridfire::solver::SolverContextBase> solver_ctx;

    void init_solver_from_engine() override;

    size_t zones;

    void set_zones(const size_t num_zones) {
        zones = num_zones;
    }

    [[nodiscard]] size_t get_zones() const {
        return zones;
    }

    int evolve(
        const double* Y_in,
        size_t num_species,
        const double* T,
        const double* rho,
        double tMax,
        double dt0,
        double* Y_out,
        double* energy_out,
        double* dEps_dT,
        double* dEps_dRho,
        double* specific_neutrino_energy_loss,
        double* specific_neutrino_flux,
        double* mass_lost,
        const GFInstantaneousOutputs* instantaneous = nullptr
    ) const;

};

std::unique_ptr<GFContext> make_gf_context(const GFContextType& type);

#endif