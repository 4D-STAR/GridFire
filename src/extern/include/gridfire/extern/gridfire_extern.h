#ifndef GF_GRIDFIRE_EXTERN_H
#define GF_GRIDFIRE_EXTERN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif
    enum GF_TYPE {
        SINGLE_ZONE = 1001,
        MULTI_ZONE = 1002
    };


    enum FDSSE_ERROR_CODES {
        FDSSE_NON_4DSTAR_ERROR = -102,
        FDSSE_UNKNOWN_ERROR = -101,
        FDSSE_SUCCESS = 1,
        FDSSE_UNKNOWN_SYMBOL_ERROR = 100,
        FDSSE_SPECIES_ERROR = 101,
        FDSSE_INVALID_COMPOSITION_ERROR = 102,
        FDSSE_COMPOSITION_ERROR = 103
    };

    enum GF_ERROR_CODES {
        GF_NON_GRIDFIRE_ERROR = -2,
        GF_UNKNOWN_ERROR = -1,
        GF_SUCCESS = 0,

        GF_INVALID_QSE_SOLUTION_ERROR = 5,
        GF_FAILED_TO_PARTITION_ENGINE_ERROR = 6,
        GF_NETWORK_RESIZED_ERROR = 7,
        GF_UNABLE_TO_SET_NETWORK_REACTIONS_ERROR = 8,
        GF_BAD_COLLECTION_ERROR = 9,
        GF_BAD_RHS_ENGINE_ERROR = 10,
        GF_STALE_JACOBIAN_ERROR = 11,
        GF_UNINITIALIZED_JACOBIAN_ERROR = 12,
        GF_UNKNOWN_JACOBIAN_ERROR = 13,
        GF_JACOBIAN_ERROR = 14,
        GF_ENGINE_ERROR = 15,

        GF_MISSING_BASE_REACTION_ERROR = 16,
        GF_MISSING_SEED_SPECIES_ERROR = 17,
        GF_MISSING_KEY_REACTION_ERROR = 18,
        GF_POLICY_ERROR = 19,

        GF_REACTION_PARSING_ERROR = 20,
        GF_REACTION_ERROR = 21,

        GF_SINGULAR_JACOBIAN_ERROR = 22,
        GF_ILL_CONDITIONED_JACOBIAN_ERROR = 23,
        GF_CVODE_SOLVER_FAILURE_ERROR = 24,
        GF_KINSOL_SOLVER_FAILURE_ERROR = 25,
        GF_SUNDIALS_ERROR = 26,
        GF_SOLVER_ERROR = 27,

        GF_HASHING_ERROR = 28,
        GF_UTILITY_ERROR = 29,

        GF_DEBUG_ERROR = 30,

        GF_GRIDFIRE_ERROR = 31,
        GF_UNINITIALIZED_INPUT_MEMORY_ERROR = 32,
        GF_UNINITIALIZED_OUTPUT_MEMORY_ERROR = 33,

        GF_INVALID_NUM_SPECIES = 34,

        GF_INVALID_TIMESTEPS = 35,
        GF_UNKNOWN_FREE_TYPE = 36,

        GF_INVALID_TYPE = 37,
    };

    char* gf_get_last_error_message(void* ptr);

    char* gf_error_code_to_string(int error_code);

    void* gf_init(const enum GF_TYPE type);

    int gf_free(const enum GF_TYPE type, void *ctx);

    int gf_set_num_zones(const enum GF_TYPE type, void* ptr, const size_t num_zones);

    int gf_register_species(void* ptr, const int num_species, const char** species_names);

    int gf_construct_engine_from_policy(void* ptr, const char* policy_name, const double *abundances, size_t num_species);

    int gf_construct_solver_from_engine(void* ptr);

    int gf_evolve(
        enum GF_TYPE type,
        void* ptr,
        const void* Y_in,
        size_t num_species,
        const void* T,
        const void* rho,
        double tMax,
        double dt0,
        void* Y_out,
        void* energy_out,
        void* dEps_dT,
        void* dEps_dRho,
        void* specific_neutrino_energy_loss,
        void* specific_neutrino_flux,
        void* mass_lost
    );

    /**
     * Same as gf_evolve, but additionally returns the instantaneous quantities at the final state. These are evaluated
     * at the same composition, temperature, and density as dEps_dT and dEps_dRho, so the four are mutually consistent:
     *   eps_nuc                   - specific nuclear energy generation rate net of reaction neutrino losses [erg/g/s]
     *   dlnEps_dlnT               - (T / eps_nuc) * dEps_dT   (0 if eps_nuc == 0)
     *   dlnEps_dlnRho             - (rho / eps_nuc) * dEps_dRho (0 if eps_nuc == 0)
     *   neutrino_energy_loss_rate - reaction neutrino energy loss rate [erg/g/s]
     * For MULTI_ZONE each of these is an array with one entry per zone.
     */
    int gf_evolve_ext(
        enum GF_TYPE type,
        void* ptr,
        const void* Y_in,
        size_t num_species,
        const void* T,
        const void* rho,
        double tMax,
        double dt0,
        void* Y_out,
        void* energy_out,
        void* dEps_dT,
        void* dEps_dRho,
        void* specific_neutrino_energy_loss,
        void* specific_neutrino_flux,
        void* mass_lost,
        void* eps_nuc,
        void* dlnEps_dlnT,
        void* dlnEps_dlnRho,
        void* neutrino_energy_loss_rate
    );

#ifdef __cplusplus
}
#endif

#endif
