#pragma once

#include "gridfire/solver/strategies/strategy_abstract.h"
#include "gridfire/engine/engine_abstract.h"
#include "gridfire/types/types.h"
#include "gridfire/exceptions/exceptions.h"
#include "gridfire/config/config.h"

#include "fourdst/atomic/atomicSpecies.h"
#include "fourdst/config/config.h"


#include <functional>
#include <any>
#include <string>
#include <vector>
#include <tuple>
#include <span>
#include <optional>

// SUNDIALS/CVODE headers
#include <cvode/cvode.h>
#include <sundials/sundials_types.h>

// Include headers for linear solvers and N_Vectors
// We will use preprocessor directives to select the correct ones
#include <sundials/sundials_context.h>
#include <sunmatrix/sunmatrix_dense.h>
#include <sunlinsol/sunlinsol_dense.h>

// These are the possible N_Vector implementations. We use the compiler defines to select the most appropriate one for the build.
// If none are defined, we default to serial.

// For precompiled binaries we will need to ensure that we have versions built for all three types (ideally with some runtime
// checks that will fail gracefully if the user tries to use an unsupported type).
#ifdef SUNDIALS_HAVE_OPENMP
    #include <nvector/nvector_openmp.h>
#endif
#ifdef SUNDIALS_HAVE_PTHREADS
    #include <nvector/nvector_pthreads.hh>
#endif
#ifndef SUNDIALS_HAVE_OPENMP
    #ifndef SUNDIALS_HAVE_PTHREADS
        #include <nvector/nvector_serial.h>
    #endif
#endif

namespace gridfire::solver {
        struct PointSolverTimestepContext final : TimestepContextBase {
            const double t;                 ///< Current integration time [s].
            const N_Vector& state;          ///< Current CVODE state vector (N_Vector).
            const double dt;                ///< Last step size [s].
            const double last_step_time;    ///< Time at last callback [s].
            const double T9;                ///< Temperature in GK.
            const double rho;               ///< Density [g cm^-3].
            const size_t num_steps;         ///< Number of CVODE steps taken so far.
            const engine::DynamicEngine& engine;    ///< Reference to the engine.
            const std::vector<fourdst::atomic::Species>& networkSpecies; ///< Species layout.
            const size_t currentConvergenceFailures; ///< Total number of convergence failures
            const size_t currentNonlinearIterations; ///< Total number of non-linear iterations
            const std::map<fourdst::atomic::Species, std::unordered_map<std::string, double>>& reactionContributionMap; ///< Map of reaction contributions for the current step
            engine::scratch::StateBlob& state_ctx; ///< Reference to the engine scratch state blob
            double current_total_energy = 0.0; ///< Current energy generation rate [erg/g/s]
            double current_neutrino_energy_loss_rate = 0.0;   ///< Current neutrino energy loss rate [erg/g/s]

            PointSolverTimestepContext(
                double t,
                const N_Vector& state,
                double dt,
                double last_step_time,
                double t9,
                double rho,
                size_t num_steps,
                const engine::DynamicEngine& engine,
                const std::vector<fourdst::atomic::Species>& networkSpecies,
                size_t currentConvergenceFailure,
                size_t currentNonlinearIterations,
                const std::map<fourdst::atomic::Species, std::unordered_map<std::string, double>> &reactionContributionMap,
                engine::scratch::StateBlob& state_ctx
            );

            [[nodiscard]] std::vector<std::tuple<std::string, std::string>> describe() const override;

            [[nodiscard]] fourdst::composition::Composition getPhysicalComposition() const;

            /**
             * @brief Zero-copy view of the raw solver state vector.
             *
             * Layout: [y_0, ..., y_{N-1}, eps] where y_i are molar abundances in
             * networkSpecies order and eps is the accumulated specific energy.
             * The view is only valid for the duration of the callback.
             */
            [[nodiscard]] std::span<const double> rawState() const;

            /// Abundance of the species at @p species_index (solver state indexing).
            [[nodiscard]] double abundance(size_t species_index) const;

            /// Abundance of @p sp, or std::nullopt if it is not in the active network.
            [[nodiscard]] std::optional<double> abundance(const fourdst::atomic::Species& sp) const;

            /// Accumulated specific energy [erg/g] (the trailing state entry).
            [[nodiscard]] double accumulatedSpecificEnergy() const;
        };

    using TimestepCallback = std::function<void(const PointSolverTimestepContext& context)>; ///< Type alias for a timestep callback function.

    struct PointSolverContext final : SolverContextBase {
        SUNContext sun_ctx = nullptr;   ///< SUNDIALS context (lifetime of the solver).
        void* cvode_mem = nullptr;      ///< CVODE memory block.
        N_Vector Y = nullptr;           ///< CVODE state vector (species + energy accumulator).
        N_Vector YErr = nullptr;        ///< Estimated local errors.
        SUNMatrix J = nullptr;          ///< Dense Jacobian matrix.
        SUNLinearSolver LS = nullptr;   ///< Dense linear solver.

        std::unique_ptr<engine::scratch::StateBlob> engine_ctx;


        std::optional<TimestepCallback> callback;      ///< Optional per-step callback.
        int num_steps = 0;              ///< CVODE step counter (used for diagnostics and triggers).

        bool stdout_logging = true; ///< If true, print per-step logs and use CV_ONE_STEP.

        N_Vector constraints = nullptr; ///< CVODE constraints vector (>= 0 for species entries).

        std::optional<double> abs_tol;        ///< User-specified absolute tolerance.
        std::optional<double> rel_tol;        ///< User-specified relative tolerance.

        bool detailed_step_logging = false;    ///< If true, log detailed step diagnostics (error ratios, Jacobian, species balance).

        size_t last_size = 0;
        size_t last_composition_hash = 0ULL;
        sunrealtype last_good_time_step = 0ULL;

        void init() override;
        void set_stdout_logging(bool enable) override;
        void set_detailed_logging(bool enable) override;

        void reset_all();
        void reset_user();
        void reset_cvode();
        void clear_context();
        void init_context();

        [[nodiscard]] bool has_context() const;

        explicit PointSolverContext(const engine::scratch::StateBlob& engine_ctx);
        ~PointSolverContext() override;

    };

    /**
     * @class PointSolver
     * @brief Stiff ODE integrator backed by SUNDIALS CVODE (BDF) for network + energy.
     *
     * Integrates the nuclear network abundances along with a final accumulator entry for specific
     * energy using CVODE's BDF method and a dense linear solver. The state vector layout is:
     *   [y_0, y_1, ..., y_{N-1}, eps], where eps is the accumulated specific energy (erg/g).
     *
     * Implementation summary:
     *  - Creates a SUNContext and CVODE memory; initializes the state from a Composition.
     *  - Enforces non-negativity on species via CVodeSetConstraints (>= 0 for all species slots).
     *  - Uses a user-provided DynamicEngine to compute RHS and to fill the dense Jacobian.
     *  - The Jacobian is assembled column-major into a SUNDenseMatrix; the energy row/column is
     *    currently set to zero (decoupled from abundances in the linearization).
     *  - An internal trigger can rebuild the engine/network; when triggered, CVODE resources are
     *    torn down and recreated with the new network size, preserving the energy accumulator.
     *  - The CVODE RHS wrapper captures exceptions::StaleEngineTrigger from the engine evaluation
     *    path as recoverable (return code 1) and stores a copy in user-data for the driver loop.
     *
     * @par Example
     * @code
     * using gridfire::solver::CVODESolverStrategy;
     * using gridfire::solver::NetIn;
     *
     * CVODESolverStrategy solver(engine);
     * NetIn in;
     * in.temperature = 1.0e9; // K
     * in.density = 1.0e6;     // g/cm^3
     * in.tMax = 1.0;          // s
     * in.composition = initialComposition;
     * auto out = solver.evaluate(in);
     * std::cout << "Final energy: " << out.energy << " erg/g\n";
     * @endcode
     */
    class PointSolver final : public SingleZoneDynamicNetworkSolver {
    public:
        /**
         * @brief Construct the CVODE strategy and create a SUNDIALS context.
         * @param engine DynamicEngine used for RHS/Jacobian evaluation and network access.
         * @throws std::runtime_error If SUNContext_Create fails.
         */
        explicit PointSolver(
            const engine::DynamicEngine& engine
        );

        PointSolver(
            const engine::DynamicEngine& engine,
            const config::GridFireConfig& config
        );


        config::GridFireConfig getConfig() const { return *m_config; }
        /**
         * @brief Integrate from t=0 to netIn.tMax and return final composition and energy.
         *
         * Implementation summary:
         *  - Converts temperature to T9, initializes CVODE memory and state (size = numSpecies + 1).
         *  - Repeatedly calls CVode in single-step or normal mode depending on stdout logging.
         *  - Wraps RHS to capture exceptions::StaleEngineTrigger as a recoverable step failure;
         *    if present after a step, it is rethrown for upstream handling.
         *  - Prints/collects diagnostics per step (step size, energy, solver iterations).
         *  - On trigger activation, rebuilds CVODE resources to reflect a changed network and
         *    reinitialized the state using the latest engine composition, preserving energy.
         *  - At the end, converts molar abundances to mass fractions and assembles NetOut,
         *    including derivatives of energy w.r.t. T and rho from the engine.
         *
         * @param solver_ctx
         * @param netIn Inputs: temperature [K], density [g cm^-3], tMax [s], composition.
         * @return NetOut containing final Composition, accumulated energy [erg/g], step count,
         *         and dEps/dT, dEps/dRho.
         * @throws std::runtime_error If any CVODE or SUNDIALS call fails (negative return codes),
         *         or if internal consistency checks fail during engine updates.
         * @throws exceptions::StaleEngineTrigger Propagated if the engine signals a stale state
         *         during RHS evaluation (captured in the wrapper then rethrown here).
         */
        NetOut evaluate(
            SolverContextBase& solver_ctx,
            const NetIn& netIn
        ) const override;

        /**
         * @brief Call to evaluate which will let the user control if the trigger reasoning is displayed
         * @param solver_ctx
         * @param netIn Inputs: temperature [K], density [g cm^-3], tMax [s], composition.
         * @param displayTrigger Boolean flag to control if trigger reasoning is displayed
         * @param forceReinitialize Boolean flag to force reinitialization of CVODE resources at the start
         * @return NetOut containing final Composition, accumulated energy [erg/g], step count,
         *         and dEps/dT, dEps/dRho.
         * @throws std::runtime_error If any CVODE or SUNDIALS call fails (negative return codes),
         *         or if internal consistency checks fail during engine updates.
         * @throws exceptions::StaleEngineTrigger Propagated if the engine signals a stale state
         *         during RHS evaluation (captured in the wrapper then rethrown here).
         */
        NetOut evaluate(
            SolverContextBase& solver_ctx,
            const NetIn& netIn,
            bool displayTrigger,
            bool forceReinitialize = false
        ) const;

    private:
        /**
         * @struct CVODEUserData
         * @brief A helper struct to pass C++ context to C-style CVODE callbacks.
         *
         * Carries pointers back to the solver instance and engine, the current thermodynamic
         * state, energy accumulator, and a slot to capture a copy of exceptions::StaleEngineTrigger
         * from RHS evaluation. The RHS wrapper treats this as a recoverable failure and returns 1
         * to CVODE, then the driver loop inspects and rethrows.
         */
        struct CVODEUserData {
            const PointSolver* solver_instance{}; // Pointer back to the class instance
            PointSolverContext* sctx;    // Pointer to the solver context
            engine::scratch::StateBlob& ctx;
            const engine::DynamicEngine* engine{};
            double abs_tol{};
            double rel_tol{};
            double T9{};
            double rho{};
            double energy_generation_rate{};
            double last_step_size{};
            double previous_accepted_energy{};
            double last_energy_change{};
            std::vector<double> previous_accepted_abundances;
            std::vector<double> last_mass_fraction_changes;
            const std::vector<fourdst::atomic::Species>* networkSpecies{};
            std::unique_ptr<exceptions::EngineError> captured_exception = nullptr;
            std::optional<std::map<fourdst::atomic::Species, std::unordered_map<std::string, double>>> reaction_contribution_map;
            double neutrino_energy_loss_rate = 0.0;
            double total_neutrino_flux = 0.0;
        };

        struct CVODERHSOutputData {
            std::optional<std::map<fourdst::atomic::Species, std::unordered_map<std::string, double>>> reaction_contribution_map;
            double neutrino_energy_loss_rate;
            double total_neutrino_flux;
        };

    private:
        fourdst::config::Config<config::GridFireConfig> m_config;
        quill::Logger* m_logger = fourdst::logging::LogManager::getInstance().getLogger("log");
        /**
         * @brief CVODE RHS C-wrapper that delegates to calculate_rhs and captures exceptions.
         * @return 0 on success; 1 on recoverable StaleEngineTrigger; -1 on other failures.
         */
        static int cvode_rhs_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, void *user_data);

        /**
         * @brief Builds CVODE error weights from runtime species and energy scales.
         */
        static int cvode_error_weight_wrapper(N_Vector y, N_Vector ewt, void* user_data);
        /**
         * @brief CVODE dense Jacobian C-wrapper that fills SUNDenseMatrix using the engine.
         *
         * Assembles J(i,j) = d(f_i)/d(y_j) for all species using engine->getJacobianMatrixEntry,
         * then zeros the last row and column corresponding to the energy variable.
         */
        static int cvode_jac_wrapper(sunrealtype t, N_Vector y, N_Vector ydot, SUNMatrix J, void *user_data, N_Vector tmp1, N_Vector tmp2, N_Vector tmp3);

        /**
         * @brief CVODE error handler that logs errors and warnings from SUNDIALS using the solver's logger.
         * @param line
         * @param func
         * @param file
         * @param msg
         * @param err_code
         * @param err_user_data
         * @param sunctx
         */
        static void cvode_error_handler(int line, const char *func, const char *file, const char *msg, SUNErrCode err_code, void *err_user_data, SUNContext sunctx);
        /**
         * @brief Compute RHS into ydot at time t from the engine and current state y.
         *
         * Converts the CVODE state to a Composition (mass fractions) and calls
         * engine.calculateRHSAndEnergy(T9, rho). Negative small abundances are clamped to zero
         * before constructing Composition. On stale engine, throws exceptions::StaleEngineTrigger.
         */
        CVODERHSOutputData calculate_rhs(sunrealtype t, N_Vector y, N_Vector ydot, const CVODEUserData *data) const;

        /**
         * @brief Allocate and initialize CVODE vectors, linear algebra, tolerances, and constraints.
         *
         * State vector m_Y is sized to N (numSpecies + 1). Species slots are initialized from Composition
         * molar abundances when present, otherwise a tiny positive value; the last slot is set to
         * accumulatedEnergy. Sets runtime-derived error weights, non-negativity constraints for species, maximum
         * step size, creates a dense matrix and dense linear solver, and registers the Jacobian.
         */
        void initialize_cvode_integration_resources(
            PointSolverContext* ctx,
            uint64_t N,
            size_t numSpecies,
            double current_time,
            const fourdst::composition::Composition& composition,
            double absTol,
            double relTol,
            double accumulatedEnergy
        ) const;


        /**
         * @brief Compute and print per-component error ratios; run diagnostic helpers.
         *
         * Gathers CVODE's estimated local errors, converts the state to a Composition, and prints a
         * sorted table of species with the highest error ratios; then invokes diagnostic routines to
         * inspect Jacobian stiffness and species balance.
         */
        void log_step_diagnostics(
            PointSolverContext* sctx_p,
            engine::scratch::StateBlob &ctx,
            const CVODEUserData& user_data,
            bool displayJacobianStiffness,
            bool displaySpeciesBalance,
            bool to_file, std::optional<std::string> filename
        ) const;
    };
}
