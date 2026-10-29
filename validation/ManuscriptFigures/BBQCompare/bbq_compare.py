import numpy as np
import pandas as pd
from IPython.core.pylabtools import figsize
from gridfire.solver import PointSolver, PointSolverContext
from gridfire.policy import MainSequencePolicy

from gridfire.engine import GraphEngine, MultiscalePartitioningEngineView, AdaptiveEngineView
from gridfire.engine import NetworkBuildDepth
from fourdst.composition.utils import buildCompositionFromMassFractions

from scipy.signal import find_peaks

from gridfire.config import GridFireConfig

from fourdst.composition import Composition
from scipy.integrate import trapezoid

from fourdst.composition import CanonicalComposition
from fourdst.atomic import Species
from gridfire.type import NetIn, NetOut

import matplotlib.pyplot as plt


## Note that my default style uses tex rendering. If you do not have tex installed
## simply comment out this line
plt.style.use("../utils/pub.mplstyle")

from scipy.interpolate import interp1d, CubicSpline

from enum import Enum

import sys
import os
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), "../utils")))

from logger import StepLogger

class ShowSave(Enum):
    SHOW="SHOW"
    SAVE="SAVE"

    def __str__(self):
        return self.value

def rescale_composition(comp_ref : Composition, ZZs : float, Y_primordial : float = 0.248) -> Composition:
    CC : CanonicalComposition  = comp_ref.getCanonicalComposition()

    dY_dZ = (CC.Y - Y_primordial) / CC.Z

    Z_new = CC.Z * (10**ZZs)
    Y_bulk_new = Y_primordial + (dY_dZ * Z_new)
    X_new = 1.0 - Z_new - Y_bulk_new

    if X_new < 0: raise ValueError(f"ZZs={ZZs} yields unphysical composition (X < 0)")

    ratio_H  = X_new / CC.X if CC.X > 0 else 0
    ratio_He = Y_bulk_new / CC.Y if CC.Y > 0 else 0
    ratio_Z  = Z_new / CC.Z if CC.Z > 0 else 0

    Y_new_list = []
    newComp : Composition = Composition()
    s: Species
    for s  in comp_ref.getRegisteredSpecies():
        Xi_ref = comp_ref.getMassFraction(s)

        if s.el() == "H":
            Xi_new = Xi_ref * ratio_H
        elif s.el() == "He":
            Xi_new = Xi_ref * ratio_He
        else:
            Xi_new = Xi_ref * ratio_Z

        Y = Xi_new / s.mass()
        newComp.registerSpecies(s)
        newComp.setMolarAbundance(s, Y)

    return newComp

def init_composition(ZZs : float = 0) -> Composition:
    X_GS98 = [0.73395, 0.00005, 0.2490, 0.00281, 0.00101, 0.00883, 0.00149, 0.00064, 0.00066, 0.00035, 0.00008, 0.00006, 0.00107]
    S_GS98 = ["H-1",   "He-3", "He-4", "C-12",  "N-14",  "O-16",  "Ne-20", "Mg-24", "Si-28", "S-32",  "Ar-36", "Ca-40", "Fe-56"]
    return buildCompositionFromMassFractions(S_GS98, X_GS98)
def init_netIn(temp: float, rho: float, time: float, comp: Composition) -> NetIn:
    n : NetIn = NetIn()
    n.temperature = temp
    n.density = rho
    n.tMax = time
    n.dt0 = 1e-12
    n.composition = comp
    return n

def years_to_seconds(years: float) -> float:
    return years * 3.1536e7

def quantify_engine_error(df_base, df_approx, r_base: NetOut, r_approx: NetOut, species_list, floor_val=1e-30):
    temporal_results = {}
    final_state_results = {}

    t_base = df_base['t'].values

    tracking_cols = ['eps'] + species_list

    for col in tracking_cols:
        if col not in df_base.columns or col not in df_approx.columns:
            continue

        y_base = df_base[col].values

        interpolator = interp1d(
            df_approx['t'],
            df_approx[col],
            kind='linear',
            bounds_error=False,
            fill_value=(df_approx[col].iloc[0], df_approx[col].iloc[-1])
        )

        y_approx_interp = interpolator(t_base)

        abs_diff = np.abs(y_approx_interp - y_base)

        rel_diff = abs_diff / np.maximum(np.abs(y_base), floor_val)

        l2_diff = np.sqrt(trapezoid(abs_diff**2, x=t_base))
        l2_base = np.sqrt(trapezoid(y_base**2, x=t_base))

        temporal_results[col] = {
            'Max Rel Error (Temporal)': np.max(rel_diff),
            'L2 Rel Error (Temporal)': l2_diff / max(l2_base, floor_val)
        }
    def calc_rel_err(val_approx, val_base):
        return abs(val_approx - val_base) / max(abs(val_base), floor_val)

    final_state_results['Energy'] = {
        'Final Rel Error': calc_rel_err(r_approx.energy, r_base.energy)
    }

    final_state_results['Neutrino Loss'] = {
        'Final Rel Error': calc_rel_err(r_approx.specific_neutrino_energy_loss, r_base.specific_neutrino_energy_loss)
    }

    for sp in species_list:
        try:
            val_base = r_base.composition[sp]
            val_approx = r_approx.composition[sp]
            final_state_results[f"Final {sp}"] = {
                'Final Rel Error': calc_rel_err(val_approx, val_base)
            }
        except (KeyError, TypeError, AttributeError):
            pass

    return pd.DataFrame(temporal_results).T, pd.DataFrame(final_state_results).T

def main(save_show):
    C = init_composition()
    netIn = init_netIn(10**7.1760912591, 10**2.2041199827, 1e18, C)

    stepLogger = StepLogger()

    engine_graph = GraphEngine(C, 4)
    blob = engine_graph.constructStateBlob()
    print(f"Gridfire Using: {len(engine_graph.getNetworkReactions(blob))} Reactions and {len(engine_graph.getNetworkSpecies(blob))} Species")
    print(engine_graph.getNetworkReactions(blob))
    print(engine_graph.getNetworkSpecies(blob))

    solver_ctx_graph = PointSolverContext(blob)
    solver_ctx_graph.stdout_logging = True
    solver_ctx_graph.callback = lambda ctx: stepLogger.log_step(ctx)

    solver_single = PointSolver(engine_graph)

    r_graph = solver_single.evaluate(solver_ctx_graph, netIn, False, False)
    df_graph : pd.DataFrame = stepLogger.df
    df_graph.to_csv("bbq_graph.csv", index=False)
    # stepLogger.reset()
    #
    # QSE_engine = MultiscalePartitioningEngineView(engine_graph)
    # solver_ctx_graph_qse = PointSolverContext(QSE_engine.constructStateBlob(engine_graph.constructStateBlob()))
    # solver_ctx_graph_qse.stdout_logging = False
    # solver_ctx_graph_qse.callback = lambda ctx: stepLogger.log_step(ctx)
    #
    # solver_QSE = PointSolver(QSE_engine)
    # r_qse = solver_QSE.evaluate(solver_ctx_graph_qse, netIn, False, False)
    #
    # df_qse : pd.DataFrame = stepLogger.df
    # df_qse.to_csv("bbq_qse.csv", index=False)
    # stepLogger.reset()


if __name__ == "__main__":
    import argparse
    app = argparse.ArgumentParser(prog="Derivative Smoothness", description="Generate of view plots of derivative smoothness")
    app.add_argument("-s", type=ShowSave, default=ShowSave.SHOW, choices=list(ShowSave), help="Whether to show or save the generated plot")

    args = app.parse_args()
    main(args.s)
