#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat August 23 17:00:00 2025

Microalgae biorefinery to produce medium chain fatty acids 
by anaerobic fermentation without external electron donor addition- TRY analysis

References
----------
[1] BioSTEAM Documentation: 
[2] Cortes-Peña et al., BioSTEAM: A Fast and Flexible Platform for the Design, 
    Simulation, and Techno-Economic Analysis of Biorefineries under Uncertainty. 
    ACS Sustainable Chem. Eng. 2020, 8 (8), 3302–3310.
[3] succinic biorefineries project:
    https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/master/biorefineries/succinic

@author: Xingdong Shi
@version: 0.0.1
"""

import os
from datetime import datetime
import numpy as np
import pandas as pd
import biosteam as bst
import matplotlib.pyplot as plt
from ._chemicals import chems
from .system import create_microalgae_MCCA_production_sys
from .tea import microalgae_tea
from .lca import create_microalgae_lca

def run_TRY_analysis(steps: int = 30, plot: bool = True):


    bst.settings.set_thermo(chems)
    sys = create_microalgae_MCCA_production_sys()
    sys.simulate()

    u = sys.flowsheet.unit
    s = sys.flowsheet.stream
    R301 = next((unit for unit in sys.units if unit.ID == 'R301'), None)
    if R301 is None:
        raise RuntimeError('no R301。')

    baseline_titer = float(getattr(R301, 'titer', 2.003))
    baseline_yield_factor = float(getattr(R301, 'caproic_acid_yield_factor', 1.0))
    baseline_abs_C6_yield = getattr(R301, 'caproic_acid_yield_absolute', None)
    if baseline_abs_C6_yield is None:
        baseline_abs_C6_yield = 0.27 * baseline_yield_factor  # fallback to base 0.27 * factor

    titers = np.linspace(0.1, 15.0, steps)
    # Scan absolute C6 yields directly (fraction of algae mass to CaproicAcid)
    yields = np.linspace(0.5 * baseline_abs_C6_yield, 1.5 * baseline_abs_C6_yield, steps)

    tea_obj = microalgae_tea(sys)
    lca = create_microalgae_lca(sys, s.caproic_acid_product, ['CaproicAcid'], u.BT601)

    MPSP = np.full((steps, steps), np.nan)
    GWP = np.full((steps, steps), np.nan)
    FEC = np.full((steps, steps), np.nan)
    EFF_TITER = np.full((steps, steps), np.nan)

    orig_titer = float(getattr(R301, 'titer', baseline_titer))
    orig_abs_yield = getattr(R301, 'caproic_acid_yield_absolute', baseline_abs_C6_yield)

    set_titer = getattr(R301, 'set_C6_titer', lambda x: setattr(R301, 'titer', float(x)))
    set_yf = getattr(R301, 'set_C6_yield', lambda x: setattr(R301, 'caproic_acid_yield_absolute', float(x)))

    # Helper to evaluate across TRY like succinic
    def _evaluate_across_TRY():
        for i, ti in enumerate(titers):
            set_titer(ti)
            for j, yf in enumerate(yields):
                set_yf(yf)
                try:
                    sys.simulate()
                    MPSP[i, j] = tea_obj.solve_price(s.caproic_acid_product)
                    GWP[i, j] = lca.GWP
                    FEC[i, j] = lca.FEC
                    EFF_TITER[i, j] = float(getattr(R301, 'effluent_titer', np.nan))
                except Exception:
                    MPSP[i, j] = np.nan
                    GWP[i, j] = np.nan
                    FEC[i, j] = np.nan
                    EFF_TITER[i, j] = np.nan
    _evaluate_across_TRY()

    set_titer(orig_titer)
    set_yf(orig_abs_yield)
    sys.simulate()

    base = os.path.dirname(__file__)
    results_dir = os.path.join(base, 'analyses', 'results')
    os.makedirs(results_dir, exist_ok=True)
    now = datetime.now()
    tag = f"microalgae_TRY_{steps}x{steps}_{now.year}.{now.month}.{now.day}-{now.hour}.{now.minute:02d}"

    np.save(os.path.join(results_dir, tag + '_MPSP.npy'), MPSP)
    np.save(os.path.join(results_dir, tag + '_GWP.npy'), GWP)
    np.save(os.path.join(results_dir, tag + '_FEC.npy'), FEC)

    pd.DataFrame(MPSP, index=np.round(titers, 4), columns=np.round(yields, 4)).to_csv(
        os.path.join(results_dir, 'MPSP-' + tag + '.csv')
    )
    pd.DataFrame(GWP, index=np.round(titers, 4), columns=np.round(yields, 4)).to_csv(
        os.path.join(results_dir, 'GWP-' + tag + '.csv')
    )
    pd.DataFrame(FEC, index=np.round(titers, 4), columns=np.round(yields, 4)).to_csv(
        os.path.join(results_dir, 'FEC-' + tag + '.csv')
    )
    pd.DataFrame(EFF_TITER, index=np.round(titers, 4), columns=np.round(yields, 4)).to_csv(
        os.path.join(results_dir, 'EFF_TITER-' + tag + '.csv')
    )

    X, Y = np.meshgrid(np.round(yields, 4), np.round(titers, 4))  # X: Yield factor, Y: Titer

    def _auto_levels(Z, n=12):
        z = Z[~np.isnan(Z)]
        if z.size == 0:
            return np.linspace(0, 1, n)
        vmin, vmax = float(np.nanmin(z)), float(np.nanmax(z))
        if vmin == vmax:
            vmax = vmin + 1e-6
        return np.linspace(vmin, vmax, n)

    def _plot_and_save(Z, title, cbar_label, fname):
        levels = _auto_levels(Z)
        fig, ax = plt.subplots(figsize=(6, 4.5), dpi=150)
        cf = ax.contourf(X, Y, Z, levels=levels, cmap='YlGn')
        c = ax.contour(X, Y, Z, levels=levels[::2], colors='k', linewidths=0.6)
        ax.clabel(c, inline=True, fontsize=8, fmt='%.2f')
        cb = fig.colorbar(cf, ax=ax)
        cb.set_label(cbar_label)
        ax.set_xlabel('Yield factor [-]')
        ax.set_ylabel('Titer [g L$^{-1}$]')
        ax.set_title(title)
        fig.tight_layout()
        fig.savefig(os.path.join(results_dir, fname), dpi=300)
        plt.close(fig)

    if plot:
        _plot_and_save(MPSP, 'MPSP contour', 'MPSP [$/kg]', f'MPSP_contour-{tag}.png')
        _plot_and_save(GWP, 'GWP$_{100}$ contour', 'kg CO$_2$-eq (kg$^{-1}$)', f'GWP_contour-{tag}.png')
        _plot_and_save(FEC, 'FEC contour', 'MJ (kg$^{-1}$)', f'FEC_contour-{tag}.png')
        _plot_and_save(EFF_TITER, 'Effluent titer contour', 'g L$^{-1}$', f'EFF_TITER_contour-{tag}.png')

    print(f"TRY saved to: {results_dir} (tag: {tag})")
    return {
        'titers': titers,
        'yields': yields,
        'MPSP': MPSP,
        'GWP': GWP,
        'FEC': FEC,
        'EFF_TITER': EFF_TITER,
        'results_dir': results_dir,
        'tag': tag,
    }

def analyze_TRY_results(results_df, save_path=None):
    """Analyze TRY results and identify optimal regions"""
    if results_df is None:
        print("No results to analyze")
        return
    
    print("\nAnalyzing TRY results...")
    
    # Find optimal conditions (minimum MFSP)
    if 'MFSP' in results_df.columns:
        min_mfsp_idx = results_df['MFSP'].idxmin()
        optimal_conditions = results_df.loc[min_mfsp_idx]
        print(f"\nOptimal conditions (minimum MFSP):")
        print(f"  Titer: {optimal_conditions['Titer']:.2f}")
        print(f"  Yield: {optimal_conditions['Yield']:.3f}")
        print(f"  MFSP: {optimal_conditions['MFSP']:.4f} $/kg")
        if 'GWP' in results_df.columns:
            print(f"  GWP: {optimal_conditions['GWP']:.4f} kg CO2-eq/kg")
    
    # Save results if path provided
    if save_path:
        results_df.to_csv(save_path, index=False)
        print(f"\nResults saved to: {save_path}")

if __name__ == '__main__':
    try:
        run_TRY_analysis(plot=True)
    except Exception as e:
        print(f"Error in TRY analysis: {e}")
        import traceback
        traceback.print_exc()