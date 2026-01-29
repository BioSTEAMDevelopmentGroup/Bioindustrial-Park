# -*- coding: utf-8 -*-
"""
Baseline Data Generation for HemDx
=================================

Compute and save baseline metrics for HemDx analysis.
"""

from warnings import filterwarnings
filterwarnings('ignore')

import argparse
import json
import os

import numpy as np
import pandas as pd

from biorefineries.prefers.v1.HemDx._models import create_model
from biorefineries.prefers.v1.HemDx._tea_config1 import PreFerSTEA
# HemDx doesn't have get_available_configs exposed yet, likely just config1
# from biorefineries.prefers.v1.HemDx.system import get_available_configs
from biorefineries.prefers.v1.utils import utils, report_generator
from biorefineries.prefers.v1.utils import sankey_utils


# =============================================================================
# CLI
# =============================================================================

def get_available_configs():
    return ['config1']

def parse_arguments():
    parser = argparse.ArgumentParser(description='HemDx baseline data generation')
    parser.add_argument('--config', type=str, default='config1',
                        choices=get_available_configs(),
                        help='Process configuration (default: config1)')
    parser.add_argument('--production', type=float, default=150,
                        help='Baseline production rate in kg/hr (default: 150)')
    parser.add_argument('--timestamp', type=str, default=None,
                        help='Custom timestamp (YYYYMMDD_HHMM). Defaults to current time.')
    parser.add_argument('--results-dir', type=str, default=None,
                        help='Existing results directory to write data into')
    args, _ = parser.parse_known_args()
    return args


# =============================================================================
# Utilities
# =============================================================================

def metric_key_to_str(key):
    if isinstance(key, tuple):
        return f"{key[0]}|{key[1]}"
    return str(key)


def serialize_index(index):
    return list(index) if isinstance(index, tuple) else index


def resolve_output_dirs(script_path, config, timestamp=None, results_dir=None):
    if results_dir:
        base_dir = os.path.abspath(results_dir)
        data_dir = os.path.join(base_dir, 'data')
        figure_dir = os.path.join(base_dir, 'figure')
        os.makedirs(data_dir, exist_ok=True)
        os.makedirs(figure_dir, exist_ok=True)
        return {
            'base': base_dir,
            'data': data_dir,
            'figure': figure_dir,
            'timestamp': os.path.basename(base_dir).split('_')[-1],
            'config': config,
        }
    return utils.get_analysis_dirs(script_path, config=config, timestamp=timestamp)


def _unit_area_from_id(unit_id):
    try:
        digits = ''.join([c for c in str(unit_id) if c.isdigit()])
        if not digits:
            return None
        value = int(digits)
        return int(value // 100) * 100
    except Exception:
        return None


def _stream_summary(system):
    rows = []
    for stream in system.streams:
        if not stream:
            continue
        rows.append({
            'ID': stream.ID,
            'Source': getattr(stream.source, 'ID', None) if stream.source else None,
            'Sink': getattr(stream.sink, 'ID', None) if stream.sink else None,
            'Mass flow [kg/hr]': stream.F_mass,
            'Temperature [K]': stream.T,
            'Pressure [Pa]': stream.P,
            'Enthalpy [kJ/hr]': stream.H,
            'LHV [kJ/kg]': getattr(stream, 'LHV', np.nan),
        })
    return pd.DataFrame(rows)


def _unit_design_summary(system):
    rows = []
    for unit in system.units:
        if not unit:
            continue
        design = getattr(unit, 'design_results', {}) or {}
        base = {
            'ID': unit.ID,
            'Class': unit.__class__.__name__,
            'Area': _unit_area_from_id(unit.ID),
        }
        if not design:
            rows.append(base)
            continue
        for key, value in design.items():
            row = base.copy()
            row['Design parameter'] = key
            row['Value'] = value
            rows.append(row)
    return pd.DataFrame(rows)


# =============================================================================
# Baseline Generation
# =============================================================================

def generate_baseline(config='config1', baseline_production_kg_hr=150, timestamp=None, results_dir=None):
    dirs = resolve_output_dirs(__file__, config, timestamp=timestamp, results_dir=results_dir)

    print("=" * 70)
    print("HEMDX - BASELINE METRICS")
    print("=" * 70)
    print(f"\nConfiguration:")
    print(f"  Process config: {config}")
    print(f"  Baseline production: {baseline_production_kg_hr} kg/hr")
    print(f"  Output data dir: {dirs['data']}")

    model = create_model(
        baseline_production_kg_hr=baseline_production_kg_hr,
        config=config,
        verbose=True,
    )

    print("\nGenerating TEA/LCA breakdown summary...")
    tea = PreFerSTEA(
        system=model.system,
        IRR=0.18,
        duration=(2024, 2044),
        depreciation='IRAS6',
        income_tax=0.17,
        operating_days=333,
        lang_factor=None,
        construction_schedule=(0.15, 0.60, 0.25),
        WC_over_FCI=0.15,
        labor_cost=10*6e4,
        fringe_benefits=0.17+0.07,
        property_tax=0.005,
        property_insurance=0.005,
        supplies=0.02,
        maintenance=0.03,
        administration=0.05,
    )
    generator = report_generator.ProcessReportGenerator(
        system=model.system,
        tea=tea,
        product_stream=model.system.flowsheet.stream.NHemDx_Product,
        config=config,
        timestamp=dirs.get('timestamp'),
        base_dir=dirs['base'],
    )
    generator.generate()

    print("\nEvaluating baseline metrics...")
    baseline_metrics = model.metrics_at_baseline()

    metrics_map = {m.index: m for m in model.metrics}
    baseline_rows = []
    for idx, value in baseline_metrics.items():
        metric = metrics_map.get(idx, None)
        if isinstance(idx, tuple):
            group, name = idx
        else:
            group, name = "", str(idx)
        baseline_rows.append({
            'Group': group,
            'Metric': name,
            'Units': metric.units if metric else "",
            'Value': value,
        })

    baseline_table = pd.DataFrame(baseline_rows)

    baseline_csv = os.path.join(dirs['data'], 'baseline_metrics.csv')
    baseline_xlsx = os.path.join(dirs['data'], 'baseline_metrics.xlsx')
    baseline_json = os.path.join(dirs['data'], 'baseline_metrics.json')

    baseline_table.to_csv(baseline_csv, index=False)
    baseline_table.to_excel(baseline_xlsx, index=False)

    baseline_serialized = {metric_key_to_str(k): v for k, v in baseline_metrics.items()}
    with open(baseline_json, 'w', encoding='utf-8') as f:
        json.dump(baseline_serialized, f, indent=2)

    print("\nGenerating Sankey data...")
    sankey_carbon = sankey_utils.generate_step_sankey_data(
        model.system,
        flow_property='carbon',
        units='kg C/hr',
    )
    sankey_carbon['title'] = 'HemDx Sankey (Carbon Flow)'
    sankey_energy = sankey_utils.generate_step_sankey_data(
        model.system,
        flow_property='energy',
        units='kJ/hr',
    )
    sankey_energy['title'] = 'HemDx Sankey (Energy Flow)'

    sankey_carbon_file = os.path.join(dirs['data'], 'sankey_carbon.json')
    sankey_energy_file = os.path.join(dirs['data'], 'sankey_energy.json')
    with open(sankey_carbon_file, 'w', encoding='utf-8') as f:
        json.dump(sankey_carbon, f, indent=2)
    with open(sankey_energy_file, 'w', encoding='utf-8') as f:
        json.dump(sankey_energy, f, indent=2)

    print(f"  [+] Saved sankey data: {sankey_carbon_file}")
    print(f"  [+] Saved sankey data: {sankey_energy_file}")

    print("\nSaving stream and unit design tables...")
    stream_table = _stream_summary(model.system)
    stream_table.to_csv(os.path.join(dirs['data'], 'stream_summary.csv'))
    
    design_table = _unit_design_summary(model.system)
    design_table.to_csv(os.path.join(dirs['data'], 'unit_design_summary.csv'))

    # -------------------------------------------------------------------------
    # Single Point Sensitivity
    # -------------------------------------------------------------------------
    print("\n" + "=" * 80)
    print("SINGLE POINT SENSITIVITY ANALYSIS")
    print("=" * 80)
    
    # Reset model to baseline first
    for param in model.parameters:
        param.setter(param.baseline)
    
    # Calculate Sensitivity
    print(f"Running single-point sensitivity for {len(model.parameters)} parameters...")
    try:
        # Note: We use a simplified local sensitivity approach where we run
        # Low, Baseline, High for each parameter.
        # This is effectively what bst.Model.spearman/sensitivity does but we want raw values
        # for the Tornado plot.
        
        # Define range (e.g. +/- 10% or using distribution bounds if available)
        # We will use the distribution bounds if they exist, else +/- 10%
        
        sensitivity_results = []
        
        baseline_metrics_values = model([p.baseline for p in model.parameters]) # Run once ensures clean state
        baseline_metrics_dict = {
            m.name_with_units: val for m, val in zip(model.metrics, baseline_metrics_values)
        }
        
        # Check output
        metric_names = [m.name_with_units for m in model.metrics]
        
        # Initialize DataFrame rows
        rows = []

        for param in model.parameters:
             # Determine bounds
             low = param.bounds[0] if param.bounds else param.baseline * 0.9
             high = param.bounds[1] if param.bounds else param.baseline * 1.1
             
             # Run Low
             param.setter(low)
             try:
                 res_low = model.system.simulate()
                 metrics_low = [m() for m in model.metrics]
             except:
                 metrics_low = [np.nan] * len(model.metrics)
             
             # Run High
             param.setter(high)
             try:
                 res_high = model.system.simulate()
                 metrics_high = [m() for m in model.metrics]
             except:
                 metrics_high = [np.nan] * len(model.metrics)
             
             # Reset
             param.setter(param.baseline)
             
             row = {
                 'Parameter': param.name,
                 'Parameter Element': param.element_name,
                 'Baseline Value': param.baseline,
                 'Low Value': low,
                 'High Value': high,
             }
             
             # Add metric results
             for i, m_name in enumerate(metric_names):
                 row[f'{m_name} Low'] = metrics_low[i]
                 row[f'{m_name} High'] = metrics_high[i]
                 row[f'{m_name} Baseline'] = baseline_metrics_dict.get(m_name, np.nan)
             
             rows.append(row)

        sensitivity_df_out = pd.DataFrame(rows)
        
        # Save to Excel
        sens_file = os.path.join(dirs['data'], 'tornado_sensitivity.xlsx')
        with pd.ExcelWriter(sens_file) as writer:
            sensitivity_df_out.to_excel(writer, sheet_name='Summary')
        
        print(f"saved sensitivity data to {sens_file}")
        
    except Exception as e:
        print(f"FAILED to run sensitivity analysis: {e}")
        import traceback
        traceback.print_exc()

    # -------------------------------------------------------------------------
    # Metadata
    # -------------------------------------------------------------------------
    unit_table = _unit_design_summary(model.system)

    streams_csv = os.path.join(dirs['data'], 'baseline_streams.csv')
    streams_xlsx = os.path.join(dirs['data'], 'baseline_streams.xlsx')
    units_csv = os.path.join(dirs['data'], 'baseline_unit_design.csv')
    units_xlsx = os.path.join(dirs['data'], 'baseline_unit_design.xlsx')
    stream_table.to_csv(streams_csv, index=False)
    stream_table.to_excel(streams_xlsx, index=False)
    unit_table.to_csv(units_csv, index=False)
    unit_table.to_excel(units_xlsx, index=False)

    metadata = {
        'config': config,
        'timestamp': dirs.get('timestamp'),
        'baseline_production_kg_hr': baseline_production_kg_hr,
        'parameter_indices': [serialize_index(p.index) for p in model.parameters],
        'metric_indices': [serialize_index(m.index) for m in model.metrics],
        'baseline_metric_key_format': 'group|metric',
    }
    metadata['sankey_carbon_json'] = os.path.basename(sankey_carbon_file)
    metadata['sankey_energy_json'] = os.path.basename(sankey_energy_file)
    metadata['baseline_streams_csv'] = os.path.basename(streams_csv)
    metadata['baseline_streams_xlsx'] = os.path.basename(streams_xlsx)
    metadata['baseline_unit_design_csv'] = os.path.basename(units_csv)
    metadata['baseline_unit_design_xlsx'] = os.path.basename(units_xlsx)
    metadata['breakdown_summary_xlsx'] = 'Breakdown_Summary.xlsx'

    metadata_file = os.path.join(dirs['data'], 'analysis_metadata.json')
    with open(metadata_file, 'w', encoding='utf-8') as f:
        json.dump(metadata, f, indent=2)

    print(f"  [+] Saved baseline metrics: {baseline_csv}")
    print(f"  [+] Saved baseline metrics: {baseline_xlsx}")
    print(f"  [+] Saved baseline metrics: {baseline_json}")
    print(f"  [+] Saved metadata: {metadata_file}")

    return {
        'dirs': dirs,
        'baseline_metrics': baseline_metrics,
        'metadata': metadata,
    }


def main():
    args = parse_arguments()
    generate_baseline(
        config=args.config,
        baseline_production_kg_hr=args.production,
        timestamp=args.timestamp,
        results_dir=args.results_dir,
    )


if __name__ == '__main__':
    main()
