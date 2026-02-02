# -*- coding: utf-8 -*-
"""
Master Runner Script for LegHb and HemDx Data Generation
=========================================================

Runs baseline and Monte Carlo data generation for both LegHb and HemDx models
in sequence with configurable pause intervals between runs, then generates
all figures in parallel.

Run Order:
1. LegHb: gen_data_base.py + gen_data_mc.py (config1) - PARALLEL
2. Pause 5 minutes
3. HemDx: gen_data_base.py + gen_data_mc.py (config1) - PARALLEL
4. Pause 5 minutes  
5. HemDx: gen_data_base.py + gen_data_mc.py (config2) - PARALLEL
6. Pause 5 minutes
7. HemDx: gen_data_base.py + gen_data_mc.py (config3) - PARALLEL
8. Pause 5 minutes
9. Figure Generation (4 parallel): LegHb c1, HemDx c1/c2/c3

Usage:
    python _RUN.py
    python _RUN.py --samples 300000 --batch-size 15000
    python _RUN.py --pause-minutes 10
    python _RUN.py --dry-run  # Show what would run without executing
"""

import argparse
import os
import subprocess
import sys
import time
from datetime import datetime, timedelta

# =============================================================================
# Configuration
# =============================================================================

# Python executable (from user rules)
PYTHON_EXE = r"C:\Users\owenp\.conda\envs\BioSTEAMML\python.exe"

# Script paths (relative to this file)
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

LEGHB_BASE = os.path.join(SCRIPT_DIR, "LegHb", "analyses", "gen_data_base.py")
LEGHB_MC = os.path.join(SCRIPT_DIR, "LegHb", "analyses", "gen_data_mc.py")
LEGHB_FIG = os.path.join(SCRIPT_DIR, "LegHb", "analyses", "gen_figure.py")
HEMDX_BASE = os.path.join(SCRIPT_DIR, "HemDx", "analyses", "gen_data_base.py")
HEMDX_MC = os.path.join(SCRIPT_DIR, "HemDx", "analyses", "gen_data_mc.py")
HEMDX_FIG = os.path.join(SCRIPT_DIR, "HemDx", "analyses", "gen_figure.py")


# =============================================================================
# CLI
# =============================================================================

def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Master runner for LegHb and HemDx data generation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    parser.add_argument('--samples', type=int, default=300000,
                        help='Number of MC samples per scenario (default: 300000)')
    parser.add_argument('--batch-size', type=int, default=15000,
                        help='Batch size for MC evaluation (default: 15000)')
    parser.add_argument('--pause-minutes', type=float, default=5.0,
                        help='Pause duration between stages in minutes (default: 5)')
    parser.add_argument('--dry-run', action='store_true',
                        help='Print commands without executing')
    parser.add_argument('--skip-leghb', action='store_true',
                        help='Skip LegHb runs (start directly with HemDx)')
    parser.add_argument('--skip-hemdx-config1', action='store_true',
                        help='Skip HemDx config1 runs')
    parser.add_argument('--only-leghb', action='store_true',
                        help='Only run LegHb, skip all HemDx')
    parser.add_argument('--only-hemdx', action='store_true',
                        help='Only run HemDx (all configs), skip LegHb')
    parser.add_argument('--timestamp', type=str, default=None,
                        help='Shared timestamp for all runs (YYYYMMDD_HHMM). Auto-generated if not provided.')
    return parser.parse_args()


# =============================================================================
# Utilities
# =============================================================================

def log(msg):
    """Print a timestamped log message."""
    ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"[{ts}] {msg}")


def format_duration(seconds):
    """Format seconds as HH:MM:SS."""
    return str(timedelta(seconds=int(seconds)))


def pause_with_countdown(minutes, dry_run=False):
    """Pause for the specified minutes with countdown display."""
    if dry_run:
        log(f"[DRY-RUN] Would pause for {minutes} minutes")
        return
    
    total_seconds = int(minutes * 60)
    log(f"Pausing for {minutes} minutes ({total_seconds} seconds)...")
    
    for remaining in range(total_seconds, 0, -60):
        mins_left = remaining // 60
        log(f"  ... {mins_left} minute(s) remaining")
        sleep_time = min(60, remaining)
        time.sleep(sleep_time)
    
    log("Pause complete. Resuming...")


def run_scripts_parallel(script_pairs, config, samples, batch_size, timestamp, dry_run=False):
    """
    Run pairs of scripts in parallel using subprocess.Popen.
    
    Args:
        script_pairs: List of (base_script, mc_script) tuples
        config: Configuration name (e.g., 'config1')
        samples: Number of MC samples
        batch_size: Batch size
        timestamp: Shared timestamp
        dry_run: If True, only print commands
        
    Returns:
        List of process exit codes (empty list for dry run)
    """
    processes = []
    commands = []
    
    for base_script, mc_script in script_pairs:
        # Build commands
        base_cmd = [
            PYTHON_EXE, base_script,
            '--config', config,
            '--timestamp', timestamp,
        ]
        
        mc_cmd = [
            PYTHON_EXE, mc_script,
            '--config', config,
            '--samples', str(samples),
            '--batch-size', str(batch_size),
            '--timestamp', timestamp,
        ]
        
        commands.append(('BASE', base_cmd))
        commands.append(('MC', mc_cmd))
    
    # Log all commands
    for label, cmd in commands:
        log(f"[{label}] Command: {' '.join(cmd)}")
    
    if dry_run:
        log("[DRY-RUN] Would execute above commands in parallel")
        return []
    
    # Start all processes
    for label, cmd in commands:
        log(f"Starting {label} process...")
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            encoding='utf-8',
            errors='replace',
        )
        processes.append((label, proc))
    
    # Wait for all processes and stream output
    log(f"Waiting for {len(processes)} processes to complete...")
    
    exit_codes = []
    for label, proc in processes:
        # Stream output as it comes
        try:
            for line in proc.stdout:
                print(f"[{label}] {line.rstrip()}")
        except Exception as e:
            log(f"[{label}] Output stream error: {e}")
        
        proc.wait()
        exit_codes.append(proc.returncode)
        
        status = "SUCCESS" if proc.returncode == 0 else f"FAILED (code {proc.returncode})"
        log(f"[{label}] Process completed: {status}")
    
    return exit_codes


# =============================================================================
# Stage Runners
# =============================================================================

def run_leghb_stage(config, samples, batch_size, timestamp, dry_run=False):
    """Run LegHb baseline and MC data generation."""
    log("=" * 80)
    log(f"STAGE: LegHb {config}")
    log("=" * 80)
    
    return run_scripts_parallel(
        [(LEGHB_BASE, LEGHB_MC)],
        config=config,
        samples=samples,
        batch_size=batch_size,
        timestamp=timestamp,
        dry_run=dry_run,
    )


def run_hemdx_stage(config, samples, batch_size, timestamp, dry_run=False):
    """Run HemDx baseline and MC data generation."""
    log("=" * 80)
    log(f"STAGE: HemDx {config}")
    log("=" * 80)
    
    return run_scripts_parallel(
        [(HEMDX_BASE, HEMDX_MC)],
        config=config,
        samples=samples,
        batch_size=batch_size,
        timestamp=timestamp,
        dry_run=dry_run,
    )


def run_figure_stage(timestamp, dry_run=False, skip_leghb=False, skip_hemdx=False):
    """
    Run figure generation for all configurations in parallel.
    
    Runs 4 scripts in parallel:
    - LegHb config1
    - HemDx config1
    - HemDx config2
    - HemDx config3
    """
    log("=" * 80)
    log("STAGE: Figure Generation (4 parallel)")
    log("=" * 80)
    
    processes = []
    commands = []
    
    # Build command list
    if not skip_leghb:
        commands.append(('LegHb-c1', [PYTHON_EXE, LEGHB_FIG, '--config', 'config1', '--timestamp', timestamp]))
    
    if not skip_hemdx:
        commands.append(('HemDx-c1', [PYTHON_EXE, HEMDX_FIG, '--config', 'config1', '--timestamp', timestamp]))
        commands.append(('HemDx-c2', [PYTHON_EXE, HEMDX_FIG, '--config', 'config2', '--timestamp', timestamp]))
        commands.append(('HemDx-c3', [PYTHON_EXE, HEMDX_FIG, '--config', 'config3', '--timestamp', timestamp]))
    
    # Log all commands
    for label, cmd in commands:
        log(f"[{label}] Command: {' '.join(cmd)}")
    
    if dry_run:
        log("[DRY-RUN] Would execute above commands in parallel")
        return []
    
    # Start all processes
    for label, cmd in commands:
        log(f"Starting {label} figure generation...")
        proc = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
            encoding='utf-8',
            errors='replace',
        )
        processes.append((label, proc))
    
    # Wait for all processes and stream output
    log(f"Waiting for {len(processes)} figure generation processes...")
    
    exit_codes = []
    for label, proc in processes:
        try:
            for line in proc.stdout:
                print(f"[{label}] {line.rstrip()}")
        except Exception as e:
            log(f"[{label}] Output stream error: {e}")
        
        proc.wait()
        exit_codes.append(proc.returncode)
        
        status = "SUCCESS" if proc.returncode == 0 else f"FAILED (code {proc.returncode})"
        log(f"[{label}] Figure generation completed: {status}")
    
    return exit_codes


# =============================================================================
# Main Entry Point
# =============================================================================

def main():
    args = parse_arguments()
    
    # Generate shared timestamp if not provided
    timestamp = args.timestamp or datetime.now().strftime("%Y%m%d_%H%M")
    
    log("=" * 80)
    log("MASTER RUNNER - LegHb & HemDx Data Generation")
    log("=" * 80)
    log(f"Configuration:")
    log(f"  Python: {PYTHON_EXE}")
    log(f"  Samples: {args.samples:,}")
    log(f"  Batch size: {args.batch_size:,}")
    log(f"  Pause between stages: {args.pause_minutes} minutes")
    log(f"  Timestamp: {timestamp}")
    log(f"  Dry run: {args.dry_run}")
    log("")
    
    # Calculate estimated total time
    # Each stage: ~5-7 hours per run, base+MC in parallel
    # We have: LegHb config1 + HemDx config1/2/3 = 4 stages
    expected_stages = []
    if not args.skip_leghb and not args.only_hemdx:
        expected_stages.append("LegHb config1")
    if not args.only_leghb:
        if not args.skip_hemdx_config1:
            expected_stages.append("HemDx config1")
        expected_stages.append("HemDx config2")
        expected_stages.append("HemDx config3")
    expected_stages.append("Figure Generation (4 parallel)")
    
    log(f"Planned stages ({len(expected_stages)}):")
    for i, stage in enumerate(expected_stages, 1):
        log(f"  {i}. {stage}")
    
    # Estimated time: 6 hours per data stage + pauses + ~30 min for figures
    n_data_stages = len(expected_stages) - 1  # Exclude figure stage
    est_hours = n_data_stages * 6 + n_data_stages * (args.pause_minutes / 60) + 0.5
    log(f"\nEstimated total time: ~{est_hours:.1f} hours")
    log("")
    
    start_time = time.time()
    all_exit_codes = []
    
    # ==========================================================================
    # Stage 1: LegHb config1
    # ==========================================================================
    if not args.skip_leghb and not args.only_hemdx:
        codes = run_leghb_stage('config1', args.samples, args.batch_size, timestamp, args.dry_run)
        all_exit_codes.extend(codes)
        
        if not args.only_leghb:
            pause_with_countdown(args.pause_minutes, args.dry_run)
    
    # ==========================================================================
    # Stage 2: HemDx config1
    # ==========================================================================
    if not args.only_leghb and not args.skip_hemdx_config1:
        codes = run_hemdx_stage('config1', args.samples, args.batch_size, timestamp, args.dry_run)
        all_exit_codes.extend(codes)
        pause_with_countdown(args.pause_minutes, args.dry_run)
    
    # ==========================================================================
    # Stage 3: HemDx config2
    # ==========================================================================
    if not args.only_leghb:
        codes = run_hemdx_stage('config2', args.samples, args.batch_size, timestamp, args.dry_run)
        all_exit_codes.extend(codes)
        pause_with_countdown(args.pause_minutes, args.dry_run)
    
    # ==========================================================================
    # Stage 4: HemDx config3
    # ==========================================================================
    if not args.only_leghb:
        codes = run_hemdx_stage('config3', args.samples, args.batch_size, timestamp, args.dry_run)
        all_exit_codes.extend(codes)
    
    # ==========================================================================
    # Pause before Figure Generation (5 minutes)
    # ==========================================================================
    log("")
    log("All data generation stages complete. Pausing before figure generation...")
    pause_with_countdown(args.pause_minutes, args.dry_run)
    
    # ==========================================================================
    # Stage 5: Figure Generation (4 parallel)
    # ==========================================================================
    fig_codes = run_figure_stage(
        timestamp=timestamp,
        dry_run=args.dry_run,
        skip_leghb=args.skip_leghb or args.only_hemdx,
        skip_hemdx=args.only_leghb,
    )
    all_exit_codes.extend(fig_codes)
    
    # ==========================================================================
    # Summary
    # ==========================================================================
    elapsed = time.time() - start_time
    
    log("")
    log("=" * 80)
    log("MASTER RUNNER COMPLETE")
    log("=" * 80)
    log(f"Total elapsed time: {format_duration(elapsed)}")
    
    if all_exit_codes:
        n_success = sum(1 for c in all_exit_codes if c == 0)
        n_failed = len(all_exit_codes) - n_success
        log(f"Processes: {n_success} succeeded, {n_failed} failed")
        
        if n_failed > 0:
            log("WARNING: Some processes failed!")
            return 1
    
    log("All stages completed successfully!")
    return 0


if __name__ == '__main__':
    sys.exit(main())
