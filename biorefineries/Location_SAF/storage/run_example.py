"""Run the coupled storage model without changing the study's default constants."""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples', type=int, default=2)
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--output', type=Path, default=Path('outputs/example.csv'))
    parser.add_argument('--plot', action='store_true')
    args = parser.parse_args()
    if args.samples < 1:
        parser.error('--samples must be positive')
    print('Loading GIS layers and coefficient interpolators; this can take several minutes.', flush=True)
    from coupled_model_mc_10 import run_coupled_mc, plot_coupled_mc
    df = run_coupled_mc(N=args.samples, seed=args.seed, verbose=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(args.output, index=False)
    print(df.groupby('combo')[['msp', 'msp_adj', 'feedstock_delivered_cost_M']].median())
    print(f'Saved {len(df)} rows to {args.output.resolve()}')
    if args.plot:
        plot_coupled_mc(df, save_path=str(args.output.with_suffix('.png')))


if __name__ == '__main__':
    main()
