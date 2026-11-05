#!/usr/bin/env python3
"""Run independent trials of the Cheimarios 2023 model and export comparisons.

No plotting dependencies are required. Never overwrites an existing experiment.
The paper's Table 2 supplies mean poisoning times, not per-trial data/error bars.
"""
import argparse
import csv
import json
import math
from pathlib import Path
import statistics
import subprocess

PAPER_TIMES = {0.1: 529.9, 0.2: 165.2, 0.3: 104.2, 0.4: 50.8, 0.5: 33.9,
               0.6: 21.4, 0.7: 17.5, 0.8: 13.0, 0.9: 11.7}

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', required=True, type=Path)
    p.add_argument('--output', required=True, type=Path)
    p.add_argument('--representation', choices=['height_field', '3d'], default='height_field')
    p.add_argument('--growth-layers', type=int, default=20)
    p.add_argument('--log-interval', type=float, default=1)
    p.add_argument('--trials', type=int, default=5)
    p.add_argument('--size', type=int, default=100)
    p.add_argument('--y', default='0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9')
    p.add_argument('--variant', choices=['both', 'growth', 'catalytic'], default='both')
    p.add_argument('--growth-duration', type=float, default=2000)
    p.add_argument('--reference-duration', type=float, default=100)
    p.add_argument('--seed', type=int, default=1234)
    p.add_argument('--reaction-selection', choices=['sites', 'pairs'], default='sites')
    p.add_argument('--pair-convention', choices=['zgb_attempts', 'ordered_pairs'], default='zgb_attempts')
    args = p.parse_args()
    try:
        values = [float(v) for v in args.y.split(',')]
    except ValueError:
        p.error('--y must be comma-separated probabilities')
    if (args.trials < 1 or args.growth_layers < 1 or args.size < 3 or args.size > 1000 or args.seed < 1 or
        any(not math.isfinite(v) or not 0 <= v <= 1 for v in values) or not values or len(set(values)) != len(values) or
        not all(math.isfinite(v) and v > 0 for v in [args.growth_duration,args.reference_duration,args.log_interval])):
        p.error('invalid dimensions, trials, seed, probabilities or durations')
    binary = args.binary.resolve()
    if not binary.is_file():
        p.error('binary not found: ' + str(binary))
    output = args.output.resolve()
    if output.exists() and any(output.iterdir()):
        p.error('output directory must be empty; choose a new experiment directory')
    output.mkdir(parents=True, exist_ok=True)
    variants = ['growth', 'catalytic'] if args.variant == 'both' else [args.variant]
    rows = []
    manifest = vars(args).copy()
    manifest.update(binary=str(binary), output=str(output), probabilities=values,
                    paper_doi='10.1088/1361-651X/ace276',
                    growth_observation='coverage at absorption; height is also averaged with exact KMC residence times until absorption; censored at growth_duration',
                    reference_observation='at reference_duration, or earlier absorption',
                    error_bars='sample standard deviation and standard error across independent trials')
    (output / 'manifest.json').write_text(json.dumps(manifest, indent=2) + '\n')
    for variant_index, variant in enumerate(variants):
        for index, y in enumerate(values):
            results = []
            for trial in range(args.trials):
                seed = args.seed + variant_index * 1000000 + index * args.trials + trial
                if seed > 2147483647:
                    raise ValueError('trial seed exceeds supported integer range')
                work = output / f'{variant}-y{y:g}-trial{trial+1}'
                work.mkdir()
                data = {
                    'schema_version': 1,
                    'model': {'type': 'cheimarios_zgb_2023', 'parameters': {
                        'variant': variant, 'y_ab': y, 'reaction_rate': 100000,
                        'pair_convention': args.pair_convention, 'reaction_selection': args.reaction_selection}},
                    'lattice': {'type': 'SimpleCubic', 'parameters': {
                        'dimensions': [args.size, args.size], 'initial_height': 10}},
                    'simulation': {'duration': args.growth_duration if variant == 'growth' else args.reference_duration,
                                   'random_seed': seed},
                    'output': {'log_interval': args.log_interval}}
                summary_name = 'zgb_summary.json'
                if args.representation == '3d':
                    data['model']['type'] = 'cheimarios_zgb_2023_3d'
                    data['lattice'] = {'representation': 'atomistic', 'type': 'SimpleCubic', 'parameters': {
                        'lattice_constant': 1, 'repetitions': [args.size,args.size,10+args.growth_layers],
                        'occupied_unit_cells': 10, 'fixed_unit_cells': 1, 'initial_species': 'substrate',
                        'periodic': [True,True,False]}}
                    data['species'] = {'substrate': {'parameters': {'element': 'X'}}}
                    data['output']['write_xyz'] = False
                    summary_name = 'zgb3d_summary.json'
                (work / 'input.json').write_text(json.dumps(data, indent=2) + '\n')
                result = subprocess.run([str(binary), '--input', 'input.json'], cwd=work,
                                        text=True, capture_output=True, timeout=1800)
                (work / 'run.log').write_text(result.stdout + result.stderr)
                if result.returncode:
                    raise RuntimeError(f'{work}: solver failed: {result.stderr}')
                summary = json.loads((work / summary_name).read_text())
                results.append(summary)
            row = {'variant': variant, 'y_ab': y, 'trials': args.trials,
                   'absorbed_trials': sum(r['poisoned'] for r in results),
                   'censored_trials': sum(not r['poisoned'] for r in results),
                   'capacity_limited_trials': sum(r.get('capacity_limited',False) for r in results)}
            for key in ['time','coverage_ab','coverage_b','coverage_vacant','mean_height','time_average_height','mean_growth_rate','rms']:
                sample = [r[key] for r in results]
                row[key + '_mean'] = statistics.mean(sample)
                row[key + '_sd'] = statistics.stdev(sample) if len(sample) > 1 else 0.0
                row[key + '_sem'] = row[key + '_sd'] / math.sqrt(len(sample))
            paper_time = PAPER_TIMES.get(round(y, 10)) if variant == 'growth' and args.size == 100 else None
            row['paper_table2_time_s'] = '' if paper_time is None else paper_time
            row['relative_table2_time_error'] = (row['time_mean'] / paper_time - 1
                if paper_time is not None and row['censored_trials'] == 0 else '')
            rows.append(row)
            with (output / 'comparison.csv').open('w', newline='') as f:
                writer = csv.DictWriter(f, fieldnames=list(rows[0]))
                writer.writeheader();writer.writerows(rows)
            print(f"{variant:9s} y={y:.3f}: AB={row['coverage_ab_mean']:.4f}, B={row['coverage_b_mean']:.4f}, "
                  f"height={row['mean_height_mean']:.4f}, time={row['time_mean']:.3f}s; "
                  f"absorbed {row['absorbed_trials']}/{args.trials}", flush=True)
    print('Saved ' + str(output / 'comparison.csv'))
    print('Compare configurations and censoring before interpreting deviations as a reproduction.')

if __name__ == '__main__':
    main()
