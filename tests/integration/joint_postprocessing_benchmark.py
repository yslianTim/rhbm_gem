"""Opt-in paired workflow/postprocessing measurements in independent processes."""
import argparse
import json
import platform
import statistics
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--executable', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--compare', type=Path)
    args = parser.parse_args()
    baseline = json.loads(args.compare.read_text()) if args.compare else None
    report = {'platform': platform.platform(), 'repetitions': 3, 'threads': 1,
              'scope': 'Fixture creation excluded from time, included in process peak RSS. Post mode uses a fixed direct-fit endpoint; workflow mode includes First initialization. Output formatting follows RSS capture.',
              'cases': {}}
    with tempfile.TemporaryDirectory(prefix='joint-post-benchmark-') as temp:
        for case in ('full', 'halo', 'multi'):
            for mode in ('workflow', 'post'):
                key = case + '-' + mode
                runs = []
                for repetition in range(3):
                    db = Path(temp) / f'{key}-{repetition}.sqlite'
                    completed = subprocess.run([str(args.executable.resolve()), case, mode, str(db)],
                                               capture_output=True, text=True, check=True, timeout=180)
                    record = json.loads(completed.stdout.splitlines()[-1])
                    record['database_bytes'] = db.stat().st_size
                    runs.append(record)
                entry = {'runs': runs, 'median_phases': {k: statistics.median(r['phases'][k] for r in runs) for k in runs[0]['phases']},
                         'median_peak_rss_bytes': statistics.median(r['process_peak_rss_bytes'] for r in runs)}
                if baseline:
                    before = baseline['cases'][key]
                    for i, r in enumerate(runs):
                        original = before['runs'][i]
                        assert r['endpoints'] == original['endpoints'], key + ': endpoint/convergence changed'
                        assert r['objective'] == original['objective'], key + ': objective changed'
                        assert r['targets'] == original['targets'], key + ': target diagnostics changed'
                    entry['target_diagnostics_exact_equal'] = True
                    entry['phase_ratios'] = {k: v / before['median_phases'][k] for k, v in entry['median_phases'].items()}
                    entry['peak_rss_ratio'] = entry['median_peak_rss_bytes'] / before['median_peak_rss_bytes']
                report['cases'][key] = entry
                print(key, entry['median_phases'], entry['median_peak_rss_bytes'], flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n')


if __name__ == '__main__':
    main()
