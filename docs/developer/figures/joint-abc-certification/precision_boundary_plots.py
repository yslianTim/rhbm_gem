"""Reproduce precision and constrained-boundary figures from certification JSON."""
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt


def read(path):
    return json.loads(path.read_text())


def main():
    root, output = map(Path, sys.argv[1:]); output.mkdir(parents=True, exist_ok=True)
    fig, axis = plt.subplots(figsize=(10, 4), layout='constrained')
    names = ('baseline', 'weak-1e-4', 'near-0.02', 'active-a')
    for offset, variant in enumerate(('legacy', 'guarded', 'guarded-log')):
        x, y = [], []
        for k, name in enumerate(names):
            for path in sorted((root/'datasets'/name/variant/'audits').glob('*.json')):
                p = read(path).get('precision', {})
                if 'maximum_scaled_precision_difference' not in p: continue
                gap = float(p['maximum_scaled_precision_difference'])
                if gap > 0: x.append(k+(offset-1)*.16); y.append(gap)
        axis.scatter(x, y, label=variant, alpha=.7)
    axis.axhline(1e-20, color='black', linestyle='--', label='Agreement threshold')
    axis.set(yscale='log', xticks=range(4), xticklabels=names, ylabel='50/100-digit scaled discrepancy', title='Independent precision agreement')
    axis.legend(ncol=2)
    for suffix in ('png', 'pdf'): fig.savefig(output/f'precision-agreement.{suffix}', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4), layout='constrained')
    for axis, name in zip(axes, ('weak-1e-4', 'near-0.02')):
        value = read(root/'datasets'/name/'legacy/audits/first-stage-double.json')
        for ladder in value['ladders']:
            candidates = [c for c in ladder['candidates'] if c['estimated_relative_error'] is not None]
            axis.loglog([.01*2**-c['index'] for c in candidates],
                        [c['estimated_relative_error'] for c in candidates], '.-', label=f"Direction {ladder['direction']}")
            selected = next((c for c in candidates if c['index'] == ladder['selected']), None)
            if selected:
                axis.scatter([.01*2**-selected['index']], [selected['estimated_relative_error']],
                             marker='*', s=100, color=f"C{ladder['direction']}", edgecolor='black', zorder=5)
        axis.axhline(1e-7, color='black', linestyle='--')
        axis.set(title=name, xlabel='Largest h in the three-step interval', ylabel='Estimated relative uncertainty')
        axis.legend(fontsize=9)
    for suffix in ('png', 'pdf'): fig.savefig(output/f'step-selection.{suffix}', dpi=160)
    plt.close(fig)
    fig, axes = plt.subplots(2, 3, figsize=(13, 7), layout='constrained')
    files = ('truth-boundary-double.json', 'legacy/audits/first-stage-double-boundary.json')
    for row, file in enumerate(files):
        scan = read(root/'datasets/active-a'/file)['precision100']['rows']
        for col, atom in enumerate((1, 5, 9)):
            selected = [r for r in scan if r['atom'] == atom]; axis = axes[row, col]
            for key, label in [('path_objective', 'Compensation path'), ('fixed_face_objective', 'Fixed face'), ('profile_objective', 'Constrained profile')]:
                valid = [r for r in selected if r[key] is not None and float(r[key]) > 0]
                axis.loglog([float(r['step']) for r in valid], [float(r[key]) for r in valid], label=label)
            axis.set(title=f'{"Generating model" if row == 0 else "Actual endpoint"}: atom {atom+1}', xlabel='Positive log-B step', ylabel='RSS / 2')
    axes[0,0].legend()
    for suffix in ('png', 'pdf'): fig.savefig(output/f'boundary-profiles.{suffix}', dpi=160)
    plt.close(fig)


if __name__ == '__main__': main()
