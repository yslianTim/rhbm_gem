"""Plot accepted and rejected widths separately for each search version."""
import json
from pathlib import Path
import sys
import matplotlib.pyplot as plt

root, output = map(Path, sys.argv[1:])
variants = ('legacy', 'guarded', 'guarded-log')
fig, axes = plt.subplots(2, 3, figsize=(13, 7), sharey=True, layout='constrained')
for row, precision in enumerate(('double', 'float32')):
    for col, variant in enumerate(variants):
        fit = json.loads((root/'datasets/near-0.02'/variant/'fits'/f'narrower-{precision}.json').read_text())
        trials = fit['trials']; axis = axes[row, col]
        axis.semilogy([t['evaluation'] for t in trials], [t['b'][0] for t in trials],
                      ':', color='0.65', label='Candidate sequence')
        accepted = [t for t in trials if t['accepted']]
        rejected = [t for t in trials if not t['accepted']]
        axis.semilogy([t['evaluation'] for t in accepted], [t['b'][0] for t in accepted],
                      'o-', color=f'C{col}', markersize=4, label='Accepted states')
        if rejected:
            axis.scatter([t['evaluation'] for t in rejected], [t['b'][0] for t in rejected],
                         color='crimson', marker='x', zorder=4, label='Rejected candidates')
        axis.set(title=f'{variant} / {precision}', xlabel='Profile evaluation', ylim=(2e-4, 2e4))
        if col == 0: axis.set_ylabel('First atom B (angstrom)')
        if row == 0: axis.legend(fontsize=8, loc='lower right')
for suffix in ('png', 'pdf'):
    fig.savefig(output/f'pathological-trials.{suffix}', dpi=160)
plt.close(fig)
