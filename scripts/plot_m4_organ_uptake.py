"""Plot saved oxygen accounting; no simulations or viability thresholds."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'reports/m4-organ-uptake'


def main():
    result = json.loads((OUT/'result.json').read_text(encoding='utf-8'))
    n = len(result['volume_names'])
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), layout='constrained')
    for e in (.5, .75, 1.):
        label = f'quarter-E{e}'
        r = next(r for r in result['cases'] if r['case'] == label and r['solver'] == 'refined')
        data = np.load(OUT/f'{label}-refined.npz'); t = data['time_s']-data['time_s'][0]
        m = len(r['organ_names'])
        axes[0].plot(t, data['oxygen_states'][:, n+2+m:].sum(axis=1), label=f'Extraction cap {e:.0%}')
    axes[0].set(xlabel='Time after arrest (s)', ylabel='Accumulated unmet demand (mL O₂)',
                title='Same hydraulics, different extraction assumptions')
    axes[0].legend(); axes[0].grid(alpha=.2)
    rows = [r for r in result['cases'] if r['solver'] == 'refined']
    labels = [r['case'].replace('quarter-', '').replace('open-', 'All open: ') for r in rows]
    axes[1].barh(labels, [r['total_unmet_mL'] for r in rows], color='#397898')
    axes[1].invert_yaxis(); axes[1].set(xlabel='Unmet demand at 8 s (mL O₂)', title='Partition and extraction sensitivity')
    axes[1].grid(axis='x', alpha=.2)
    fig.suptitle('Exploratory oxygen-deficit proxies — not injury or recovery thresholds', fontsize=13)
    fig.savefig(OUT/'uptake-sensitivity.png', dpi=160)
    fig.savefig(OUT/'uptake-sensitivity.svg')
    plt.close(fig)


if __name__ == '__main__': main()
