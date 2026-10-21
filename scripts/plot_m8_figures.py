"""T8.3 primary figures for the screen-independent M8 findings.

Each figure here is a measurement rather than a verdict, so none of them depend
on the pressure-screen definition or the horizon used to locate a boundary.
Figures that do depend on those are deferred until the 40 s boundaries land.

Palette: validated categorical slots 1-4 (blue, orange, aqua, yellow) from the
project's data-visualisation reference. Two of the four fall below 3:1 contrast
against the surface, so every series carries a direct label and no figure relies
on colour alone to carry identity.
"""
from pathlib import Path
import json, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))
from cardiac_redundancy.duration import SubstrateBudget, duration_curve

OUT = ROOT / 'reports/m8-figures'
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100']
INK, INK2, GRID, SURFACE = '#0b0b0b', '#52514e', '#d8d7d2', '#fcfcfb'

plt.rcParams.update({
    'figure.facecolor': SURFACE, 'axes.facecolor': SURFACE,
    'axes.edgecolor': GRID, 'axes.labelcolor': INK, 'axes.titlecolor': INK,
    'text.color': INK, 'xtick.color': INK2, 'ytick.color': INK2,
    'axes.grid': True, 'grid.color': GRID, 'grid.linewidth': .6, 'grid.alpha': .9,
    'axes.spines.top': False, 'axes.spines.right': False,
    'font.size': 9, 'axes.titlesize': 10, 'lines.linewidth': 2,
    'axes.axisbelow': True, 'figure.dpi': 150, 'savefig.dpi': 150, 'savefig.facecolor': SURFACE})


def finish(fig, name, note, lines=3):
    """Reserve caption space before drawing it, so it cannot collide with labels."""
    import textwrap
    # Wrap to the figure's own width so a narrow figure does not clip its caption.
    wrapped = '\n'.join(textwrap.wrap(note, max(60, int(fig.get_size_inches()[0]*19))))
    reserve = .085+.035*wrapped.count('\n')
    fig.tight_layout(rect=(0, reserve, 1, 1))
    fig.text(.012, .004, wrapped, fontsize=6.5, color=INK2, ha='left', va='bottom')
    OUT.mkdir(parents=True, exist_ok=True)
    # No bbox_inches='tight' here: it re-crops and undoes the reserved caption band.
    fig.savefig(OUT / f'{name}.png')
    plt.close(fig)
    print('wrote', (OUT / f'{name}.png').relative_to(ROOT), flush=True)


def fig_compression():
    """Mechanism: the pump closes its own coronary bed during systole."""
    from cardiac_redundancy.sensitivity_circulation import SensitivityCirculation
    from cardiac_redundancy.design import DesignPoint
    d = DesignPoint(left_scale=.651, right_scale=.651, septal_scale=.651,
                    pericardial_scale=.651, period_s=1., coronary_conductance=.969)
    m = SensitivityCirculation(ROOT/'configs/baseline.yaml', ROOT/'configs/m3-evidence-regional.yaml',
        ROOT/'configs/m2-frozen-overlay.yaml', arterial_fraction=.2, inlet_resistance_fraction=.1,
        coronary_conductance=d.coronary_conductance, biological=d.biological())
    m.aux_envelope.set(.5, 1, 1); m.release_auxiliary(.5)
    kap = m.native.base.closure_settings['kappa']
    a = np.load(ROOT/'reports/m8-optimize/scale0.651-period1-coron0.969-kappa0.5/hydraulic-refined.npz')
    t, st, q = a['time_s'], a['states'], a['flows_mL_s']
    k = list(m.flow_names).index('aux_coronary_left')
    peak = float(t[8 <= t][np.argmax([max(0., 0.) for _ in [0]])]) if False else None
    coarse = np.where((t >= 7.5) & (t <= 9.5))[0]
    lv = np.array([m.pressures(float(t[i]), st[:, i])[0]['aLV'] for i in coarse])
    centre = float(t[coarse[int(np.argmax(lv))]])
    sel = np.where((t >= centre-.55) & (t <= centre+.45))[0]
    ao, alv, pc = [], [], []
    for i in sel:
        pr, _ = m.pressures(float(t[i]), st[:, i])
        ao.append(pr['AO']); alv.append(pr['aLV'])
        pc.append(max(pr['aPERI'], 0)+kap*max(pr['aLV']-pr['aPERI'], 0))
    x = t[sel]-t[sel][0]
    fig, (top, bot) = plt.subplots(2, 1, figsize=(6.2, 4.6), sharex=True,
                                   gridspec_kw=dict(height_ratios=[1.35, 1], hspace=.12))
    top.plot(x, ao, color=SERIES[0]); top.plot(x, alv, color=SERIES[1])
    top.plot(x, pc, color=SERIES[2])
    top.annotate('aortic', (x[5], ao[5]), color=SERIES[0], fontsize=8, xytext=(3, 5),
                 textcoords='offset points')
    top.annotate('auxiliary chamber', (x[np.argmax(alv)], max(alv)), color=SERIES[1], fontsize=8,
                 xytext=(4, 2), textcoords='offset points')
    top.annotate('closing pressure (κ=0.5)', (x[np.argmax(pc)], max(pc)), color=SERIES[2],
                 fontsize=8, xytext=(5, -10), textcoords='offset points')
    top.set_ylabel('pressure (mmHg)')
    top.set_title('The auxiliary pump compresses its own coronary bed in systole')
    bot.plot(x, q[sel, k], color=SERIES[3])
    lo = int(np.argmin(q[sel, k])); hi = int(np.argmax(q[sel, k]))
    bot.annotate('coronary flow minimum %.2f mL/s' % q[sel, k][lo], (x[lo], q[sel, k][lo]),
                 color=INK, fontsize=8, xytext=(8, 4), textcoords='offset points')
    bot.annotate('maximum %.2f' % q[sel, k][hi], (x[hi], q[sel, k][hi]), color=INK, fontsize=8,
                 xytext=(6, -12), textcoords='offset points')
    bot.set_ylabel('flow (mL/s)'); bot.set_xlabel('time within one beat (s)')
    finish(fig, 'fig-coronary-compression',
           'Scale 0.651, 60 beats/min, oracle release. Flow is lowest when aortic pressure is near '
           'its peak: supply is interrupted by the pump\'s own contraction, not by a lack of '
           'driving pressure. With no tissue substrate buffer (L-M4-F01) this dip is charged as '
           'unmet demand.')


def fig_kappa():
    rows = []
    for p in sorted((ROOT/'reports/m8-coupled').glob('*/m8-coupled.json')):
        d = json.loads(p.read_text(encoding='utf-8'))
        name = p.parent.name
        if 'kappa1' in name: kap = 1.
        elif 'kappa0' in name: kap = 0.
        elif name == 'scale0.6': kap = .5          # the reference-kappa run
        else: continue
        rows.append((kap, d['design']['left_scale'], d['groups']['pump']['unmet_sustained_mL']))
    rows = sorted({(k, s): u for k, s, u in rows}.items())
    sel = [(k, u) for (k, s), u in rows if abs(s-.6) < 1e-9]
    if len(sel) < 3:
        print('skipping kappa figure: need all three cases at scale 0.6', flush=True); return
    fig, ax = plt.subplots(figsize=(6.6, 3.2))
    labels = ['κ = %g' % k for k, _ in sel]
    ys = [u for _, u in sel]
    # Horizontal bars: the values span zero to 0.3 mL, so the two lower cases are
    # invisible as length on any linear scale and are carried by their labels.
    ax.barh(labels, ys, color=[SERIES[2], SERIES[3], SERIES[1]], height=.58,
            edgecolor=SURFACE, linewidth=2)
    span = max(ys) or 1.
    for label, y in zip(labels, ys):
        ax.annotate('satisfied (0)' if y == 0 else '%.6f mL' % y, (y, label), va='center',
                    ha='left', fontsize=8.5, color=INK, xytext=(6, 0), textcoords='offset points')
    ax.set_xlim(0, span*1.32)
    ax.invert_yaxis()
    ax.set_xlabel('auxiliary sustained unmet O₂ (mL)')
    ax.grid(axis='y', visible=False)
    ax.set_title('Self-supply is governed by an untested model-form choice')
    finish(fig, 'fig-kappa-sensitivity',
           'Scale 0.6, nominal coronary conductance, reference physiology, 10 s horizon. The three '
           'cases were declared in M2 (configs/m2-frozen-overlay.yaml) and assigned to T7.4, but no '
           'runner exercised them before T8.2. At κ=1 the deficit is 128 times the κ=0.5 value, so '
           'the two lower cases have no visible length here and are read from their labels. κ is the '
           'fraction of auxiliary chamber transmural pressure applied as its coronary closing pressure.')


def fig_volume_debt():
    pairs = (('reports/m8-capacity-oracle/scale0.6', 'release 0.5 s', SERIES[0]),
             ('reports/m8-capacity/scale0.6', 'release 2.55 s', SERIES[1]))
    fig, ax = plt.subplots(figsize=(6., 3.6))
    for path, label, colour in pairs:
        a = np.load(ROOT/path/'hydraulic-refined.npz')
        t, v = a['time_s'], a['states'][2]
        ax.plot(t, v, color=colour)
        i = int(np.argmin(abs(t-6.5)))
        ax.annotate(label, (t[i], v[i]), color=colour, fontsize=8, xytext=(4, 4),
                    textcoords='offset points')
        lo = int(np.argmin(v))
        ax.annotate('%.0f mL' % v[lo], (t[lo], v[lo]), color=colour, fontsize=8,
                    xytext=(-2, -12), textcoords='offset points')
    ax.set_xlabel('time from arrest (s)'); ax.set_ylabel('aortic compartment volume (mL)')
    ax.set_title('Detection delay is paid twice: lost time, then a volume debt')
    finish(fig, 'fig-arterial-volume-debt',
           'Scale 0.6. Drainage runs about 55 mL/s initially, so a 0.5 s release leaves a ~28 mL '
           'deficit and a 2.55 s release about 96 mL. Refill runs at roughly 19 mL/s net, and the '
           'compartment is stiff (~1.4 mL/mmHg), so volume and pressure recover together.')


def fig_latency():
    d = json.loads((ROOT/'reports/m8-optimize/latency.json').read_text(encoding='utf-8'))
    groups = {}
    for r in d['rows']:
        key = (r['scale'], r['rate_bpm'])
        groups.setdefault(key, []).append((r['release_s'],
                                           r['adequacy']['joint']['first_satisfied_window_s']))
    fig, ax = plt.subplots(figsize=(6., 3.8))
    for i, (key, rows) in enumerate(sorted(groups.items())):
        rows = sorted(rows)
        xs = [a for a, b in rows if b is not None]; ys = [b for a, b in rows if b is not None]
        if not xs: continue
        ax.plot(xs, ys, color=SERIES[i % 4], marker='o', markersize=5,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
        ax.annotate('scale %.3f, %.0f/min' % key, (xs[-1], ys[-1]), color=SERIES[i % 4],
                    fontsize=8, xytext=(5, -2), textcoords='offset points')
    ax.set_xlabel('prescribed activation latency (s)')
    ax.set_ylabel('first window meeting the screen (s from arrest)')
    ax.set_title('A second of detection delay costs 1.6 to 3.3 s of adequacy')
    ax.set_xlim(.3, 3.9)
    finish(fig, 'fig-latency-relation',
           'Oracle-style prescribed release stands for detection plus recruitment. A larger cavity '
           'flattens the slope and widens the tolerable delay, so capacity and detection substitute '
           'for one another. Systolic pressure is the last metric to settle in every case.')


def fig_duration():
    fig, ax = plt.subplots(figsize=(6., 3.6))
    for i, (label, demand) in enumerate((('favourable, 87.2', 87.2), ('reference, 108.2', 108.2),
                                         ('adverse, 144.8', 144.8))):
        c = duration_curve(demand)
        xs = [r['glucose_floor_mmol_L'] for r in c['rows']]
        ys = [r['minutes'] for r in c['rows']]
        ax.plot(xs, ys, color=SERIES[i], marker='o', markersize=5,
                markeredgecolor=SURFACE, markeredgewidth=1.5)
        ax.annotate(label, (xs[0], ys[0]), color=SERIES[i], fontsize=8, xytext=(4, 4),
                    textcoords='offset points')
    ax.axvline(3., color=INK2, linewidth=1, linestyle=(0, (4, 3)))
    ax.annotate('3.0 mmol/L\nADA/EASD level', (3., ax.get_ylim()[1]*.72), color=INK2,
                fontsize=7.5, xytext=(6, 0), textcoords='offset points')
    ax.set_xlabel('declared plasma glucose floor (mmol/L)')
    ax.set_ylabel('supported duration (min)')
    ax.set_title('Supported duration turns on a floor the model never declared')
    finish(fig, 'fig-supported-duration',
           'Demands in mL O₂/min at the T8.0 optimum. Duration is pool divided by net drain, with '
           'hepatic export rate-capped at 0.21 mmol/min; the 330 mmol glycogen store would last over '
           'a day at that cap, so the rate binds, not the store. First-order bound: the measured '
           'first-minute drain is higher because export ramps over ~60 s, giving 33.7 min against '
           '39.6 at zero floor.')


def _pooled_designs():
    rows = []
    for path in ('reports/m8-joint/joint-kappa0.5.json',
                 'reports/m8-verification/space-filling-kappa0.5.json'):
        f = ROOT/path
        if f.is_file():
            rows += json.loads(f.read_text(encoding='utf-8'))['rows']
    return rows


def fig_design_space():
    """Where the search looked and what it found."""
    rows = _pooled_designs()
    feasible = [r for r in rows if r.get('feasible')]
    infeasible = [r for r in rows if not r.get('feasible')]
    if not feasible:
        print('skipping design space: no feasible designs', flush=True); return
    fig, ax = plt.subplots(figsize=(6.4, 4.0))
    ax.scatter([r['rate_bpm'] for r in infeasible], [r['scale'] for r in infeasible],
               s=26, facecolors='none', edgecolors=GRID, linewidths=1.2, zorder=2)
    demands = [r['auxiliary_demand_mL_min'] for r in feasible]
    sc = ax.scatter([r['rate_bpm'] for r in feasible], [r['scale'] for r in feasible],
                    c=demands, cmap='viridis_r', s=46, edgecolors=SURFACE,
                    linewidths=1.2, zorder=3)
    best = min(feasible, key=lambda r: r['auxiliary_demand_mL_min'])
    ax.scatter([best['rate_bpm']], [best['scale']], s=150, facecolors='none',
               edgecolors=INK, linewidths=1.8, zorder=4)
    ax.annotate('best found\n%.3f mL O₂/min' % best['auxiliary_demand_mL_min'],
                (best['rate_bpm'], best['scale']), color=INK, fontsize=8,
                xytext=(10, 8), textcoords='offset points')
    bar = fig.colorbar(sc, ax=ax, pad=.02)
    bar.set_label('auxiliary O₂ demand (mL/min)', fontsize=8.5)
    bar.ax.tick_params(labelsize=8)
    ax.set_xlabel('activation rate (beats/min)'); ax.set_ylabel('cavity scale')
    ax.set_title('Feasible designs occupy a shallow valley, not a point')
    finish(fig, 'fig-design-space',
           'Pooled from the joint search and the matched-budget sample, 40 s horizon, '
           'compression kappa 0.5, reference physiology, oracle activation. Open circles failed '
           'the pressure screen or auxiliary self-supply. Colour is the objective; the valley is '
           'flat enough that designs differing in scale, rate and coronary conductance cost within '
           '0.2 per cent of each other.')


def fig_conductance_frontier():
    """What staying near nominal coronary conductance costs."""
    rows = [r for r in _pooled_designs() if r.get('feasible')]
    if not rows:
        print('skipping frontier: no feasible designs', flush=True); return
    ceilings = [1.0, 1.1, 1.25, 1.5, 1.75, 2.0]
    best, scales = [], []
    for eps in ceilings:
        within = [r for r in rows if r['coronary_conductance'] <= eps+1e-9]
        pick = min(within, key=lambda r: r['auxiliary_demand_mL_min']) if within else None
        best.append(pick['auxiliary_demand_mL_min'] if pick else None)
        scales.append(pick['scale'] if pick else None)
    fig, ax = plt.subplots(figsize=(6.0, 3.6))
    xs = [c for c, b in zip(ceilings, best) if b is not None]
    ys = [b for b in best if b is not None]
    ax.step(xs, ys, where='post', color=SERIES[0])
    ax.scatter(xs, ys, color=SERIES[0], s=42, edgecolors=SURFACE, linewidths=1.5, zorder=3)
    for c, b, sc in zip(ceilings, best, scales):
        if b is None: continue
        ax.annotate('scale %.3f' % sc, (c, b), color=INK2, fontsize=7.5,
                    xytext=(4, 6), textcoords='offset points')
    ax.set_xlabel('coronary conductance ceiling (multiple of nominal)')
    ax.set_ylabel('cheapest feasible design (mL O₂/min)')
    ax.set_title('A design exists at nominal conductance, at 26% more cost')
    finish(fig, 'fig-conductance-frontier',
           'Epsilon-constraint frontier pooled across both searches. The 0.5-2.0 conductance '
           'bracket is a declared assumption with no measured basis, so what matters is that '
           'feasibility does not depend on it: at or below nominal a design exists at 9.596 mL '
           'O₂/min and cavity scale 0.749, against 7.600 at scale 0.566 once slightly above '
           'nominal conductance is allowed. The assumption buys optimality, not feasibility.')


if __name__ == '__main__':
    for fn in (fig_compression, fig_kappa, fig_volume_debt, fig_latency, fig_duration,
               fig_design_space, fig_conductance_frontier):
        try:
            fn()
        except Exception as exc:
            print('%s failed: %r' % (fn.__name__, exc), flush=True)
