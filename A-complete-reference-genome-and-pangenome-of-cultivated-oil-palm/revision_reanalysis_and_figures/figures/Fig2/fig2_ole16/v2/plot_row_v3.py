#!/usr/bin/env python3
"""Fig. 2 bottom row, final (v3): f tree, g OLE16a RNA + FL allele origin, h oleosin:LDAP, i OLE16a across
tissues/materials.  Changes vs v2: tree labels kept clear of the bracket (bracket placed after measuring the
labels), h value labels moved off the y axis, i widened with two-level x labels (material; tissue bracket)."""
import importlib.util, inspect
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('pp', HERE.parent / 'plot_panels.py'); pp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pp)
plt = pp.plt
plt.rcParams['mathtext.sf'] = 'Arial:italic'   # $\mathsf{...}$ = italic gene symbol incl. digits
FL, TN, EOL, DARK, GREY = pp.FL, pp.TN, pp.EOL, pp.DARK, pp.GREY
EGC, BCC, KC = '#8C8C8C', '#7A6FA8', '#B08D57'
W, H = 510.0, 108.0

# tree: label size 5.0 and a bracket position set from the rendered label extent
_src = inspect.getsource(pp.draw_tree).replace('xb = 141 / 510', 'xb = XB / 510').replace('fontsize=5.2', 'fontsize=5.0')
_ns = dict(vars(pp)); _ns['XB'] = 0.0; exec(_src, _ns)

GROUPS = [('EG kernel', 'EG', KC), ('EG mesocarp', 'EG', EGC), ('EO mesocarp', 'EO', EOL),
          ('Backcross mesocarp', 'BC', BCC), ('TN mesocarp', 'TN', TN), ('FL mesocarp', 'FL', FL)]


def draw_i(ax):
    d = pd.read_csv(HERE / 'Fig2i_OLE16a_public_RNA.tsv', sep='\t')
    d = d[d.Plotted.astype(str).str.startswith('yes')]
    rng = np.random.default_rng(3)
    FLOOR = 0.1
    for k, (g, lab, col) in enumerate(GROUPS):
        sub = d[d.Group == g]
        v = sub.OLE16a_RPM.values.astype(float)
        early = sub.Stage.astype(str).str.startswith('1.5').values
        x = k + rng.uniform(-0.2, 0.2, len(v))
        ax.scatter(x[~early], np.maximum(v[~early], FLOOR), s=5, color=col, lw=0, alpha=0.9, zorder=3, clip_on=False)
        ax.scatter(x[early], np.maximum(v[early], FLOOR), s=5, facecolor='white', edgecolor=col, lw=0.5, zorder=3, clip_on=False)
        if len(v):
            ax.plot([k - 0.3, k + 0.3], [max(np.median(v), FLOOR)] * 2, color=DARK, lw=0.7, zorder=4)
    ax.set_yscale('log')
    ax.yaxis.set_minor_locator(pp.matplotlib.ticker.NullLocator())
    ax.set_ylim(0.07, 3e4)
    ax.set_yticks([0.1, 10, 1000], ['0', '10', '10$^{3}$'])
    ax.set_xlim(-0.6, len(GROUPS) - 0.4)
    ax.set_xticks(range(len(GROUPS)), [g[1] for g in GROUPS], fontsize=5.5)
    ax.tick_params(axis='x', length=0, pad=1.5)
    ax.set_ylabel(r'$\mathsf{OLE16a}$ (RPM)', labelpad=1)
    pp.clean(ax)
    import matplotlib.transforms as mt
    tr = mt.blended_transform_factory(ax.transData, ax.transAxes)
    yb, yt = -0.17, -0.21
    for a, b, lab in ((-0.3, 0.3, 'Kernel'), (0.7, 5.3, 'Mesocarp')):
        ax.plot([a, b], [yb, yb], color=GREY, lw=0.6, clip_on=False, transform=tr)
        ax.text((a + b) / 2, yt, lab, ha='center', va='top', fontsize=5.5, color='#555555', clip_on=False, transform=tr)


def draw_ratio(ax):
    pp.draw_ratio(ax)
    for t in list(ax.texts):          # value labels of the first point: move right of the point, off the y axis
        t.remove()
    for ln, lab in zip(ax.get_lines()[:2], ('0.46', '0.005')):
        x, y = ln.get_xdata()[0], ln.get_ydata()[0]
        ax.annotate(lab, (x, y), xytext=(3, 3 if lab == '0.46' else -3), textcoords='offset points', ha='left',
                    va='bottom' if lab == '0.46' else 'top', fontsize=5.5, color=ln.get_color())


def row():
    fig = plt.figure(figsize=(W * pp.PT, H * pp.PT))
    gx, hx, ix = (174, 264), (300, 352), (382, 508)
    fx = lambda a, b: (a / W, b / W)  # noqa: E731
    at = fig.add_axes([0.004, 0.07, 52 / W, 0.86])
    _ns['draw_tree'](at)
    # measure tip labels, place bracket 4 pt right of the longest
    fig.canvas.draw(); rend = fig.canvas.get_renderer()
    brk = [t for t in fig.texts + at.texts if t.get_text() in ('L-oleosins', 'H-oleosins')]
    tips = [t for t in at.texts if t not in brk and t.get_text() != '0.5']
    xmax = max(t.get_window_extent(rend).x1 for t in tips) / fig.dpi * 72
    xb = (xmax + 4) / W
    for ln in at.get_lines():
        xd = ln.get_xdata()
        if len(xd) == 2 and ln.get_transform() is not at.transData and xd[0] == xd[1] == 0.0:
            ln.set_xdata([xb, xb])
    for t in brk:
        t.set_x(xb + 3 / W)
    x0, x1 = fx(*gx)
    ag = fig.add_axes([x0, 0.47, x1 - x0, 0.46]); agb = fig.add_axes([x0, 0.26, x1 - x0, 0.17])
    pp.draw_rna(ag, agb)
    ag.set_ylabel(r'$\mathsf{OLE16a}$ RNA' + '\nlog$_{10}$(count + 1)', labelpad=1)
    lg = agb.get_legend(); hs, ls = lg.legend_handles, [t.get_text() for t in lg.get_texts()]; lg.remove()
    agb.legend(hs, ["FL-Hap1 (" + pp.EO + "-derived)", "FL-Hap2"], loc="upper left", bbox_to_anchor=(-0.01, -0.42),
               ncol=1, frameon=False, handlelength=0.9, handleheight=0.8, labelspacing=0.25, borderaxespad=0,
               fontsize=5.5)
    x0, x1 = fx(*hx)
    ah = fig.add_axes([x0, 0.30, x1 - x0, 0.63]); draw_ratio(ah)
    ah.set_xticks([0, 2, 4, 6], ['185d', '24h', '48h', '72h'], rotation=45, ha='right')
    ah.legend(loc='center left', bbox_to_anchor=(0.02, 0.40), frameon=False, handlelength=1.2, borderaxespad=0.1, ncol=2, columnspacing=0.6)
    x0, x1 = fx(*ix)
    ai = fig.add_axes([x0, 0.30, x1 - x0, 0.63]); draw_i(ai)
    pp.letter(fig, 0, 1, 'f', H); pp.letter(fig, gx[0] - 30, 1, 'g', H)
    pp.letter(fig, hx[0] - 30, 1, 'h', H); pp.letter(fig, ix[0] - 26, 1, 'i', H)
    fig.canvas.draw()
    print('bracket x (pt):', round(xb * W, 1), 'longest label end (pt):', round(xmax, 1))
    fig.savefig(HERE / 'row_ole16_v3.pdf'); fig.savefig(HERE / 'row_ole16_v3.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    row()
