#!/usr/bin/env python3
"""Fig. 2 v2 bottom row: f tree, g OLE16a RNA + FL allele origin, h oleosin:LDAP, i OLE16a across tissues/materials."""
import sys, importlib.util
from pathlib import Path
import numpy as np, pandas as pd
HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('pp', HERE.parent / 'plot_panels.py'); pp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(pp)
import inspect
_src = inspect.getsource(pp.draw_tree).replace('xb = 141 / 510', 'xb = XB / 510').replace('fontsize=5.2', 'fontsize=5.0')
_ns = dict(vars(pp)); _ns['XB'] = 104; exec(_src, _ns); draw_tree = _ns['draw_tree']
plt = pp.plt
FL, TN, EOL, DARK, GREY = pp.FL, pp.TN, pp.EOL, pp.DARK, pp.GREY
EGC, BCC = '#8C8C8C', '#7A6FA8'

GROUPS = [('EG kernel', 'EG kernel', '#B08D57'), ('EG mesocarp', 'EG meso.', EGC),
          ('EO mesocarp', 'EO meso.', EOL), ('Backcross mesocarp', 'BC meso.', BCC),
          ('TN mesocarp', 'TN meso.', TN), ('FL mesocarp', 'FL meso.', FL)]


def draw_i(ax):
    d = pd.read_csv(HERE / 'Fig2i_OLE16a_public_RNA.tsv', sep='\t')
    d = d[d.Plotted.astype(str).str.startswith('yes')]
    rng = np.random.default_rng(3)
    FLOOR = 0.1
    for k, (g, lab, col) in enumerate(GROUPS):
        sub = d[d.Group == g]
        v = sub.OLE16a_RPM.values.astype(float)
        early = sub.Stage.astype(str).str.startswith('1.5').values   # kernel/mesocarp 1.5 months: open symbols
        x = k + rng.uniform(-0.22, 0.22, len(v))
        ax.scatter(x[~early], np.maximum(v[~early], FLOOR), s=5, color=col, lw=0, alpha=0.9, zorder=3, clip_on=False)
        ax.scatter(x[early], np.maximum(v[early], FLOOR), s=5, facecolor='white', edgecolor=col, lw=0.5, zorder=3, clip_on=False)
        if len(v):
            ax.plot([k - 0.3, k + 0.3], [max(np.median(v), FLOOR)] * 2, color=DARK, lw=0.7, zorder=4)
    ax.set_yscale('log')
    ax.yaxis.set_minor_locator(pp.matplotlib.ticker.NullLocator())
    ax.set_ylim(0.07, 3e4)
    ax.set_yticks([0.1, 10, 1000], ['0', '10', '10$^{3}$'])
    ax.set_xlim(-0.6, len(GROUPS) - 0.4)
    ax.set_xticks(range(len(GROUPS)), [g[1] for g in GROUPS], fontsize=5.5, rotation=45, ha='right', rotation_mode='anchor')
    ax.tick_params(axis='x', length=0, pad=1.5)
    ax.set_ylabel('OLE16a (RPM)', labelpad=1)
    pp.clean(ax)


def row():
    W, H = 510.0, 108.0
    fig = plt.figure(figsize=(W * pp.PT, H * pp.PT))
    tx, gx, hx, ix = (0, 150), (180, 284), (326, 388), (420, 508)
    fx = lambda a, b: (a / W, b / W)  # noqa: E731
    at = fig.add_axes([tx[0] / W + 0.004, 0.07, (tx[1] - tx[0]) / W * 0.36, 0.86])
    draw_tree(at)
    x0, x1 = fx(*gx)
    ag = fig.add_axes([x0, 0.47, x1 - x0, 0.46]); agb = fig.add_axes([x0, 0.26, x1 - x0, 0.17])
    pp.draw_rna(ag, agb)
    x0, x1 = fx(*hx)
    ah = fig.add_axes([x0, 0.30, x1 - x0, 0.63]); pp.draw_ratio(ah)
    ah.set_xticks([0, 2, 4, 6], ['185d', '24h', '48h', '72h'], rotation=45, ha='right')
    ah.legend(loc='center left', bbox_to_anchor=(0.02, 0.40), frameon=False, handlelength=1.2, borderaxespad=0.1, ncol=2, columnspacing=0.6)
    x0, x1 = fx(*ix)
    ai = fig.add_axes([x0, 0.30, x1 - x0, 0.63]); draw_i(ai)
    pp.letter(fig, 0, 1, 'f', H); pp.letter(fig, gx[0] - 28, 1, 'g', H)
    pp.letter(fig, hx[0] - 28, 1, 'h', H); pp.letter(fig, ix[0] - 26, 1, 'i', H)
    fig.savefig(HERE / 'row_ole16_v2.pdf'); fig.savefig(HERE / 'row_ole16_v2.png', dpi=300)
    plt.close(fig)


if __name__ == '__main__':
    row()
