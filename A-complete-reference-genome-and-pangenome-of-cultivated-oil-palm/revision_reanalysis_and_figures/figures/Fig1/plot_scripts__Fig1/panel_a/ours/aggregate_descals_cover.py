#!/usr/bin/env python3
"""Aggregate Descals et al. 2019 10-m oil-palm classification (Zenodo 4473715,
oil_palm_map.zip, 634 GeoTIFF tiles L2_2019b_XXXX.tif, EPSG:4326) to
(1) a global 0.05-degree lattice and (2) the 634 ~100-km processing tiles.

Classes in the tiles: 1 = industrial closed-canopy oil palm,
2 = smallholder closed-canopy oil palm, 3 = other land cover.
Oil-palm cover (%) = 100 * (n1 + n2) / (n1 + n2 + n3) within the cell.
Each raster carries a ~0.0046-degree (~51 px) buffer beyond its 100-km grid
cell; pixels are first clipped to the tile's own grid-cell bounds (grid.zip /
grid_withOP.shp, matched by centre) so that overlaps are not double counted.
Pixels are assigned to 0.05-degree cells by pixel-centre coordinates;
cells split between neighbouring tiles are summed across tiles.
"""
import glob, os, sys
import numpy as np
import rasterio
from multiprocessing import Pool

SRC = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'src', 'opmap')
OUT = os.path.dirname(os.path.abspath(__file__))
STEP = 0.05


GRID = os.path.join(OUT, 'grid_withOP_bounds.tsv')


def load_grid():
    import csv
    rows = list(csv.DictReader(open(GRID), delimiter='\t'))
    return [(float(r['xmin']), float(r['ymin']), float(r['xmax']), float(r['ymax']), r['descals_tile_id']) for r in rows]


def proc(path):
    tid = os.path.basename(path).replace('.tif', '')
    with rasterio.open(path) as r:
        a = r.read(1)
        t = r.transform
        b = r.bounds
    cx, cy = (b.left + b.right) / 2, (b.bottom + b.top) / 2
    g = min(load_grid(), key=lambda q: ((q[0] + q[2]) / 2 - cx) ** 2 + ((q[1] + q[3]) / 2 - cy) ** 2)
    if abs((g[0] + g[2]) / 2 - cx) > 0.01 or abs((g[1] + g[3]) / 2 - cy) > 0.01:
        # raster without a row in grid_withOP table: trim the standard buffer
        buf = ((b.right - b.left) - 0.898315) / 2
        g = (b.left + buf, b.bottom + buf, b.right - buf, b.top - buf, 'NA_not_in_grid_table')
    h, w = a.shape
    xc = t.c + (np.arange(w) + 0.5) * t.a
    yc = t.f + (np.arange(h) + 0.5) * t.e
    kx = (xc >= g[0]) & (xc < g[2])
    ky = (yc >= g[1]) & (yc < g[3])
    a = a[np.ix_(ky, kx)]
    xc = xc[kx]; yc = yc[ky]
    b = type('B', (), dict(left=g[0], bottom=g[1], right=g[2], top=g[3]))
    tid = tid + '|grid_id=' + g[4]
    ci = np.floor(xc / STEP + 1e-9).astype(np.int64)
    ri = np.floor(yc / STEP + 1e-9).astype(np.int64)
    cstarts = np.r_[0, np.where(np.diff(ci))[0] + 1]
    rstarts = np.r_[0, np.where(np.diff(ri))[0] + 1]
    res = {}
    for cls in (1, 2, 3):
        m = (a == cls).astype(np.int32)
        s = np.add.reduceat(np.add.reduceat(m, rstarts, axis=0), cstarts, axis=1)
        res[cls] = s
    n1, n2, n3 = int((a == 1).sum()), int((a == 2).sum()), int((a == 3).sum())
    cells = []
    for i, rr in enumerate(ri[rstarts]):
        for j, cc in enumerate(ci[cstarts]):
            cells.append((int(rr), int(cc), int(res[1][i, j]), int(res[2][i, j]), int(res[3][i, j])))
    tile = (tid, b.left, b.bottom, b.right, b.top, n1, n2, n3)
    return tile, cells


if __name__ == '__main__':
    files = sorted(glob.glob(os.path.join(SRC, 'L2_2019b_*.tif')))
    # L2_2019b_3008.tif and L2_2019b_10008.tif cover the same grid cell
    # (descals_tile_id 10008, Guinea; <0.001% oil palm). Keep 10008, whose id
    # matches grid_withOP, so the cell is not counted twice.
    files = [f for f in files if not f.endswith('L2_2019b_3008.tif')]
    print(len(files), 'tiles', file=sys.stderr)
    acc = {}
    owners = {}
    tiles = []
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 6) as pool:
        for k, (tile, cells) in enumerate(pool.imap_unordered(proc, files, chunksize=2)):
            tiles.append(tile)
            for rr, cc, a1, a2, a3 in cells:
                v = acc.setdefault((rr, cc), [0, 0, 0])
                v[0] += a1; v[1] += a2; v[2] += a3
                owners.setdefault((rr, cc), []).append(tile[0])
            if k % 50 == 0:
                print(k, file=sys.stderr)
    tiles.sort()
    with open(os.path.join(OUT, 'Fig1a_cover_tiles_634.tsv'), 'w') as f:
        f.write('tile_id\tlon_min\tlat_min\tlon_max\tlat_max\tlon_centre\tlat_centre\t'
                'n_industrial_px\tn_smallholder_px\tn_other_px\toil_palm_cover_pct\t'
                'industrial_cover_pct\tsmallholder_cover_pct\toil_palm_area_km2_approx\n')
        for tid, x0, y0, x1, y1, n1, n2, n3 in tiles:
            n = n1 + n2 + n3
            f.write(f'{tid}\t{x0:.6f}\t{y0:.6f}\t{x1:.6f}\t{y1:.6f}\t{(x0+x1)/2:.6f}\t{(y0+y1)/2:.6f}\t'
                    f'{n1}\t{n2}\t{n3}\t{100*(n1+n2)/n:.4f}\t{100*n1/n:.4f}\t{100*n2/n:.4f}\t'
                    f'{(n1+n2)*1e-4:.2f}\n')
    with open(os.path.join(OUT, 'Fig1a_cover_grid.tsv'), 'w') as f:
        f.write('cell_id\tlon_min\tlat_min\tlon_max\tlat_max\tlon_centre\tlat_centre\t'
                'n_industrial_px\tn_smallholder_px\tn_other_px\toil_palm_cover_pct\tsource_tiles\n')
        for (rr, cc) in sorted(acc, key=lambda k: (-k[0], k[1])):
            a1, a2, a3 = acc[(rr, cc)]
            if a1 + a2 == 0:
                continue
            n = a1 + a2 + a3
            x0, y0 = cc * STEP, rr * STEP
            f.write(f'c005_{rr:+05d}_{cc:+05d}\t{x0:.2f}\t{y0:.2f}\t{x0+STEP:.2f}\t{y0+STEP:.2f}\t'
                    f'{x0+STEP/2:.3f}\t{y0+STEP/2:.3f}\t{a1}\t{a2}\t{a3}\t{100*(a1+a2)/n:.4f}\t'
                    f'{",".join(sorted(set(owners[(rr, cc)])))}\n')
    print('done', file=sys.stderr)
