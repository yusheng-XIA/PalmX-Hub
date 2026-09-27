#!/usr/bin/env python3
"""Aggregate Descals et al. 2021 10-m oil-palm tiles (classes 1 industrial, 2 smallholder, 3 other)
into a global 0.1-degree grid: closed-canopy oil-palm pixel fraction per cell (and per class).
Output: descals_density_0p1deg.npz with arrays palm_ind, palm_sm, total (pixel counts) + grid edges."""
import glob, sys, time
import numpy as np
import rasterio

RES = float(sys.argv[1]) if len(sys.argv) > 1 else 0.1
LON0, LON1, LAT0, LAT1 = -120.0, 180.0, -30.0, 30.0
ncol = int(round((LON1 - LON0) / RES)); nrow = int(round((LAT1 - LAT0) / RES))
ind = np.zeros((nrow, ncol), dtype=np.int64)
sm = np.zeros((nrow, ncol), dtype=np.int64)
tot = np.zeros((nrow, ncol), dtype=np.int64)

def bin_edges(coords, origin):
    """coords: pixel-centre coordinates (1-D, monotonic). returns bin index per pixel and reduceat starts."""
    idx = np.floor((coords - origin) / RES).astype(np.int64)
    change = np.flatnonzero(np.diff(idx)) + 1
    starts = np.concatenate([[0], change])
    return idx[starts], starts

files = sorted(glob.glob("descals_tiles/tiles/*.tif"))
t0 = time.time()
for i, f in enumerate(files, 1):
    with rasterio.open(f) as ds:
        a = ds.read(1)
        h, w = a.shape
        tr = ds.transform
        lons = tr.c + tr.a * (np.arange(w) + 0.5)
        lats = tr.f + tr.e * (np.arange(h) + 0.5)          # tr.e < 0: descending
    ci, cstart = bin_edges(lons, LON0)
    # rows: lats descending -> bin indices descending; reduceat needs contiguous runs, which they are
    ri, rstart = bin_edges(lats, LAT0)
    m_ind = (a == 1); m_sm = (a == 2)
    def agg(mask):
        s = np.add.reduceat(mask.astype(np.int32), cstart, axis=1)
        s = np.add.reduceat(s, rstart, axis=0)
        return s
    s_ind = agg(m_ind); s_sm = agg(m_sm)
    s_tot = np.outer(np.diff(np.append(rstart, h)), np.diff(np.append(cstart, w)))
    rr = np.clip(ri, 0, nrow - 1); cc = np.clip(ci, 0, ncol - 1)
    ok_r = (ri >= 0) & (ri < nrow); ok_c = (ci >= 0) & (ci < ncol)
    rsel = np.ix_(rr[ok_r], cc[ok_c])
    ind[rsel] += s_ind[np.ix_(ok_r, ok_c)]
    sm[rsel] += s_sm[np.ix_(ok_r, ok_c)]
    tot[rsel] += s_tot[np.ix_(ok_r, ok_c)]
    if i % 25 == 0 or i == len(files):
        print(f"{i}/{len(files)} tiles  {time.time()-t0:6.0f} s  palm px so far: {(ind.sum()+sm.sum())/1e6:.1f} M", flush=True)

np.savez_compressed(f"descals_density_{str(RES).replace('.', 'p')}deg.npz", ind=ind, sm=sm, tot=tot,
                    lon_edges=np.linspace(LON0, LON1, ncol + 1), lat_edges=np.linspace(LAT0, LAT1, nrow + 1), res=RES)
frac = (ind + sm) / np.maximum(tot, 1)
print("done. cells with palm:", int((frac > 0).sum()), " max fraction %.3f" % frac.max(),
      " total palm area ~%.2f Mha" % ((ind.sum() + sm.sum()) * 100 / 1e10))
