"""Minimal re-implementation of EMMAX (Kang et al. 2010; emmax-intel64 2012-02-10) used for the PC-adjusted
SNP-GWAS, sensitivity scans and conditional tests.

Model: y = X b + g beta + u + e, u ~ N(0, vg K), e ~ N(0, ve I); delta = ve/vg estimated once per phenotype by REML
under the null (EMMA grid of 100 log-delta intervals on [-10, 10] plus root refinement of the derivative), then GLS
with V = K + delta I for every marker; two-sided t test with n - q - 1 d.f. Missing genotypes are mean-imputed over
the analysed individuals (as in EMMAX)."""
import numpy as np
from scipy import stats, optimize

NGRID, LLIM, ULIM = 100, -10.0, 10.0


def reml(y, X, K):
    n, q = X.shape
    S = np.eye(n) - X @ np.linalg.solve(X.T @ X, X.T)
    w, U = np.linalg.eigh(S @ (K + np.eye(n)) @ S)
    order = np.argsort(w)[::-1][: n - q]
    lam = w[order] - 1.0
    eta = U[:, order].T @ y
    e2 = eta ** 2
    m = n - q

    def ll(ld):
        d = np.exp(ld)
        return 0.5 * (m * np.log(m / (2 * np.pi)) - m - m * np.log(np.sum(e2 / (lam + d))) - np.sum(np.log(lam + d)))

    def dll(ld):
        d = np.exp(ld)
        a = np.sum(e2 / (lam + d) ** 2); b = np.sum(e2 / (lam + d))
        return 0.5 * d * (m * a / b - np.sum(1.0 / (lam + d)))

    grid = np.linspace(LLIM, ULIM, NGRID + 1)
    dv = np.array([dll(g) for g in grid])
    cands = [(ll(LLIM), LLIM), (ll(ULIM), ULIM)]
    for i in range(NGRID):
        if dv[i] > 0 and dv[i + 1] < 0:
            r = optimize.brentq(dll, grid[i], grid[i + 1], xtol=1e-12)
            cands.append((ll(r), r))
    best_ll, best = max(cands)
    delta = float(np.exp(best))
    vg = float(np.sum(e2 / (lam + delta)) / m)
    return dict(delta=delta, vg=vg, ve=vg * delta, ll=best_ll, h2=1.0 / (1.0 + delta))


class GLS:
    """Pre-computed transform for one phenotype/covariate/kinship set."""

    def __init__(self, y, X, K):
        self.n, self.q = X.shape
        self.r = reml(y, X, K)
        e, U = np.linalg.eigh(K)
        self.T = (U / np.sqrt(np.clip(e, 0, None) + self.r["delta"])).T      # T' T = (K + delta I)^-1
        ys = self.T @ y
        Xs = self.T @ X
        self.Q, _ = np.linalg.qr(Xs)
        self.ry = ys - self.Q @ (self.Q.T @ ys)
        self.yy = float(self.ry @ self.ry)
        self.df = self.n - self.q - 1

    def scan(self, G):
        """G: (m, n) float genotype dosages (NaN = missing). Returns beta, se, p."""
        G = np.array(G, dtype=np.float64, copy=True)
        miss = np.isnan(G)
        if miss.any():
            mu = np.nanmean(np.where(miss, np.nan, G), axis=1)
            mu = np.where(np.isnan(mu), 0.0, mu)
            G[miss] = np.take(mu, np.nonzero(miss)[0])
        Gs = G @ self.T.T
        rG = Gs - (Gs @ self.Q) @ self.Q.T
        gg = np.einsum("ij,ij->i", rG, rG)
        gy = rG @ self.ry
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = gy / gg
            rss = self.yy - beta * gy
            se = np.sqrt(rss / self.df / gg)
            t = beta / se
        p = 2 * stats.t.sf(np.abs(t), self.df)
        bad = ~np.isfinite(p) | (gg < 1e-10)
        p[bad] = 1.0; beta[bad] = 0.0; se[bad] = np.nan
        return beta, se, p


# --- genotype readers -------------------------------------------------------------------------------------
_LUT = np.full((256, 4), np.nan, dtype=np.float32)
for b in range(256):
    for k in range(4):
        code = (b >> (2 * k)) & 3
        _LUT[b, k] = {0: 2.0, 1: np.nan, 2: 1.0, 3: 0.0}[code]


def read_bed_rows(bed_path, n_ind, start, stop):
    """Dosage (count of bim A1) for SNP rows [start, stop) of a SNP-major PLINK .bed."""
    nb = (n_ind + 3) // 4
    with open(bed_path, "rb") as fh:
        fh.seek(3 + start * nb)
        raw = np.frombuffer(fh.read((stop - start) * nb), dtype=np.uint8).reshape(stop - start, nb)
    return _LUT[raw].reshape(stop - start, nb * 4)[:, :n_ind]


def read_matrix(path):
    return np.loadtxt(path)


def read_pheno(path, ids):
    v = {}
    for line in open(path):
        f = line.split()
        if len(f) >= 3:
            v[f[1]] = np.nan if f[2] in ("NA", "-9", "nan") else float(f[2])
    return np.array([v.get(i, np.nan) for i in ids])


def tped_dosage(fields):
    """tped genotype fields -> dosage of allele '1' counts (EMMAX-like; '0' = missing)."""
    a = np.array(fields, dtype="U8").reshape(-1, 2)
    miss = (a == "0").any(1) | (a == "N").any(1)
    alle = [x for x in np.unique(a) if x not in ("0", "N")]
    if not alle:
        return np.full(len(a), np.nan)
    ref = alle[0]
    d = (a == ref).sum(1).astype(float)
    d[miss] = np.nan
    return d
