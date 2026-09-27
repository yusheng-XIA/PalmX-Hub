import itertools, numpy as np, importlib.util, sys, io, contextlib
# brute-force check of kbreak_dp and subset_dp_values on small random matrices
src = open('enh_C.py').read()
ns = {}
exec(src[src.index('def kbreak_dp'):src.index('rows = []')], {'np': np}, ns)
exec(src[src.index('def subset_dp_values'):src.index('def subset_dp_path')], {'np': np}, ns)
rng = np.random.default_rng(1); ok = 0; n_t = 0
for t in range(300):
    n, k = rng.integers(2, 7), rng.integers(2, 5)
    C = rng.integers(0, 6, (n, k)).astype(np.int64)
    for K in (0, 1, 2):
        best = min(sum(C[i, p[i]] for i in range(n)) for p in itertools.product(range(k), repeat=n)
                   if sum(p[i] != p[i+1] for i in range(n-1)) <= K)
        path = ns['kbreak_dp'](C, K)
        n_t += 1; ok += (C[np.arange(n), path].sum() == best and (np.diff(path) != 0).sum() <= K)
    if k >= 3:
        subs = np.array(list(itertools.combinations(range(k), 2)))
        v = ns['subset_dp_values'](C, subs, 2)
        for j, s in enumerate(subs):
            best = min(sum(C[i, p[i]] for i in range(n)) for p in itertools.product(list(s), repeat=n)
                       if sum(p[i] != p[i+1] for i in range(n-1)) <= 2)
            n_t += 1; ok += (v[j] == best)
print(f"brute-force tests passed {ok}/{n_t}")
