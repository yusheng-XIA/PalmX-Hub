#!/usr/bin/env python3
"""ED1e redraw2: KR-normalise the FL Pore-C contact matrix exactly as HapHiC_plot.py (v1.0.x) does.
Input : juicerbox/contact_matrix.pkl (raw 500-kb matrix written by haphic plot, 2025-09-26; vmax logged 0.008289517522942717)
        juicerbox/04.build/scaffolds_reordered.agp (final assembly completed_genome_v2, chr01A..chr16B)
Output: ED1e_KR_500kb.npz (float32 KR matrix, bin edges, chromosome table, vmax)
bnewt()/normalize_matrix() copied verbatim from HapHiC_plot.py (user/tools/HapHiC, read only).
"""
import pickle, sys, logging
from math import ceil
import numpy as np
logger = logging.getLogger("x"); logging.basicConfig(level=logging.INFO)
def bnewt(A, tol=1e-6, x0=None, delta=0.1, Delta=3, fl=0):

    # A python implemention of the Knight-Ruiz (KR) normalization algorithm described in:
    # https://academic.oup.com/imajna/article/33/3/1029/659457

    n = A.shape[0]
    e = np.ones(n)
    res = []

    if x0 is None:
        x0 = e

    g = 0.9
    etamax = 0.1
    eta = etamax
    stop_tol = tol * 0.5

    x = x0
    rt = tol**2
    v = x * (A @ x)
    rk = 1 - v
    rho_km1 = rk @ rk
    rout = rho_km1
    rold = rout

    MVP = 0
    i = 0

    if fl == 1:
        # print('it in. it res')
        pass

    nn = 0
    max_nn, max_mm = 1000, 10000
    error_message = (
            'Unable to converge. Maybe the matrix is too sparse (too few Hi-C links). '
            'You can try another normalization method.')

    while rout > rt:

        # to avoid endless loop
        nn += 1
        mm = 0
        if nn > max_nn:
            logger.info(error_message)
            raise RuntimeError(error_message)

        i += 1
        k = 0
        y = e
        innertol = max([eta**2 * rout, rt])

        while rho_km1 > innertol:

            # to avoid endless loop
            mm += 1
            if mm > max_mm:
                logger.info(error_message)
                raise RuntimeError(error_message)

            k += 1
            if k == 1:
                Z = rk / v
                p = Z
                rho_km1 = rk @ Z
            else:
                beta = rho_km1 / rho_km2
                p = Z + beta * p

            w = x * (A @ (x * p)) + v * p
            alpha = rho_km1 / (p @ w)
            ap = alpha * p

            ynew = y + ap
            if min(ynew) <= delta:
                if delta == 0:
                    break
                ind = np.where(ap < 0)
                gamma = min((delta - y[ind]) / ap[ind])
                y = y + gamma * ap
                break
            if max(ynew) >= Delta:
                ind = np.where(ynew > Delta)
                gamma = min((Delta - y[ind]) / ap[ind])
                y = y + gamma * ap
                break
            y = ynew
            rk = rk - alpha * w
            rho_km2 = rho_km1
            Z = rk / v
            rho_km1 = rk @ Z

        x = x * y
        v = x * (A @ x)
        rk = 1 - v
        rho_km1 = rk @ rk
        rout = rho_km1
        MVP += k + 1

        rat = rout / rold
        rold = rout
        res_norm = np.sqrt(rout)
        eta_o = eta
        eta = g * rat
        if g * eta_o**2 > 0.1:
            eta = max([eta, g * eta_o**2])
        eta = max([min([eta, etamax]), stop_tol / res_norm])

        if fl == 1:
            # print(f'{i:3d} {k:6d} {res_norm:.3e}')
            res.append(res_norm)

    # print(f'Matrix-vector products = {MVP:6d}')
    return x, res


def normalize_matrix(contact_matrix, group_list, group_size_dict, bin_size, normalization, vmax_coef, manual_vmax):

    if normalization == 'KR':

        # a dict used to save intra-scaffold matrix for each scaffold
        normalized_intra_matrix_dict = dict()

        logger.info('Normalizing contact mattrix using the Knight-Ruiz (KR) balancing algorithm')
        group_start_bin = 0

        # save indices for zero-value elements
        zero_indices = np.argwhere(contact_matrix == 0)
        # to avoid divide-by-zero warning/error
        contact_matrix = contact_matrix + 0.00001

        for group in group_list:
            group_bin_num = ceil(group_size_dict[group]/bin_size)
            group_end_bin = group_start_bin + group_bin_num
            intra_matrix = contact_matrix[group_start_bin:group_end_bin,group_start_bin:group_end_bin]
            group_start_bin += group_bin_num

            # intra-scaffold KR normalization for each scaffold
            x, _ = bnewt(intra_matrix)
            d = np.diag(x)
            normalized_intra_matrix = d @ intra_matrix @ d
            normalized_intra_matrix_dict[group] = normalized_intra_matrix

        # inter-scaffold KR normalization
        x, _ = bnewt(contact_matrix)
        d = np.diag(x)
        normalized_inter_matrix = d @ contact_matrix @ d

        # combine inter- and intra- matrices, and calculate vmax
        non_diagonal_list = list()
        group_start_bin = 0

        for group in group_list:
            group_bin_num = ceil(group_size_dict[group]/bin_size)
            group_end_bin = group_start_bin + group_bin_num
            normalized_inter_matrix[group_start_bin:group_end_bin,group_start_bin:group_end_bin] = normalized_intra_matrix_dict[group]
            for n, row in enumerate(normalized_intra_matrix_dict[group]):
                for m, i in enumerate(row):
                    if n != m:
                        non_diagonal_list.append(i)
            group_start_bin += group_bin_num

        # retrieve zero values
        normalized_inter_matrix[zero_indices[:, 0], zero_indices[:, 1]] = 0

        if manual_vmax < 0:
            vmax = np.median(non_diagonal_list) * vmax_coef
            logger.info('The vmax for the KR-normalized matrix is calculated to be {} ({} * median)'.format(vmax, vmax_coef))
        else:
            vmax = manual_vmax
            logger.info('The vmax for the KR-normalized matrix is manually designated as {})'.format(vmax))

        return normalized_inter_matrix, vmax

    else:

        if normalization == 'log10':
            logger.info('Normalizing contact matrix using log10...')
            normalized_matrix = np.log10(contact_matrix + 1)
        else:
            logger.info('Normalization is disabled')
            normalized_matrix = contact_matrix

        # calculate vmax
        non_diagonal_list = list()
        group_start_bin = 0

        for group in group_list:
            group_bin_num = ceil(group_size_dict[group]/bin_size)
            group_end_bin = group_start_bin + group_bin_num
            intra_matrix = normalized_matrix[group_start_bin:group_end_bin,group_start_bin:group_end_bin]
            for n, row in enumerate(intra_matrix):
                for m, i in enumerate(row):
                    if n != m:
                        non_diagonal_list.append(i)
            group_start_bin += group_bin_num

        if manual_vmax < 0:
            vmax = np.median(non_diagonal_list) * vmax_coef
        else:
            vmax = manual_vmax

        if normalization == 'log10':
            if manual_vmax < 0:
                logger.info('The vmax for the log-normalized matrix is calculated to be {} ({} * median)'.format(vmax, vmax_coef))
            else:
                logger.info('The vmax for the log-normalized matrix is manually designated as {}'.format(vmax))
        else:
            if manual_vmax < 0:
                logger.info('The vmax for the raw matrix is calculated to be {} ({} * median)'.format(vmax, vmax_coef))
            else:
                logger.info('The vmax for the raw matrix is manually designated as {}'.format(vmax))

        return normalized_matrix, vmax



J = "${ANALYSIS_DIR}/01_seedless_results1/18_T2Tassmbly_evaluate/juicerbox"
mat, args = pickle.load(open(J + "/contact_matrix.pkl", "rb"))
print("args", args, "shape", mat.shape, "sum", mat.sum(), flush=True)
bin_size = args.bin_size * 1000   # --bin_size is in kbp (HapHiC main)
gsize = {}
for l in open(J + "/04.build/scaffolds_reordered.agp"):
    f = l.split("	")
    if l.startswith("#") or len(f) < 5: continue
    gsize[f[0]] = max(gsize.get(f[0], 0), int(f[2]))
group_list = [g for g, s in gsize.items() if s >= args.min_len * 1000000]
print("groups", len(group_list), "bins by size//bs+1", sum(gsize[g] // bin_size + 1 for g in group_list), flush=True)
norm, vmax = normalize_matrix(mat, group_list, gsize, bin_size, "KR", getattr(args, "vmax_coef", 4.0), -1)
print("vmax", repr(vmax), flush=True)
np.savez_compressed(sys.argv[1], kr=norm.astype(np.float32), vmax=vmax, bin_size=bin_size,
                    chroms=np.array(group_list), sizes=np.array([gsize[g] for g in group_list]), raw_sum=mat.sum())
print("done", flush=True)
