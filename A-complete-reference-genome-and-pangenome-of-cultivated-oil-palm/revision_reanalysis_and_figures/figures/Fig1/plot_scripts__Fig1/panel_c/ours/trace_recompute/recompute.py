#!/usr/bin/env python3
# Fig1c independent recomputation (pure python). Output: work/trace/Fig1c/out/
import os, re, csv, collections, sys
A = '${ANALYSIS_DIR}'
W = '${CLUSTER_WORK}/trace/Fig1c'
OUT = W + '/out'; os.makedirs(OUT, exist_ok=True)
SD = W + '/sd'
FLFA = A + '/08_hifi_chromosome/05_fianl_all_genome_chrom_only/Africa_hap2.fasta'
EGFA = A + '/20_results/Figure1/04_genome_rename/EG11.fa'
FLGFF = A + '/20_results/10_database/03_genes/Africa_hap2.gff3'
EGGFF = A + '/05_GWAS/00_analysis/00_data/public_data/EG11_genomic.gff'
TEGFF = A + '/11_TE/02_wuzi_EDTA/Africa_hap2_results/Africa_hap2.fa.mod.EDTA.TEanno.gff3'
BISER = A + '/12_SDs/04_BISER/seedless_hap2_SD'
TRF = A + '/12_SDs/03_TRF/Africa_hap2/Africa_hap2.fasta.2.6.6.80.10.50.2000.dat'
PAF = A + '/20_results/Figure1/10_new_figure/02_genomesyn/4genome_run/output/minimap2/1.EG11vsAfrica_hap2.paf'

def norm(c):
    m = re.search(r'chr(\d+)', c)
    return 'chr%02d' % int(m.group(1)) if m else c

def fai(p):
    d = {}
    for l in open(p + '.fai'):
        f = l.split('\t'); d[norm(f[0])] = int(f[1])
    return d

def merge(iv, gap=0):
    iv = sorted(iv); out = []
    for s, e in iv:
        if out and s <= out[-1][1] + gap:
            if e > out[-1][1]: out[-1][1] = e
        else: out.append([s, e])
    return out

def read_sd(name):
    rows = list(csv.reader(open(SD + '/' + name + '.tsv'), delimiter='\t'))
    return rows[0], rows[1:]

rep = open(OUT + '/summary.txt', 'w')
def P(*a):
    print(*a); print(*a, file=rep); rep.flush()

step = sys.argv[1] if len(sys.argv) > 1 else 'all'
fl_len = fai(FLFA); eg_len = fai(EGFA)
P('FL lengths', sum(fl_len.values()), 'EG11', sum(eg_len.values()))

# ---------- 1. gene density ----------
if step in ('all', 'gene'):
    WIN = 200000
    fl = collections.Counter(); nfl = 0
    for l in open(FLGFF):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9 or f[2] != 'gene': continue
        s, e = int(f[3]) - 1, int(f[4]); mid = (s + e) // 2
        fl[(norm(f[0]), mid // WIN)] += 1; nfl += 1
    # EG11: map NC_ to chr by length
    reg = {}; eg = collections.Counter(); egpc = collections.Counter(); neg = 0; negpc = 0; nall = collections.Counter()
    lines = []
    for l in open(EGGFF):
        if l.startswith('#'): continue
        f = l.rstrip('\n').split('\t')
        if len(f) < 9: continue
        if f[2] == 'region' and f[0].startswith('NC_'):
            reg[f[0]] = int(f[4])
        if f[2] == 'gene':
            lines.append((f[0], int(f[3]) - 1, int(f[4]), f[8]))
    bylen = {v: k for k, v in eg_len.items()}
    m = {nc: bylen.get(L) for nc, L in reg.items()}
    P('EG11 NC map', sorted((v, k) for k, v in m.items() if v))
    for c, s, e, a in lines:
        ch = m.get(c)
        nall[(ch is not None, 'protein_coding' in a)] += 1
        if ch is None: continue
        mid = (s + e) // 2
        eg[(ch, mid // WIN)] += 1; neg += 1
        if 'gene_biotype=protein_coding' in a:
            egpc[(ch, mid // WIN)] += 1; negpc += 1
    P('FL genes', nfl, 'EG11 genes on chr (all biotypes)', neg, 'protein_coding', negpc, 'breakdown(onchr,pc)', dict(nall))
    h, rows = read_sd('Fig1c_gene_density')
    diff = collections.Counter(); ex = []
    sdtot = collections.Counter()
    for r in rows:
        smp, ch, st = r[0], r[1], int(float(r[2])); v = int(float(r[4]))
        k = (ch, st // WIN)
        if smp.startswith('EG11'):
            sdtot['EG11'] += v
            for nm, src in (('EG_all', eg), ('EG_pc', egpc)):
                if src[k] != v: diff[nm] += 1
            if egpc[k] != v and len(ex) < 5: ex.append((smp, ch, st, v, egpc[k], eg[k]))
        else:
            sdtot['FL'] += v
            if fl[k] != v:
                diff['FL'] += 1
                if len(ex) < 10: ex.append((smp, ch, st, v, fl[k]))
    P('SD gene totals', dict(sdtot), 'windows differing', dict(diff), 'examples', ex)
    P('FL windows recomputed', len([k for k in fl]), 'EG windows', len(eg))

# ---------- 2. TE density ----------
if step in ('all', 'te'):
    WIN = 1000000
    iv = collections.defaultdict(list); types = collections.Counter(); names = set()
    for l in open(TEGFF):
        if l.startswith('#'): continue
        f = l.split('\t')
        if len(f) < 9: continue
        names.add(f[0]); types[f[2]] += 1
        iv[norm(f[0])].append((int(f[3]) - 1, int(f[4])))
    P('TE gff seqnames sample', sorted(names)[:20], 'types', types.most_common(30))
    cov = {}
    tot = 0
    for ch, L in iv.items():
        mg = merge(L)
        tot += sum(e - s for s, e in mg)
        c = collections.Counter()
        for s, e in mg:
            w = s // WIN
            while s < e:
                we = (w + 1) * WIN; ee = min(e, we)
                c[w] += ee - s; s = ee; w += 1
        cov[ch] = c
    P('TE merged total bp', tot, 'genome frac', tot / sum(fl_len.values()))
    h, rows = read_sd('Fig1c_TE_density')
    mx = 0; nd = 0; ex = []
    for r in rows:
        ch, st, en, v = r[0], int(float(r[1])), int(float(r[2])), float(r[3])
        mine = cov.get(ch, {}).get(st // WIN, 0) / (en - st)
        d = abs(mine - v); mx = max(mx, d)
        if d > 1e-4:
            nd += 1
            if len(ex) < 8: ex.append((ch, st, v, round(mine, 6)))
    P('TE windows', len(rows), 'n diff>1e-4', nd, 'maxdiff', mx, ex)

# ---------- 3. EG11 gaps ----------
if step in ('all', 'gap'):
    gaps = collections.defaultdict(list)
    def flush(ch, buf):
        seq = ''.join(buf)
        for mm in re.finditer(r'[Nn]+', seq): gaps[ch].append((mm.start(), mm.end()))
    ch = None; buf = []
    for l in open(EGFA):
        if l.startswith('>'):
            if ch: flush(ch, buf)
            ch = norm(l[1:].split()[0]); buf = []
        else: buf.append(l.rstrip('\n'))
    if ch: flush(ch, buf)
    # dedupe (runs spanning lines were appended at close only)
    allg = {c: merge(v) for c, v in gaps.items()}
    n = sum(len(v) for v in allg.values()); nb = sum(e - s for v in allg.values() for s, e in v)
    P('EG11 all N-runs', n, 'bp', nb, 'chrs with gaps', len([c for c in allg if allg[c]]), {c: len(v) for c, v in sorted(allg.items())})
    WIN = 200000
    nw = collections.Counter()
    for c, v in allg.items():
        for s, e in v:
            w = s // WIN
            while s < e:
                ee = min(e, (w + 1) * WIN); nw[(c, w)] += ee - s; s = ee; w += 1
    keep = set()
    for (c, w), b in nw.items():
        wl = min(eg_len[c], (w + 1) * WIN) - w * WIN
        if b / wl >= 0.01: keep.add((c, w))
    disp = []
    for c, v in allg.items():
        for s, e in v:
            ws = range(s // WIN, (e - 1) // WIN + 1)
            if any((c, w) in keep for w in ws): disp.append((c, s, e))
    h, rows = read_sd('Fig1c_EG11_gaps')
    sdset = set((r[0], int(float(r[1])), int(float(r[2]))) for r in rows)
    ds = set(disp)
    P('display gaps recomputed', len(ds), 'bp', sum(e - s for c, s, e in ds), 'SD rows', len(sdset), 'SD bp', sum(e - s for c, s, e in sdset),
      'intersect', len(ds & sdset), 'only_mine', list(sorted(ds - sdset))[:5], 'only_SD', list(sorted(sdset - ds))[:5])
    P('EG11 chrs in SD gaps', sorted(set(r[0] for r in rows)))

# ---------- 4. SD (BISER) ----------
if step in ('all', 'sd'):
    iv = collections.defaultdict(list); npair = 0
    for l in open(BISER):
        f = l.split('\t')
        if len(f) < 6: continue
        npair += 1
        iv[norm(f[0])].append((int(f[1]), int(f[2])))
        iv[norm(f[3])].append((int(f[4]), int(f[5])))
    u = {c: merge(v) for c, v in iv.items()}
    n = sum(len(v) for v in u.values()); b = sum(e - s for v in u.values() for s, e in v)
    h, rows = read_sd('Fig1c_segmental_dups')
    sdset = set((r[0], int(float(r[1])), int(float(r[2]))) for r in rows)
    mine = set((c, s, e) for c, v in u.items() for s, e in v)
    P('BISER pairs', npair, 'union intervals', n, 'bp', b, 'SD rows', len(sdset), 'SD bp', sum(e - s for c, s, e in sdset), 'exact match', len(mine & sdset))
    # alternative: pairs with len>=1kb only
    for thr in (1000, 5000):
        iv2 = collections.defaultdict(list)
        for l in open(BISER):
            f = l.split('\t')
            if len(f) < 6: continue
            for c, s, e in ((f[0], int(f[1]), int(f[2])), (f[3], int(f[4]), int(f[5]))):
                if e - s >= thr: iv2[norm(c)].append((s, e))
        u2 = set((c, s, e) for c, v in iv2.items() for s, e in merge(v))
        P(' thr', thr, 'union', len(u2), 'bp', sum(e - s for c, s, e in u2), 'exact match', len(u2 & sdset))

# ---------- 5. TRF ----------
if step in ('all', 'trf'):
    iv = collections.defaultdict(list); ch = None; n = 0
    for l in open(TRF):
        if l.startswith('Sequence:'):
            ch = norm(l.split()[1]); continue
        f = l.split()
        if len(f) >= 14 and f[0].isdigit() and f[1].isdigit():
            iv[ch].append((int(f[0]) - 1, int(f[1]))); n += 1
    h, rows = read_sd('Fig1c_tandem_repeats')
    sdset = set((r[0], int(float(r[1])), int(float(r[2]))) for r in rows)
    P('TRF records', n, 'SD blocks', len(sdset), 'SD bp', sum(e - s for c, s, e in sdset))
    for gap in (1000,):
        for base in ('start-1', 'start'):
            blk = set()
            for c, v in iv.items():
                vv = v if base == 'start-1' else [(s + 1, e) for s, e in v]
                for s, e in merge(vv, gap):
                    if e - s >= 5000: blk.add((c, s, e))
            P(' merge gap', gap, base, 'blocks>=5kb', len(blk), 'bp', sum(e - s for c, s, e in blk), 'exact', len(blk & sdset),
              'only_mine', sorted(blk - sdset)[:4], 'only_SD', sorted(sdset - blk)[:4])

# ---------- 6. synteny ----------
if step in ('all', 'syn'):
    rows_m = []; tot = 0
    for l in open(PAF):
        f = l.split('\t'); tot += 1
        q, qs, qe, strand, t, ts, te, blk, mq = f[0], int(f[2]), int(f[3]), f[4], f[5], int(f[7]), int(f[8]), int(f[10]), int(f[11])
        if norm(q) != norm(t) or mq < 20 or blk < 50000: continue
        rows_m.append((norm(t), ts, te, norm(q), qs, qe, strand, blk, mq))
    h, rows = read_sd('Fig1c_synteny')
    sdset = set((r[0], int(float(r[1])), int(float(r[2])), r[3], int(float(r[4])), int(float(r[5])), r[6], int(float(r[7])), int(float(r[8]))) for r in rows)
    ms = set(rows_m)
    P('PAF lines', tot, 'filtered', len(ms), 'SD rows', len(sdset), 'exact', len(ms & sdset), 'only_mine', sorted(ms - sdset)[:3], 'only_SD', sorted(sdset - ms)[:3])
    plus = sum(r[7] for r in ms if r[6] == '+'); minus = sum(r[7] for r in ms if r[6] == '-')
    P(' + blocks', sum(1 for r in ms if r[6] == '+'), plus, ' - blocks', sum(1 for r in ms if r[6] == '-'), minus)
    # EG11 coverage by filtered same-chr alignments
    cov = collections.defaultdict(list)
    for r in ms: cov[r[0]].append((r[1], r[2]))
    cb = sum(e - s for c, v in cov.items() for s, e in merge(v))
    P(' EG11 bp covered by filtered same-chr blocks', cb, cb / sum(eg_len.values()))
    # large inversion candidates: merge '-' blocks per chr on EG11 within 1 Mb, size>=1Mb
    inv = []
    for c in sorted(set(r[0] for r in ms)):
        v = [(r[1], r[2]) for r in ms if r[0] == c and r[6] == '-']
        for s, e in merge(v, 1000000):
            if e - s >= 1000000: inv.append((c, s, e))
    P(' inverted segments (- blocks merged within 1Mb, >=1Mb on EG11):', len(inv), [(c, round(s/1e6,1), round(e/1e6,1)) for c, s, e in inv])
rep.close()
