#!/usr/bin/env python3
"""Terminal telomere screen, read-only on FASTA via .fai offsets.

For every chromosome-level sequence and each end it reports
  * tidk-style 10-kb windows (anchored at position 0, as `tidk search -w 10000`):
    TTTAGGG / CCCTAAA counts in the 5 outermost windows (= outermost <=50 kb);
    SD2 rule: max(TTTAGGG, CCCTAAA) in any of those windows >= 3
  * V19 rule (09_genome_eva/05_ragtag_6haps/scripts/fasta_qc.py):
    TTTAGGG + CCCTAAA in outermost 100 kb >= 10
usage: telo_ends.py manifest.tsv out.tsv   (manifest: label<TAB>version<TAB>fasta)
"""
import re, sys

CHR_RE = re.compile(r'^(chr|Chr|CHR|scaffold_)?0*(\d+)[AB]?(_RagTag)?$')
W = 10000


def read_fai(fa):
    out = []
    for l in open(fa + '.fai'):
        p = l.rstrip('\n').split('\t')
        out.append((p[0], int(p[1]), int(p[2]), int(p[3]), int(p[4])))
    return out


def fetch(fh, ent, s, e):
    name, L, off, lb, lw = ent
    s = max(0, s); e = min(L, e)
    a = off + (s // lb) * lw + s % lb
    b = off + ((e - 1) // lb) * lw + (e - 1) % lb + 1
    fh.seek(a)
    return fh.read(b - a).replace(b'\n', b'').replace(b'\r', b'').upper().decode()


def cnt(seq):
    return seq.count('TTTAGGG'), seq.count('CCCTAAA')


def main(man, outp):
    out = open(outp, 'w')
    out.write('label\tversion\tfasta\tseq\tchr_idx\tlength\tend\t'
              'win_counts_F\twin_counts_R\tmax_win\tmax_motif\tmax_win_end\t'
              'sd2_call_ge3\ttotal100kb\tv19_call_ge10\n')
    for line in open(man):
        if not line.strip() or line.startswith('#'):
            continue
        label, ver, fa = line.rstrip('\n').split('\t')
        fai = read_fai(fa)
        chroms = [e for e in fai if CHR_RE.match(e[0])]
        chroms = [e for e in chroms if e[1] > 20_000_000][:16]
        with open(fa, 'rb') as fh:
            for i, ent in enumerate(chroms, 1):
                L = ent[1]
                nwin = (L + W - 1) // W
                for end in ('L', 'R'):
                    idxs = range(0, 5) if end == 'L' else range(nwin - 5, nwin)
                    F, R, best = [], [], (-1, '', 0)
                    for k in idxs:
                        s, e = k * W, min((k + 1) * W, L)
                        f, r = cnt(fetch(fh, ent, s, e))
                        F.append(f); R.append(r)
                        for v, m in ((f, 'TTTAGGG'), (r, 'CCCTAAA')):
                            if v > best[0]:
                                best = (v, m, e)
                    t = fetch(fh, ent, 0, 100000) if end == 'L' else fetch(fh, ent, L - 100000, L)
                    tf, tr = cnt(t)
                    out.write('\t'.join(map(str, [label, ver, fa, ent[0], i, L, end,
                        ','.join(map(str, F)), ','.join(map(str, R)), best[0], best[1], best[2],
                        int(best[0] >= 3), tf + tr, int(tf + tr >= 10)])) + '\n')
                out.flush()
        print(label, ver, len(chroms), 'chromosomes', flush=True)


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2])
