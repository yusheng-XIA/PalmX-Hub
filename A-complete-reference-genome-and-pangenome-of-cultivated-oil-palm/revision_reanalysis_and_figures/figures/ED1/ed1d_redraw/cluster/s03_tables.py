#!/usr/bin/env python3
"""Per-read and per-base tables for ED1d from region BAMs (pysam)."""
import sys, pysam, collections
PL = sys.argv[1]
W = "${CLUSTER_WORK}/ed1d_redraw/"
LOCI = [("chr12A_left", "chr12A", 1, 130000, 75515),        # junction: last added base 75,515 | draft-derived seq from 75,516
        ("chr07B_right", "chr07B", 129352881, 129452880, 129402880)]  # draft-derived seq ends 129,402,880 | added from 129,402,881
FLANK = 1000  # read must be aligned >=1 kb on both sides of the junction to count as spanning
for L, C, A, B, J in LOCI:
    bam = pysam.AlignmentFile(W + f"region/{L}.{PL}.bam")
    rd = open(W + f"tables/{L}.{PL}.reads.tsv", "w")
    rd.write("read\tflag\tprimary\tsupplementary\tstrand\tmapq\tref_start_1b\tref_end_1b\taligned_ref_len\tquery_len\t"
             "left_softclip\tright_softclip\tspans_junction_1kb\n")
    for r in bam.fetch(C, A - 1, B):
        if r.is_unmapped or r.is_secondary:
            continue
        cig = r.cigartuples
        lc = cig[0][1] if cig[0][0] in (4, 5) else 0
        rc = cig[-1][1] if cig[-1][0] in (4, 5) else 0
        s1, e1 = r.reference_start + 1, r.reference_end
        span = int(s1 <= J - FLANK + 1 and e1 >= J + FLANK)
        rd.write(f"{r.query_name}\t{r.flag}\t{int(not r.is_supplementary)}\t{int(r.is_supplementary)}\t"
                 f"{'-' if r.is_reverse else '+'}\t{r.mapping_quality}\t{s1}\t{e1}\t{r.reference_length}\t"
                 f"{r.infer_read_length() or 0}\t{lc}\t{rc}\t{span}\n")
    rd.close()
    # per-base depth (all primary+supplementary, and MAPQ>=20), 100-bp bins
    out = open(W + f"tables/{L}.{PL}.depth100.tsv", "w")
    out.write("bin_start_1b\tbin_end_1b\tmean_depth_all\tmean_depth_mapq20\n")
    cov = {}
    for q in (0, 20):
        arr = bam.count_coverage(C, A - 1, B, quality_threshold=0,
                                 read_callback=lambda r, q=q: (not r.is_unmapped and not r.is_secondary
                                                               and r.mapping_quality >= q))
        cov[q] = [sum(x) for x in zip(*arr)]
    n = B - A + 1
    for i in range(0, n, 100):
        j = min(n, i + 100)
        out.write(f"{A + i}\t{A + j - 1}\t{sum(cov[0][i:j]) / (j - i):.2f}\t{sum(cov[20][i:j]) / (j - i):.2f}\n")
    out.close()
print("tables done", PL)
