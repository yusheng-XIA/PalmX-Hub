#!/usr/bin/env python3
"""Parse samtools mpileup (no -f, --no-output-ins/del/ends) from stdin -> chrom pos depth A C G T."""
import sys, re
strip = re.compile(r"\^.")
out = sys.stdout
out.write("chrom\tpos\tA\tC\tG\tT\n")
for line in sys.stdin:
    t = line.rstrip("\n").split("\t")
    b = strip.sub("", t[4]).replace("$", "").upper()
    out.write(f"{t[0]}\t{t[1]}\t{b.count('A')}\t{b.count('C')}\t{b.count('G')}\t{b.count('T')}\n")
