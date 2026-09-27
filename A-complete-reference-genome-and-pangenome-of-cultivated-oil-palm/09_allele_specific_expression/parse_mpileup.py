#!/usr/bin/env python3
"""Parse samtools mpileup (no -f; --no-output-ins/del/ends) from stdin -> chrom pos A C G T."""
import re
import sys

strip = re.compile(r"\^.")
sys.stdout.write("chrom\tpos\tA\tC\tG\tT\n")
for line in sys.stdin:
    t = line.rstrip("\n").split("\t")
    b = strip.sub("", t[4]).replace("$", "").upper()
    sys.stdout.write(f"{t[0]}\t{t[1]}\t{b.count('A')}\t{b.count('C')}\t{b.count('G')}\t{b.count('T')}\n")
