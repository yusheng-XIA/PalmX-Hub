#!/bin/bash
# Transposable elements: EDTA v2.1.0 (LTR_FINDER, LTRharvest, LTR_retriever, TIR-Learner) and
# whole-genome annotation with RepeatMasker v4.1.5 in sensitive mode
set -euo pipefail
threads=64

EDTA.pl --genome ${sample}.fa --species others --anno 0 --threads ${threads}
RepeatMasker -pa ${threads} -s -gff -no_is -lib ${sample}.fa.mod.EDTA.TElib.fa ${sample}.fa
# TE superfamilies follow the unified classification (Wicker et al.)

# Insertion time of intact LTR-RTs: T = D / 2mu, mu = 6.5e-9 substitutions per site per year
LTR_retriever -genome ${sample}.fa -inharvest ${sample}.rawLTR.scn -threads ${threads} -u 6.5e-9
# -> ${sample}.fa.pass.list (column "Insertion_Time")

# The TN-Hap2 EDTA library was used to annotate the 29 assemblies used for TE density profiles
# (the TN haplotypes use their own EDTA annotation; TK, NS, Nigerian and E. oleifera use transferred annotations).
