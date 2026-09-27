#!/bin/bash
# Resistance gene analogues (RGAugury): NBS, RLK, RLP, TM-CC and other classes
# Counts use the 33 biological materials as denominator: for two-haplotype materials, allelic loci are merged
# using orthology and local collinearity (JCVI), while non-collinear paralogues, tandem expansions and
# haplotype-specific candidates are retained.
set -euo pipefail
threads=32

for g in $(cat assemblies_39.txt); do
    RGAugury.pl -p ${g}.pep.fa -n ${g}.cds.fa -g ${g}.fa -gff ${g}.gff3 -c ${threads} -pfx ${g}
done

# Allelic RGA loci between the two haplotypes of one material
python -m jcvi.formats.gff bed --type=mRNA --key=ID ${m}_Hap1.gff3 -o ${m}_Hap1.bed
python -m jcvi.formats.gff bed --type=mRNA --key=ID ${m}_Hap2.gff3 -o ${m}_Hap2.bed
python -m jcvi.compara.catalog ortholog ${m}_Hap1 ${m}_Hap2 --cscore=.99 --no_strip_names
