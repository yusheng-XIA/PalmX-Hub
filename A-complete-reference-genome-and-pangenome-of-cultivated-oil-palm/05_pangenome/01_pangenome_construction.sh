#!/bin/bash
# Gene-family pangenome (39 genome/haplotype assemblies, 33 biological materials) and graph pangenome
set -euo pipefail
threads=64

# ---- 1. Gene families: OrthoFinder with DIAMOND, FAMSA, FastTree and MCL inflation 1.2 ------------------------
orthofinder -f proteomes_39/ -S diamond -M msa -A famsa -T fasttree -I 1.2 -t ${threads} -a ${threads} -o of39
# 48,920 orthogroups + 40,747 unassigned genes; of the latter, 15,657 retained as GO-supported singleton families
# (normalised amino-acid sequence identical to an eggNOG-annotated protein with >= 1 GO term) -> 64,577 families.
# Occupancy classes over 39 assemblies: core 39, soft-core 38, shell 2-37, cloud 1.
python 02_pangenome_curves.py --families pan_families.tsv --out pan_curves.tsv --perm 1000 --seed 42

# ---- 2. Ks and Ka/Ks within families (Fig. 4h): families reclassified by occurrence across the 33 materials,
#         up to 2,000 families with 2-700 genes per class (seed 42); wgd with MAFFT v7.525 and codeml (PAML v4.10.10)
while read fam; do
    wgd ksd --pairwise -o ksd_${fam} families/${fam}.tsv cds_39.fa
done < sampled_families.txt

# ---- 3. Graph pangenome: minigraph-cactus, 39 assemblies, FL-Hap2 reference, chr01-chr16 only ------------------
# seqfile.txt: "<name>\t<chromosome-only fasta>" for the 39 assemblies (unplaced contigs removed)
cactus-pangenome ./js seqfile.txt --outDir mc39 --outName oilpalm39 --reference FL_Hap2 \
    --giraffe clip filter --gbz clip filter full --gfa clip filter full --vcf \
    --permissiveContigFilter --haplo --chrom-vg clip filter --chrom-og full --viz \
    --maxCores ${threads}
vg stats -lz mc39/oilpalm39.full.gbz > oilpalm39.graph_stats.txt
# Graph sequence beyond the reference path = total node length - FL-Hap2 path length
