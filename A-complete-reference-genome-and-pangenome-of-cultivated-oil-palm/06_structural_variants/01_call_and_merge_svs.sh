#!/bin/bash
# Linear SV atlas: each of the 38 non-reference assembly paths vs FL-Hap2 (chr01B-chr16B)
# minimap2 v2.30-r1287, SyRI v1.7.1, SVIM-asm v1.0.3, cuteSV v2.1.3, SURVIVOR v1.0.7, SAMtools v1.23.1
set -euo pipefail
threads=32
REF=FL_Hap2.chr16.fa
REPEAT=FL_Hap2.fa.mod.EDTA.TEanno.bed
FAI=FL_Hap2.fa.fai

# paths.tsv: sample  assembly_fasta  hifi_reads   (both haplotypes of one individual share its HiFi call set)
while IFS=$'\t' read -r sample query reads; do
    out=callers/${sample}; mkdir -p ${out}/syri ${out}/svim ${out}/cutesv

    # --- assembly-based: SyRI and SVIM-asm on the same asm5 alignment ---
    minimap2 -ax asm5 -t ${threads} --eqx ${REF} ${query} -o ${out}/syri/alignment.sam
    (cd ${out}/syri && syri -c alignment.sam -r ../../../${REF} -q ${query} -k -F S --nc ${threads} \
        --samplename ${sample} --lf syri.log)
    samtools sort -@ ${threads} -o ${out}/svim/${sample}.asm5.bam ${out}/syri/alignment.sam
    samtools index ${out}/svim/${sample}.asm5.bam
    svim-asm haploid ${out}/svim/svim_output ${out}/svim/${sample}.asm5.bam ${REF} \
        --min_sv_size 30 --max_sv_size 100000 --sample ${sample}.SVIM_asm

    # --- read-based: cuteSV on HiFi reads ---
    minimap2 -ax map-hifi -t ${threads} --MD -R "@RG\tID:${sample}\tSM:${sample}\tPL:PACBIO" ${REF} ${reads} \
        | samtools sort -@ 4 -o ${out}/cutesv/${sample}.hifi.bam
    samtools index ${out}/cutesv/${sample}.hifi.bam
    cuteSV ${out}/cutesv/${sample}.hifi.bam ${REF} ${out}/cutesv/${sample}.cuteSV.vcf ${out}/cutesv/tmp \
        --sample ${sample} --min_size 30 --max_size 100000 \
        --max_cluster_bias_INS 1000 --diff_ratio_merging_INS 0.9 \
        --max_cluster_bias_DEL 1000 --diff_ratio_merging_DEL 0.5 \
        --min_support 3 -t ${threads} --genotype

    # --- normalise to chr01B-chr16B; drop length-bearing events < 50 bp (breakpoint-only TRA kept) ---
    mkdir -p norm
    for caller in syri svimasm cutesv; do
        case ${caller} in
            syri)    src=${out}/syri/syri.vcf ;;
            svimasm) src=${out}/svim/svim_output/variants.vcf ;;
            cutesv)  src=${out}/cutesv/${sample}.cuteSV.vcf ;;
        esac
        python normalize_oilpalm_sv_vcf.py --caller ${caller} --in ${src} --out norm/${sample}.${caller}.norm.vcf \
            --sample ${sample} --contig-map contig_rename.tsv --min-svlen 50 --stats norm/${sample}.${caller}.stats.tsv
    done

    # --- within-path merge of the three callers: same type, breakpoints <= 1 kb, SV >= 50 bp ---
    mkdir -p merged
    ls norm/${sample}.{syri,svimasm,cutesv}.norm.vcf > merged/${sample}.caller_vcfs.txt
    SURVIVOR merge merged/${sample}.caller_vcfs.txt 1000 1 1 0 0 50 merged/${sample}.3caller.vcf
done < paths.tsv

# Dual-evidence records (cuteSV AND SyRI or SVIM-asm) from the 38 merged VCFs
python 02_extract_dual_evidence_records.py --merged-dir merged --out results/highconf_sv.records.tsv --workers 16

# Cross-path clustering (INS: breakpoints <= 500 bp; interval SVs: endpoints <= 1 kb or reciprocal overlap >= 50%;
# both with >= 50% length similarity), TE context, carrier frequencies and saturation
python 03_cluster_population_catalog.py --records results/highconf_sv.records.tsv --repeat-bed ${REPEAT} \
    --fai ${FAI} --out-dir results/population_repeat --saturation-iterations 200 --seed 20260808

# INV / DUP / TRA from an independent cross-path clustering of SyRI rearrangements
python 04_extract_syri_rearrangements.py --norm-dir norm --out-dir results/syri_rearrangements
python 05_cluster_rearrangement_catalog.py \
    --records results/syri_rearrangements/syri_large_rearrangements.records.tsv \
    --repeat-bed ${REPEAT} --fai ${FAI} --out-dir results/syri_rearrangements_population \
    --layer-label 'SyRI large rearrangements (INV/DUP/TRA), hap38' --saturation-iterations 200 --seed 20260808

# Hybrid catalogue (dual-evidence DEL/INS + SyRI INV/DUP/TRA; 178,314 clusters) and stringent catalogue (160,396)
python 06_build_hybrid_catalog.py --tier1 results/population_repeat/sv_population_catalog.tsv \
    --syri-large results/syri_rearrangements_population/sv_population_catalog.tsv \
    --fai ${FAI} --out-dir results/final_catalog
