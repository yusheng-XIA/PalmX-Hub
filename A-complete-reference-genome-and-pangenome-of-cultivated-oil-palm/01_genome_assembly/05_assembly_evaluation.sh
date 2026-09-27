#!/bin/bash
# Assembly quality assessment (Supplementary Data 1-3)
set -euo pipefail
threads=64

# 1. BUSCO v6.0.0, embryophyta_odb12
busco -m genome -i ${sample}.fa -l embryophyta_odb12 -o ${sample}_busco -c ${threads} --offline

# 2. Telomeres (tidk v0.2.65): TTTAGGG repeats counted in 10-kb windows; a chromosome end is telomere-positive
#    when at least one window within the terminal 50 kb contains >= 3 copies (same rule for EG11 and EO12)
tidk search -s TTTAGGG -w 10000 -o ${sample} -d tidk_${sample} ${sample}.fa
awk -F '\t' 'NR==FNR{len[$1]=$2; next}
    FNR>1 { n = $3 + $4; if (n < 3) next
            if ($2 <= 50000) left[$1] = 1
            if ($2 >= len[$1] - 50000) right[$1] = 1 }
    END { for (c in len) print c "\t" (c in left ? 1 : 0) "\t" (c in right ? 1 : 0) }' \
    ${sample}.fa.fai tidk_${sample}/${sample}_telomeric_repeat_windows.tsv | sort > ${sample}.telomere_ends.tsv

# 3. LTR Assembly Index (LTR_retriever v3.0.2; maximum intact-LTR length 15 kb)
gt suffixerator -db ${sample}.fa -indexname ${sample} -tis -suf -lcp -des -ssp -sds -dna
gt ltrharvest -index ${sample} -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 \
    -similar 85 -vic 10 -seed 20 -seqids yes > ${sample}.harvest.scn
LTR_FINDER_parallel -seq ${sample}.fa -threads ${threads} -harvest_out
cat ${sample}.harvest.scn ${sample}.fa.finder.combine.scn > ${sample}.rawLTR.scn
LTR_retriever -genome ${sample}.fa -inharvest ${sample}.rawLTR.scn -threads ${threads} -maxlenltr 15000
LAI -genome ${sample}.fa -intact ${sample}.fa.pass.list -all ${sample}.fa.out

# 4. QV and k-mer completeness (Merqury v1.3, HiFi meryl database)
k=$(best_k.sh 1800000000 | tail -1 | cut -d. -f1)
meryl k=${k} count output ${sample}.hifi.meryl ${sample}.hifi.fq.gz
merqury.sh ${sample}.hifi.meryl ${sample}.hap1.fa ${sample}.hap2.fa ${sample}_merqury

# 5. Long-read remapping: mapping rate and covered fraction
minimap2 -ax map-hifi -t ${threads} ${sample}.fa ${sample}.hifi.fq.gz | samtools sort -@ ${threads} -o ${sample}.hifi.bam
minimap2 -ax map-ont  -t ${threads} ${sample}.fa ${sample}.ont.fq.gz  | samtools sort -@ ${threads} -o ${sample}.ont.bam
for b in hifi ont; do
    samtools index ${sample}.${b}.bam
    samtools flagstat ${sample}.${b}.bam > ${sample}.${b}.flagstat
    samtools coverage ${sample}.${b}.bam > ${sample}.${b}.coverage.tsv
done
