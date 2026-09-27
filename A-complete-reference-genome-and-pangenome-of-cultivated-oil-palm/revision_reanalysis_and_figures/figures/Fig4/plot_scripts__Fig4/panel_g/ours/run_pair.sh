#!/usr/bin/env bash
# Fig. 4g / ST19 allele classification rerun on the 2026-08-12 final annotations.
# Replicates Fig4_d_i_pan39_material33_redraw_20260808/scripts/02_run_allele_pair.sh (steps 1-5)
# with the attempt7 unique-ID fix (P__/A__ prefixes) built into the GeneTribe step.
# Differences from the original: inputs = final annotations/genomes; GMAP index in ${SCRATCH} (RAM),
# removed on exit; own conda env fig4g_g2 (python 3.9, blast 2.16.0, gmap 2025-07-31, jcvi, gffread).
set -euo pipefail

PAIR_ID="$1"; P_SAMPLE="$2"; A_SAMPLE="$3"; THREADS="${4:-14}"
FINAL=${ANALYSIS_DIR}/08_hifi_chromosome/11_final_corrected_39genomes_annotations_20260812/final
GENETRIBE=${DATA_DIR}/youzong/software/genetribe/genetribe
GET_UNMAPPED=${ANALYSIS_DIR}/10_genome_ann_contigs/13_allele/04_step03_cds/get_unmapped_pavs_fix.pl
ENV=${CLUSTER_HOME}/miniconda3/envs/fig4g_g2
BASE=${CLUSTER_WORK}/G2
WD="${BASE}/${PAIR_ID}"
OUT="${BASE}/summary/${PAIR_ID}.allele_summary.final_annot_20260923.tsv"
SHM="${SCRATCH}/fig4g_${PAIR_ID}_$$"
export PATH="${ENV}/bin:${PATH}"; export TMPDIR="${BASE}/tmp_${PAIR_ID}"; mkdir -p "${TMPDIR}"

cleanup() { [[ "${SHM}" == ${SCRATCH}/fig4g_* ]] && rm -rf -- "${SHM}"; }
trap cleanup EXIT

mkdir -p "${WD}/inputs" "${BASE}/summary"
exec > >(tee -a "${WD}/pipeline.log") 2>&1
echo "START=$(date --iso-8601=seconds) HOST=$(hostname) PAIR=${PAIR_ID} P=${P_SAMPLE} A=${A_SAMPLE} THREADS=${THREADS}"
blastp -version | head -1; gmap --version 2>&1 | grep -m1 -i version || true
python -c 'import jcvi, sys; print("jcvi", jcvi.__version__, "python", sys.version.split()[0])'
echo "GENETRIBE $(${GENETRIBE} --help 2>/dev/null | grep -m1 Version || true)"

# ---- STEP0 inputs from the final annotations (gffread -E -V, as in 11_new_pan/scripts/01_prepare_proteomes.sh)
cd "${WD}/inputs"
for side in p a; do
  s=$([[ ${side} == p ]] && echo "${P_SAMPLE}" || echo "${A_SAMPLE}")
  g=$(readlink -f "${FINAL}/01_genomes/${s}.fa"); gff="${FINAL}/02_annotations/${s}.gene.gff3"
  ln -sfn "${g}" "${side}.genome"; ln -sfn "${g}.fai" "${side}.genome.fai"; ln -sfn "${gff}" "${side}.src.gff3"
  sha256sum "${gff}" >> input_sha256.tsv
  stat -Lc '%n\t%s\t%Y' "${g}" >> genome_identity.tsv
  if [[ ! -s "${side}.pep" ]]; then
    gffread -E -V -g "${side}.genome" -x "${side}.cds" -y "${side}.pep" -o "${side}.normalized.gff3" "${side}.src.gff3" \
      2> "${side}.gffread.stderr"
  fi
  echo "INPUT ${side}=${s} mRNA=$(awk -F'\t' '$3=="mRNA"' "${gff}" | wc -l) pep=$(grep -c '^>' "${side}.pep") cds=$(grep -c '^>' "${side}.cds")"
done

# ---- STEP1 GeneTribe with haplotype-specific ID prefixes (attempt7 fix)
echo "STEP1 GeneTribe $(date --iso-8601=seconds)"
mkdir -p "${WD}/GeneTribe"; cd "${WD}/GeneTribe"
python -m jcvi.formats.gff bed --type=mRNA --key=ID ../inputs/p.src.gff3 -o pctg.bed.tmp
python -m jcvi.formats.gff bed --type=mRNA --key=ID ../inputs/a.src.gff3 -o actg.bed.tmp
sed 's/evm.TU/evm.model/g' pctg.bed.tmp > pctg.plain.bed
sed 's/evm.TU/evm.model/g' actg.bed.tmp > actg.plain.bed
awk -v p='P__' '/^>/{sub(/^>/, ">" p)} {print}' ../inputs/p.pep > pctg.fa
awk -v p='A__' '/^>/{sub(/^>/, ">" p)} {print}' ../inputs/a.pep > actg.fa
awk -v OFS='\t' -v p='P__' '{$4=p $4; print}' pctg.plain.bed > pctg.bed
awk -v OFS='\t' -v p='A__' '{$4=p $4; print}' actg.plain.bed > actg.bed
echo N > pctg.chrlist; echo N > actg.chrlist
[[ "$(wc -l < pctg.bed)" -gt 25000 && "$(wc -l < actg.bed)" -gt 25000 ]] || { echo "FATAL small BED"; exit 12; }
overlap="$(comm -12 <(cut -f4 pctg.bed | sort -u) <(cut -f4 actg.bed | sort -u) | wc -l)"
[[ "${overlap}" -eq 0 ]] || { echo "FATAL prefixed ID overlap=${overlap}"; exit 12; }
echo "PREFIXED_ID_OVERLAP=${overlap} RAW_SHARED_IDS=$(comm -12 <(cut -f4 pctg.plain.bed | sort -u) <(cut -f4 actg.plain.bed | sort -u) | wc -l)"
"${GENETRIBE}" core -l pctg -f actg -s : -n "${THREADS}" 2>&1 | tee genetribe.log
[[ -s pctg_actg.RBH ]] || { echo "FATAL empty RBH"; exit 13; }

# ---- STEP2 gene-region BLASTN (a genes vs p genes)
echo "STEP2 gene BLAST $(date --iso-8601=seconds)"
mkdir -p "${WD}/gene_blast_gene"; cd "${WD}/gene_blast_gene"
bedtools getfasta -fi ../inputs/p.genome -bed ../GeneTribe/pctg.plain.bed -fo p-gene.fa.tmp -name
bedtools getfasta -fi ../inputs/a.genome -bed ../GeneTribe/actg.plain.bed -fo a-gene.fa.tmp -name
awk -F ':' '{print $1}' p-gene.fa.tmp > p-gene.fa; awk -F ':' '{print $1}' a-gene.fa.tmp > a-gene.fa
rm -f p-gene.fa.tmp a-gene.fa.tmp
makeblastdb -in p-gene.fa -dbtype nucl -parse_seqids -hash_index -out p-gene > makeblastdb.log
blastn -db p-gene -query a-gene.fa -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore' \
  -out a-gene_to_p-gene_blast.out -num_threads "${THREADS}" -dust no
[[ -s a-gene_to_p-gene_blast.out ]] || { echo "FATAL empty gene BLAST"; exit 14; }

# ---- STEP3 a-CDS vs p-genome BLASTN -> CDS-unmapped
echo "STEP3 CDS BLAST $(date --iso-8601=seconds)"
mkdir -p "${WD}/cds_blast"; cd "${WD}/cds_blast"
ln -sfn ../inputs/p.genome p-contig.fasta; ln -sfn ../inputs/a.cds a-cds.fa
makeblastdb -in p-contig.fasta -dbtype nucl -parse_seqids -hash_index -out p-contig > makeblastdb.log
blastn -db p-contig -query a-cds.fa -evalue 1e-5 \
  -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qlen slen evalue bitscore' \
  -out a-cds_to-p-contig.blast.out -num_threads "${THREADS}" -dust no
perl "${GET_UNMAPPED}" a-cds_to-p-contig.blast.out a-cds.fa cds_unique.list
rm -f p-contig.n*

# ---- STEP4 GMAP (index in RAM, ${SCRATCH})
echo "STEP4 GMAP $(date --iso-8601=seconds)"
mkdir -p "${WD}/cds_gmap" "${SHM}/gmap_db"; cd "${WD}/cds_gmap"
gmap_build -d p-contig "${WD}/inputs/p.genome" -D "${SHM}/gmap_db" > gmap_build.log 2>&1
gmap -t "${THREADS}" -d p-contig -D "${SHM}/gmap_db" -f samse "${WD}/inputs/a.cds" > "${SHM}/a-cds_to-p-contig.sam" 2> gmap.log
[[ -s "${SHM}/a-cds_to-p-contig.sam" ]] || { echo "FATAL empty GMAP SAM"; exit 15; }
grep 'No paths found for' gmap.log | awk '{print $5}' > gmap_unique.list || true
gzip -c "${SHM}/a-cds_to-p-contig.sam" > a-cds_to-p-contig.sam.gz
rm -rf -- "${SHM}"

# ---- STEP5 classification (identical rules to 02_run_allele_pair.sh / attempt7)
echo "STEP5 classify $(date --iso-8601=seconds)"
mkdir -p "${WD}/final"; cd "${WD}/final"
sort ../cds_blast/cds_unique.list -o cds_unique.sorted
sort ../cds_gmap/gmap_unique.list -o gmap_unique.sorted
comm -12 cds_unique.sorted gmap_unique.sorted > unique_intersection.list
awk '{p=$1; a=$2; sub(/^P__/, "", p); sub(/^A__/, "", a); print a "\t" p}' ../GeneTribe/pctg_actg.RBH > rbh_pairs.tsv
awk 'NR==FNR {keys[$1 SUBSEP $2]=1; next} {key=$1 SUBSEP $2; if(key in keys && !seen[key]++) print}' \
  rbh_pairs.tsv ../gene_blast_gene/a-gene_to_p-gene_blast.out > rbh_blast.out
awk '$3==100 && $4==$11 {print $1"\t"$2}' rbh_blast.out | sort -u > homozygous.list
awk '$3<100 || $4!=$11 {print $1"\t"$2}' rbh_blast.out | sort -u > heterozygous.list
awk 'NR==FNR {hit[$1 SUBSEP $2]=1; next} {if(!(($1 SUBSEP $2) in hit)) print $1"\t"$2}' rbh_blast.out rbh_pairs.tsv > no_blast_hit.list
rbh_n=$(wc -l < rbh_pairs.tsv); homo_n=$(wc -l < homozygous.list); hetero_n=$(wc -l < heterozygous.list)
nohit_n=$(wc -l < no_blast_hit.list); cds_n=$(wc -l < cds_unique.sorted); gmap_n=$(wc -l < gmap_unique.sorted)
inter_n=$(wc -l < unique_intersection.list)
[[ $((homo_n + hetero_n + nohit_n)) -eq "${rbh_n}" ]] || { echo "FATAL RBH accounting mismatch"; exit 16; }
printf 'Sample\tp_gene\ta_gene\tRBH\tCDS_unique\tGMAP_unique\tUnique_inter\tHomozygous\tHeterozygous\tNo_hit\tp_mRNA_gff\ta_mRNA_gff\n' > "${OUT}"
printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "${PAIR_ID}" \
  "$(grep -c '^>' ../inputs/p.pep)" "$(grep -c '^>' ../inputs/a.pep)" "${rbh_n}" "${cds_n}" "${gmap_n}" "${inter_n}" \
  "${homo_n}" "${hetero_n}" "${nohit_n}" \
  "$(awk -F'\t' '$3=="mRNA"' ../inputs/p.src.gff3 | wc -l)" "$(awk -F'\t' '$3=="mRNA"' ../inputs/a.src.gff3 | wc -l)" >> "${OUT}"
cat "${OUT}"
echo "DONE=$(date --iso-8601=seconds)"
touch "${WD}/SUCCESS"
