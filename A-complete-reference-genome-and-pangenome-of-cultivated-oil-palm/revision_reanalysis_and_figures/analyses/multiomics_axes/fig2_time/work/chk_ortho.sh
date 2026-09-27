B=${ANALYSIS_DIR}/20_results/Figure2/07_new_figure/02_comparative_orthofinder
R=$B/OrthoFinder_Results/Results_Feb05
P=$B/phylo_divtime_cafe
echo "== md5 alignment vs SpeciesTreeAlignment"; md5sum $R/WorkingDirectory/Alignments_ids/SpeciesTreeAlignment.fa $P/01_iqtree/alignment.fa $P/01_iqtree/alignment_trim.fa $P/02_mcmctree/alignment_trim.fa
echo "== untrimmed/trimmed lengths"; for f in $P/01_iqtree/alignment.fa $P/01_iqtree/alignment_trim.fa; do awk '/^>/{if(s){print length(s);exit}s="";next}{s=s$0}END{if(s)print length(s)}' $f; done
echo "== single-copy OGs"; wc -l < $R/Orthogroups/Orthogroups_SingleCopyOrthologues.txt; ls $R/Single_Copy_Orthologue_Sequences | wc -l
echo "== orthofinder log species tree lines"; grep -i -n "species tree\|single-copy\|orthogroups with\|STAG\|Alignments_ids\|minimum" $B/run_orthofinder.log | head -20
ls $R/WorkingDirectory/Alignments_ids | head; ls $R/WorkingDirectory/Alignments_ids | wc -l
grep -i -n "single\|orthogroups" $R/Log.txt | head
echo "== trimal/iqtree SIF"; grep -n "^SIF\|SIF=" $B/run_phylo_divtime_cafe5.sh
echo "== iqtree best model"; grep -n "Best-fit model\|Model of substitution" $P/01_iqtree/species_tree.iqtree | head
echo "== SpeciesTree_Gene_Duplications / statistics"; ls $R/Species_Tree $R/Comparative_Genomics_Statistics
grep -i "single-copy\|Number of orthogroups" $R/Comparative_Genomics_Statistics/Statistics_Overall.tsv | head
