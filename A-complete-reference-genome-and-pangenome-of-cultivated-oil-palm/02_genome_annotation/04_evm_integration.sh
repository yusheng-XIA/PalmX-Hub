#!/bin/bash
# EvidenceModeler integration of gene predictions
threads=64

# ============================================================
# Prepare weight file
# ============================================================
cat > weights.txt << 'EOF'
# Ab initio gene predictions
ABINITIO_PREDICTION	AUGUSTUS	2
ABINITIO_PREDICTION	FGENESH	2
ABINITIO_PREDICTION	EviAnn	1
ABINITIO_PREDICTION	GeneMark.hmm3	2
ABINITIO_PREDICTION	Helixer	2

# Protein alignments
PROTEIN	miniprot_protAln	5

# Transcript alignments
TRANSCRIPT	transdecoder	7
EOF

# ============================================================
# Run EvidenceModeler
# ============================================================
partition_EVM_inputs.pl \
    --genome ${sample}.nucleus.masked.fasta \
    --gene_predictions ${sample}.gene_prediction.gff3 \
    --protein_alignments ${sample}.protein_alignments.gff3 \
    --transcript_alignments ${sample}.pasa_assemblies.gff3 \
    --segmentSize 5000000 \
    --overlapSize 100000 \
    --partition_listing partitions_list.out

write_EVM_commands.pl \
    --genome ${sample}.nucleus.masked.fasta \
    --gene_predictions ${sample}.gene_prediction.gff3 \
    --protein_alignments ${sample}.protein_alignments.gff3 \
    --transcript_alignments ${sample}.pasa_assemblies.gff3 \
    --weights weights.txt \
    --output_file_name evm.out \
    --partitions partitions_list.out > evm_commands.list

execute_EVM_commands.pl evm_commands.list
recombine_EVM_partial_outputs.pl \
    --partitions partitions_list.out --output_file_name evm.out
convert_EVM_outputs_to_GFF3.pl \
    --partitions partitions_list.out \
    --output evm.out \
    --genome ${sample}.nucleus.masked.fasta

# ============================================================
# Two rounds of PASA update
# ============================================================
for round in 1 2; do
    Launch_PASA_pipeline.pl \
        -c annotCompare.config \
        -A -g ${sample}.nucleus.masked.fasta \
        -t all_transcripts.fasta \
        --CPU ${threads}
done

# ============================================================
# Extract protein and CDS sequences
# ============================================================
gff3_file_to_proteins.pl ${sample}.evm.gff3 ${sample}.nucleus.fasta prot > ${sample}.protein.faa
gff3_file_to_proteins.pl ${sample}.evm.gff3 ${sample}.nucleus.fasta CDS > ${sample}.cds.fna
