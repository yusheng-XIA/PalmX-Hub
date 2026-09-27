#!/bin/bash
# Graph-based complex SVs with SWave v1.5.2 on the completed 39-assembly minigraph-cactus graph
# Reference: FL-Hap2. Inputs per run: the clipped GFA from minigraph-cactus, the chromosome-split raw decomposed
# VCF, the reference genome and the list of the 38 non-reference assemblies.
# "SWave call" was run separately for chr01B-chr16B with default parameters (see the SWave documentation for the
# argument names of the installed version).
set -euo pipefail
for chr in $(seq -f "chr%02gB" 1 16); do
    bcftools view -r ${chr} oilpalm39.raw.vcf.gz -Oz -o raw_split/oilpalm39.raw.${chr}.vcf.gz
done
# After calling: merge chromosome-level sample VCFs and split VCFs; check that every chromosome contains all
# 38 non-reference paths. PASS records are classified from the decomposed SVTYPE and the BKPS component
# annotation as simple (one component), complex (>= 2 components, a hyperCPX event or a composite SVTYPE)
# or higher-order (>= 3 components).
