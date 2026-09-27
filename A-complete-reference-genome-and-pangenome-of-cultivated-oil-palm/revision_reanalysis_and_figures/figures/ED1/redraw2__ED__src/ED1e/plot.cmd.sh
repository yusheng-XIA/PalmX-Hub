#!/usr/bin/bash

min_quality=1
cool_binsize=10k
heatmap_binsize=500k


cphasing pairs2cool ../20241120-NPL2404383-P6-PBA75842-sup.pass.pairs.gz ../mixed_assembly.bp.p_utg.contigsizes 20241120-NPL2404383-P6-PBA75842-sup.pass.q${min_quality}.${cool_binsize}.cool -q ${min_quality} -bs ${cool_binsize}
cphasing plot -a ../4.scaffolding/groups.agp -m 20241120-NPL2404383-P6-PBA75842-sup.pass.q${min_quality}.${cool_binsize}.cool -o groups.q${min_quality}.${heatmap_binsize}.wg.png -bs ${heatmap_binsize} -oc
    