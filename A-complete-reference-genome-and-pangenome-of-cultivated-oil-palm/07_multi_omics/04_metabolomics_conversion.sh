#!/bin/bash
# Untargeted metabolomics (Waters): raw -> centroided, compressed mzML with ProteoWizard msconvert (v3.0.26121)
# Vendor peak picking at all MS levels; lock-mass calibration scans ignored.
set -euo pipefail
for raw in raw/*.raw; do
    msconvert ${raw} --mzML --zlib --64 --outdir mzML \
        --filter "peakPicking vendor msLevel=1-" --filter metadataFixer --ignoreCalibrationScans
done
# Then, separately for each ionization mode:
#   Rscript 05_metabolomics_xcms.R neg mzML_neg.list sample_manifest.tsv xcms_neg     (R v4.5.3, xcms v4.8.0)
#   Rscript 06_metabolomics_qc.R  neg xcms_neg qc_neg
