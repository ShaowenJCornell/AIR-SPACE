#!/bin/bash
## cluster 3
wdir="/workdir/sj657/lymph_node/nanopore_runs_data/ONT_concensus_reads"
cd $wdir

declare -a Sample=("Mock" "D3" "D7" "D10" "D14" "D21")

for tp in "${Sample[@]}"; do
  echo $tp
  mkdir "Mixcr_results_${tp}"
  cd "Mixcr_results_${tp}"
  pwd
  # do Mixcr
  cp /workdir/sj657/lymph_node/R2C2_pipeline/mixcr_ont_preset_mod_new_param.yaml ./
  cp ../backup/Mixcr_results_${tp}/${tp}_blaze_cDNA.fastq ./
  mixcr analyze local:mixcr_ont_preset_mod_new_param \
      -Xmx40g -f --verbose \
      --species mmu \
      --not-aligned-R1 ./${tp}_unaligned.fastq \
      -Malign.parameters.saveOriginalReads=true \
      ./${tp}_blaze_cDNA.fastq \
      ./${tp}_blaze_mixcr \
      1> mixcr_mod_ps.log 2> mixcr_mod_ps.err

  mixcr exportAlignments -f -descrsR1 ./${tp}_blaze_mixcr.vdjca ./${tp}_blaze_bulk_align.tsv
  mixcr exportAlignments -f -descrsR1 -cloneId ./${tp}_blaze_mixcr.clna ./${tp}_blaze_bulk_clone.tsv
  # Try on exportAlign with nMutation
  # new version of mixcr
  mixcr exportAlignments -f -descrsR1 -cloneId -allNMutations -allNMutationsCount -nMutationsRate VRegion -chains -isotype subclass -targetSequences -isProductive VRegion \
    --impute-germline-on-export -allNLength ${tp}_blaze_mixcr.clna ./${tp}_bulk_clone_mut.tsv
  cd $wdir
done
