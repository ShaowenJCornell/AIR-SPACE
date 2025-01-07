#!/bin/bash

/programs/ont-guppy/bin/guppy_basecaller \
--recursive \
-c dna_r10.4.1_e8.2_400bps_sup.cfg \
--input_path /workdir/Data/sj657/Nanopore_R2C2_pLND10_08172023/pod5/pod5 \
--save_path /workdir/sj657/lymph_node/nanopore_runD10 \
-q 0 \
--gpu_runners_per_device 20 \
--cpu_threads_per_caller 2 \
--min_qscore 7 \
--device cuda:2,3 \
 1> pipeline.log 2> pipeline.err
