# C3POa & BC1 pipeline

### example: D21


#mkdir analysis

# All pipeline transfer to /workdir/sj657/lymph_node/R2C2_pipeline/
# run preprocessing step
python3 ../../../software/C3POa/C3POa.py -r ../Data/Longread/Nanopore_basecalling/nanopore_runD21/analysis/post_basecalling_reads.fastq \
                  -o ../Data/Longread/Nanopore_basecalling/nanopore_runD21/analysis/final_run_c3poa/c3poa_output \
                  -s ../Data/Longread/Nanopore_basecalling/nanopore_runD21/analysis/final_run_c3poa/UMI_Splints1.fasta \
                  -l 500 \
                  -n 32 \
                  1> c3poa_pre.log 2> c3poa_pre.err
# postprocessing step
python3 ../../../software/C3POa/C3POa_postprocessing.py -i ../Data/Longread/Nanopore_basecalling/nanopore_runD21/analysis/final_run_c3poa/c3poa_output \
                 -a ../Data/Longread/Nanopore_basecalling/R2C2_10x_Adapters.fasta \
                 -n 32 \
                 1> c3poa_post.log 2> c3poa_post.err


# Try BC1
python3 ../../../software/BC1/bc1.py -p ../Data/Longread/Nanopore_basecalling/nanopore_runD21/analysis/final_run_c3poa/BC1_output_rc/ \
-i ../Data/Longread/Nanopore_basecalling/nanopore_runD21/analysis/final_run_c3poa/c3poa_output/UMI_Splint_1/R2C2_full_length_consensus_reads.fasta \
-o BC1_consensus -s ../Data/Longread/Nanopore_basecalling/nanopore_runD21/analysis/final_run_c3poa/c3poa_output/UMI_Splint_1/R2C2_Subreads.fastq -f -t 30 \
-u 3.0:1..NNNNNNNN.TCTTCAGCGTTCCCGAGA,3.8:9.TCTTCAGCGTTCCCGAGA.NNNNNNNNNNNNN.VV \
1> bc1_rc.log 2> bc1_rc.err
