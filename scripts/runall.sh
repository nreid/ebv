#!/bin/bash

jid1=$(sbatch --parsable 01_symlinkrawdata.sh )
jid2=$(sbatch --parsable --dependency=afterok:$jid1 02_combine_lanes.sh )
jid3=$(sbatch --parsable --dependency=afterok:$jid2 03a_fastqc_small.sh )
jid4=$(sbatch --parsable --dependency=afterok:$jid3 03b_fastqc_total.sh )
jid5=$(sbatch --parsable --dependency=afterok:$jid4 04_multiqc.sh )
jid6=$(sbatch --parsable --dependency=afterok:$jid5 05a_trimmomatic_smallRNA.sh )
jid7=$(sbatch --parsable --dependency=afterok:$jid6 05b_trimmomatic_totalRNA.sh )
jid8=$(sbatch --parsable --dependency=afterok:$jid7 06a_fastqc_small_trimmed.sh )
jid9=$(sbatch --parsable --dependency=afterok:$jid8 06b_fastqc_total_trimmed.sh )
jid10=$(sbatch --parsable --dependency=afterok:$jid9 07_multiqc_trimmed.sh)
jid11=$(sbatch --parsable --dependency=afterok:$jid10 08_get_genomes.sh)
jid12=$(sbatch --parsable --dependency=afterok:$jid11 total_analysis/01_kallisto_index.sh)
jid13=$(sbatch --parsable --dependency=afterok:$jid12 total_analysis/02_kallisto_counts.sh)

total_analysis/01_kallisto_index.sh
total_analysis/02_kallisto_counts.sh
