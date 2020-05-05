# ======
# Needed Modules on St. Jude HPC LSF cluster
# ------

module load node
module load toil
module load igvtools/2.3.2
module load fastqc/0.11.5
module load bowtie/1.2.2
module load macs/041014
module load ucsc/041619
module load bedtools/2.25.0
module load meme/4.11.2
module load bedops/2.4.2
module load java/1.8.0_60
module load BAM2GFF/1.1.0
module load ROSE/1.1.0
module load SICER2/1.0.1
module load samtools/1.9
module load R/3.6.1 

export R_LIBS_USER=$R_LIBS_USER:/rgs01/project_space/abrahgrp/Software_Dev_Sandbox/common/madetunj/R #local SPP

#module load phantompeakqualtools/1.2.1.1

# ======
