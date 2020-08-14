# nf-gatk

An Illumina sequencing pipeline using BWA, SAMTools, Picard, and GATK within Nextflow that is configured to be executed using Slurm.


## Contents
   1. advice.txt - file containing information that would have been very useful to me when developing this pipeline. Hopefully this can benefit future interns in some way. 
   2. nextflow - the nextflow executable.
   3. nextflow.config - configuration file for running workflows on grace/farnam clusters via Slurm.
   4. old_work - directory containing all intermediate files and final pipeline output for the sfBacteria dataset..
   5. sequence.nf - main Nextflow script.
   6. sfBacteria - dataset containing the raw reads (fastq files) and the reference genome. 
   7. sfBacteria_copy - just a copy of the dataset.

## Usage:

To execute the pipeline: 
  1. Change the directory variables at the top of sequence.nf to the corresponding paths in your local directory. 
  2. All processes are configured to be ran on Grace. If running on Farnam, change all label directives in all processes in sequence.nf from "grace" to "farnam".
  3. To have all intermediate files be saved to your work directory, uncomment the publishDir directive in each process. 
  4. If the fastq filenames do not match the pattern "*_LXXX_R[1/2]_XXX.fastq", then change the glob pattern to match the filenames on line 20 of sequence.nf that creates the input channel to the first process. 
  5. If a "vcf.gz" folder exists, make sure their isn't any previous final.vcf.gz file in that directory. If so, be sure to remove it from that folder. If you run the pipeline without doing so, the last process will be SKIPPED and a new final.vcf.gz file will NOT be made. 
  6. Finally, run the following command:
    `./nextflow run sequence.nf`
    OR 
    if you are restarting the pipeline (either after an error or change in the code) run 
    `./nextflow run sequence.nf -resume`
  7. Final output should be located in a "vcf.gz" folder in your work directory.


