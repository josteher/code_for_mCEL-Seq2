# code_for_mCEL-Seq2

- **Load_data.R**: Loading routine for data generated with the plate-based single cell sequencing protocol 'mCEL-Seq2'. Count tables with numbered columns 1-192 are loaded and the file name is used as identifier.
- **cresc_mapping_erv_plus_count.py**: Helper python script, that creates slurm jobs for multiple fastq files. UMI information is extracted from fastq files using umitools. The fastq files are mapped using STAR, while allowing for multimappers. Unique mapped reads for transposable elements are counted using featureCounts and the corresponding gtf file from the following website: http://hammelllab.labsites.cshl.edu/software/#TEtranscripts.
- **example_slurm_job.sh**: Example of  a slurm job that gets created by "cresc_mapping_erv_plus_count.py"
- **runs.sh**: Submit multiple slurm jobs.
