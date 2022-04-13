#!/usr/bin/env python

import os

fastdir = "/data/gruen/herman/H2AFY/JSH_23__CELSeq2_ERVmap/Data"
outdir  = os.getcwd()


files_R1_fastdir = list()
files_R2_fastdir = list()

# Find all files in directory
for file in os.listdir(fastdir):
    #print(file)
    if file.endswith("_R1.fastq.gz"):
        files_R1_fastdir.append(file)
    elif file.endswith("_R2.fastq.gz"):
        files_R2_fastdir.append(file)
    else:
        pass

file_names = [i.split("_R1.fastq.gz")[0] for i in files_R1_fastdir]
fastq_paths_R1 = [fastdir + "/" + i for i in file_names]
fastq_paths_R2 = [fastdir + "/" + i for i in file_names]


bash_string = """#! /bin/bash
#SBATCH -p bioinfo
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --export=ALL
#SBATCH --mail-user=herman@ie-freiburg.mpg.de
#SBATCH --mail-type=FAIL

d="%s"
rn1="%s"
rn2="%s"
outdir="%s"


mkdir -p $outdir/umitools

/package/anaconda3/envs/umi_tools-0.5.1/bin/umi_tools extract --bc-pattern=NNNNNNCCCCCC \
                  --stdin $d/$rn1'.fastq.gz' \
                  --stdout $outdir/umitools/$rn1'_extracted.fastq.gz' \
                  --read2-in $d/$rn2'.fastq.gz' \
                  --read2-out=$outdir/umitools/$rn2'_extracted.fastq.gz' \
                  --filter-cell-barcode \
                  --whitelist=/data/gruen/group/herman/resources/umitools_celseq_barcodes.192.txt






mkdir -p $outdir/STAR_map

/package/STAR-2.5.3a/bin/STAR --runThreadN 4 \
       --genomeDir /data/repository/organisms/GRCm38_ensembl/STARIndex \
       --readFilesIn $outdir/umitools/$rn2'_extracted.fastq.gz' \
       --readFilesCommand zcat \
       --outFilterMultimapNmax 100 \
       --winAnchorMultimapNmax 100 \
       --outSAMtype BAM SortedByCoordinate \
       --outFileNamePrefix $outdir/STAR_map/$rn2






mkdir -p $outdir/featureCounts


/package/subread-1.5.3/bin/featureCounts -t exon -g gene_id -T 4 -R BAM -a /data/gruen/herman/tools/TEToolkit/GRCm38_rmsk_TE.gtf -o $outdir/featureCounts/$rn2'_gene_assigned.txt' $outdir/STAR_map/$rn2'Aligned.sortedByCoord.out.bam'

/package/samtools-1.6.0/bin/samtools sort $outdir/featureCounts/$rn2'Aligned.sortedByCoord.out.bam.featureCounts.bam' -o $outdir/featureCounts/$rn2'_assigned_sorted.bam'
/package/samtools-1.6.0/bin/samtools index $outdir/featureCounts/$rn2'_assigned_sorted.bam'


/package/anaconda3/envs/umi_tools-0.5.1/bin/umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I $outdir/featureCounts/$rn2'_assigned_sorted.bam' -S $outdir/featureCounts/$rn2'_counts.tsv'



"""
# Get exported variables from bash


for i in range(len(fastq_paths_R1)):
    rn1 = file_names[i] + "_R1"
    rn2 = file_names[i] + "_R2"
    n_bash = bash_string %(fastdir, rn1 ,rn2, outdir)
    # write to file
    run_sh = file_names[i] + '_run_mapping.sh'
    f = open(run_sh,'w')
    f.write(n_bash)
    f.close()

