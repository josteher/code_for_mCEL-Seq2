#! /bin/bash
#SBATCH -p bioinfo
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=20G
#SBATCH --export=ALL
#SBATCH --mail-user=herman@ie-freiburg.mpg.de
#SBATCH --mail-type=FAIL

d="/data/gruen/herman/H2AFY/JSH_23__CELSeq2_ERVmap/Data"
rn1="WT_P14_CLP_1_R1"
rn2="WT_P14_CLP_1_R2"
outdir="/data/gruen/herman/H2AFY/JSH_23__CELSeq2_ERVmap"


mkdir -p $outdir/umitools

/package/anaconda3/envs/umi_tools-0.5.1/bin/umi_tools extract --bc-pattern=NNNNNNCCCCCC                   --stdin $d/$rn1'.fastq.gz'                   --stdout $outdir/umitools/$rn1'_extracted.fastq.gz'                   --read2-in $d/$rn2'.fastq.gz'                   --read2-out=$outdir/umitools/$rn2'_extracted.fastq.gz'                   --filter-cell-barcode                   --whitelist=/data/gruen/group/herman/resources/umitools_celseq_barcodes.192.txt






mkdir -p $outdir/STAR_map

/package/STAR-2.5.3a/bin/STAR --runThreadN 4        --genomeDir /data/repository/organisms/GRCm38_ensembl/STARIndex        --readFilesIn $outdir/umitools/$rn2'_extracted.fastq.gz'        --readFilesCommand zcat        --outFilterMultimapNmax 100        --winAnchorMultimapNmax 100        --outSAMtype BAM SortedByCoordinate        --outFileNamePrefix $outdir/STAR_map/$rn2







mkdir -p $outdir/featureCounts


/package/subread-1.5.3/bin/featureCounts -t exon -g gene_id -T 4 -R BAM -a /data/gruen/herman/tools/TEToolkit/GRCm38_rmsk_TE.gtf -o $outdir/featureCounts/$rn2'_gene_assigned.txt' $outdir/STAR_map/$rn2'Aligned.sortedByCoord.out.bam'

/package/samtools-1.6.0/bin/samtools sort $outdir/featureCounts/$rn2'Aligned.sortedByCoord.out.bam.featureCounts.bam' -o $outdir/featureCounts/$rn2'_assigned_sorted.bam'
/package/samtools-1.6.0/bin/samtools index $outdir/featureCounts/$rn2'_assigned_sorted.bam'


/package/anaconda3/envs/umi_tools-0.5.1/bin/umi_tools count --per-gene --gene-tag=XT --per-cell --wide-format-cell-counts -I $outdir/featureCounts/$rn2'_assigned_sorted.bam' -S $outdir/featureCounts/$rn2'_counts.tsv'



