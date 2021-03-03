#!/bin/bash 

###############################################################################
#
# important settings for sbatch, bismark alignment to execute
#
###############################################################################
# --cpus-per-task=50
# --time=96:00:00
# --mem=156g
# --gres=lscratch:800
###############################################################################

function do_cut_and_tag_analysis () {
	basedir="/data/NHLBI_BCB/Tisdale_Lab/08_Bjorg_analysis_FS/01-POGZ_cut_and_tag"
	datadir="$basedir/01-fastqs"
	outdir="$basedir/02_cut_and_tag_analysis"
	bowtie2index_hg38_ucsc="/fdb/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/"
	bowtie2index_Ecoli="/data/NHLBI_BCB/bin/bowtie2_indexes/"
	numcpus=8

	while read SAMPLE READ1 READ2; do
	    echo $SAMPLE,$READ1,$READ2
#		mkdir -p "$outdir/$SAMPLE"
		sbatch  --job-name="${SAMPLE}" \
			--partition=norm \
			--time=6:00:00 \
			--mem=256g \
			--cpus-per-task=$numcpus \
			--gres=lscratch:800 \
			--error="$outdir/$SAMPLE.POGZ.cut.and.tag.err.txt" \
			--output="$outdir/$SAMPLE.POGZ.cut.and.tag.out.txt" \
			$outdir/cut-and-tag-analysis.sh $SAMPLE $datadir $READ1 $READ2 $bowtie2index_hg38_ucsc $bowtie2index_Ecoli $outdir $numcpus
	done < "/data/NHLBI_BCB/Tisdale_Lab/08_Bjorg_analysis_FS/01-POGZ_cut_and_tag/POGZ_cut_and_tag_sampleIDs_read1_read2.txt"
}

function do_delete_sam_files () {
	basedir="/data/NHLBI_BCB/Sean_MethylSeq/23_MKJ5564"
	datadir="$basedir/01-fastqs"
	outdir="$basedir/02_methylseq_analysis_pipeline"

	while read SAMPLE READ1 READ2; do
		echo $SAMPLE,$READ1,$READ2,${SAMPLE}.bismark_pe.sam
		rm -f $outdir/${SAMPLE}/bismark_alignment/${SAMPLE}.bismark_pe.sam
	done < "/data/NHLBI_BCB/Sean_MethylSeq/23_MKJ5564/covid19_batch4_sampleIDs_read1_read2_filenames.txt"
}

do_cut_and_tag_analysis
# do_delete_sam_files
