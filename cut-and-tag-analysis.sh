#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

SAMPLE=$1
DATAPATH=$2
READ1=$3
READ2=$4
bowtie2index_hg38_ucsc=$5
bowtie2index_Ecoli=$6
outdir=$7
numcpus=$8

trimthreads=4
alignthreads=2

# module load trimmomatic
# module load trimgalore
module load fastqc
module load bowtie
module load picard
module load samtools
module load bedtools
module load R
module load homer
module load macs

function do_fastqc () {
	date
	########################################################################################################################
	trimmed=$1

	if [[ "$trimmed" == "no" ]]; 
	then
		DATAPATH=$DATAPATH
		READ1=$READ1
		READ2=$READ2
		out_dir="fastqc"
		#statements
	else
		# if trimming, change DATAPATH (scratch OR saving trimmed fastq files in local directory)
#		DATAPATH="/lscratch/${SLURM_JOBID}"
		DATAPATH="$outdir/$SAMPLE/TrimGalore"

		# if trimming, change $READ1 and $READ2, using Trimmomatic
		#	READ1="${SAMPLE}_1P.fastq.gz"
		#	READ2="${SAMPLE}_2P.fastq.gz"

		# if trimming, change $READ1 and $READ2, using TrimGalore
		READ1="${SAMPLE}_val_1.fq.gz"
		READ2="${SAMPLE}_val_2.fq.gz"

		# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
		out_dir="fastqc_posttrimming"
	fi

#	fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

	mkdir -p $outdir/$SAMPLE/$out_dir

	fastqc -o "$outdir/$SAMPLE/$out_dir"  \
	--nogroup \
	"$DATAPATH/$READ1"  \
	"$DATAPATH/$READ2"  \
	|| fail "fastqc failed"

	echo "fastqc done"
	########################################################################################################################
	date
}

function do_get_fastqc_stats () {
	date
	########################################################################################################################
	trimmed=$1

	if [[ "$trimmed" == "no" ]]; 
	then
		fastqc_READ1=`echo $READ1 | sed 's/.fastq.gz/_fastqc.zip/g'`
		fastqc_READ2=`echo $READ2 | sed 's/.fastq.gz/_fastqc.zip/g'`
		fastqc_READ1_dir=`echo $READ1 | sed 's/.fastq.gz/_fastqc/g'`
		out_dir="fastqc"
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ1 -d $outdir/$SAMPLE/$out_dir/
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ2 -d $outdir/$SAMPLE/$out_dir/
		total_sequences=`grep "Total Sequences" $outdir/$SAMPLE/$out_dir/$fastqc_READ1_dir/fastqc_data.txt`
		total_sequences_num=`echo $total_sequences | cut -d" " -f3`
		echo -e $SAMPLE'\t'$total_sequences_num >> "$outdir/total_sequences_num_notrim.txt"
	else
		# if trimming, change $READ1 and $READ2, using TrimGalore
		READ1="${SAMPLE}_val_1.fq.gz"
		READ2="${SAMPLE}_val_2.fq.gz"
		fastqc_READ1=`echo $READ1 | sed 's/.fq.gz/_fastqc.zip/g'`
		fastqc_READ2=`echo $READ2 | sed 's/.fq.gz/_fastqc.zip/g'`
		fastqc_READ1_dir=`echo $READ1 | sed 's/.fq.gz/_fastqc/g'`
		# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
		out_dir="fastqc_posttrimming"
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ1 -d $outdir/$SAMPLE/$out_dir/
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ2 -d $outdir/$SAMPLE/$out_dir/
		total_sequences=`grep "Total Sequences" $outdir/$SAMPLE/$out_dir/$fastqc_READ1_dir/fastqc_data.txt`
		total_sequences_num=`echo $total_sequences | cut -d" " -f3`
		echo -e $SAMPLE'\t'$total_sequences_num >> "$outdir/total_sequences_num.txt"
	fi
	########################################################################################################################
	date
}

function do_trimmomatic () {
	date
	########################################################################################################################
	java -jar $TRIMMOJAR PE \
	            -threads $numcpus \
	            "$DATAPATH/$READ1" \
	            "$DATAPATH/$READ2" \
	 	    -baseout "/lscratch/${SLURM_JOBID}/${SAMPLE}.fastq.gz" \
	           ILLUMINACLIP:"/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
	           MINLEN:50 \
		   HEADCROP:10
	
	# if trimming bases from the start of reads, change HEADCROP
	# HEADCROP:10	\
	echo "trimmomatic done"
	########################################################################################################################
	date
}

function do_trimgalore () {
	date
	########################################################################################################################
	out_dir="TrimGalore"
	mkdir -p $outdir/$SAMPLE/$out_dir

	trim_galore 	--paired \
			--cores $trimthreads \
			--basename ${SAMPLE} \
			--output_dir "$outdir/$SAMPLE/$out_dir" \
			--clip_R1 10 \
			--clip_R2 10 \
			--three_prime_clip_R1 10 \
			--three_prime_clip_R2 10 \
			--length 50 \
			--gzip \
			"$DATAPATH/$READ1" "$DATAPATH/$READ2"

	echo "TrimGalore done"
	########################################################################################################################
	date
#	--output_dir "/lscratch/${SLURM_JOBID}/" \
}

function do_bowtie2_alignment () {
	date
	########################################################################################################################
	trimmed=$1

	if [[ "$trimmed" == "no" ]]; 
	then
		DATAPATH=$DATAPATH
		READ1=$READ1
		READ2=$READ2
		out_dir="bowtie2"
		#statements
	else
		echo "processing TRIMMED fastq files..."
		# if trimming, change DATAPATH (scratch OR saving trimmed fastq files in local directory)
#		DATAPATH="/lscratch/${SLURM_JOBID}"
#		DATAPATH="$outdir/$SAMPLE/TrimGalore"

		# if trimming, change $READ1 and $READ2, using Trimmomatic
#		READ1="${SAMPLE}_1P.fastq.gz"
#		READ2="${SAMPLE}_2P.fastq.gz"

		# if trimming, change $READ1 and $READ2, using TrimGalore
#		READ1="${SAMPLE}_val_1.fq.gz"
#		READ2="${SAMPLE}_val_2.fq.gz"

		# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
#		out_dir="bowtie2_posttrimming"
	fi

	mkdir -p $outdir/$SAMPLE/$out_dir

	bowtie2 --end-to-end \
	--very-sensitive \
	--no-mixed \
	--no-discordant \
	--phred33 \
	-I 10 \
	-X 700 \
	-p ${numcpus} \
	-x "${bowtie2index_hg38_ucsc}/genome" \
	-1 $DATAPATH/$READ1 \
	-2 $DATAPATH/$READ2 \
	-S "$outdir/$SAMPLE/$out_dir/${SAMPLE}_bowtie2.sam"

	echo "bowtie2 done"
	########################################################################################################################
	date

}

function do_bowtie2_build_index () {
	date
	########################################################################################################################

	bowtie2-build /data/NHLBI_BCB/bin/bowtie2_indexes/Ecoli.sequence.fasta /data/NHLBI_BCB/bin/bowtie2_indexes/Ecoli

	########################################################################################################################
	date
}

function do_bowtie2_alignment_to_spikein_genome () {
	date
	########################################################################################################################

	##== linux command ==##
#	spikeInRef="/shared/ngs/illumina/henikoff/Bowtie2/Ecoli"
#	chromSize="/fh/fast/gottardo_r/yezheng_working/SupplementaryData/hg38/chromSize/hg38.chrom.size"

	trimmed=$1

	if [[ "$trimmed" == "no" ]]; 
	then
		DATAPATH=$DATAPATH
		READ1=$READ1
		READ2=$READ2
		out_dir="bowtie2_Ecoli"
		#statements
	else
		echo "processing TRIMMED fastq files..."
		# if trimming, change DATAPATH (scratch OR saving trimmed fastq files in local directory)
#		DATAPATH="/lscratch/${SLURM_JOBID}"
#		DATAPATH="$outdir/$SAMPLE/TrimGalore"

		# if trimming, change $READ1 and $READ2, using Trimmomatic
#		READ1="${SAMPLE}_1P.fastq.gz"
#		READ2="${SAMPLE}_2P.fastq.gz"

		# if trimming, change $READ1 and $READ2, using TrimGalore
#		READ1="${SAMPLE}_val_1.fq.gz"
#		READ2="${SAMPLE}_val_2.fq.gz"

		# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
#		out_dir="bowtie2_Ecoli_posttrimming"
	fi

	mkdir -p $outdir/$SAMPLE/$out_dir

	## bowtie2-build path/to/Ecoli/fasta/Ecoli.fa /path/to/bowtie2Index/Ecoli
	bowtie2 --end-to-end \
	--very-sensitive \
	--no-mixed \
	--no-discordant \
	--phred33 \
	-I 10 \
	-X 700 \
	-p ${numcpus} \
	-x "${bowtie2index_Ecoli}/Ecoli" \
	-1 $DATAPATH/$READ1 \
	-2 $DATAPATH/$READ2 \
	-S "$outdir/$SAMPLE/$out_dir/${SAMPLE}_bowtie2_Ecoli.sam"

#	seqDepthDouble=`samtools view -F 0x04 seqDepth=$((seqDepthDouble/2))
#	echo $seqDepth >$projPath/alignment/sam/bowtie2_summary/${histName}_bowtie2_spikeIn.seqDept
	date
	########################################################################################################################
}

function do_picard_remove_PCR_duplicates () {
	date
	########################################################################################################################

	out_dir="picard_summary"

	##== linux command ==##
	## depending on how you load picard and your server environment, the picardCMD can be different. Adjust accordingly.
#	picardCMD="java -jar picard.jar"
	mkdir -p $outdir/$SAMPLE/$out_dir

	TMP_DIR=/lscratch/$SLURM_JOBID

	## Sort by coordinate
	java -Xmx156g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar SortSam I=$outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.sam O=$TMP_DIR/${SAMPLE}_bowtie2.sorted.sam SORT_ORDER=coordinate

	## mark duplicates
	java -Xmx156g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates I=$TMP_DIR/${SAMPLE}_bowtie2.sorted.sam O=$TMP_DIR/${SAMPLE}_bowtie2.sorted.dupMarked.sam METRICS_FILE=$outdir/$SAMPLE/$out_dir/${SAMPLE}_picard.dupMark.txt

	## remove duplicates
	java -Xmx156g -XX:ParallelGCThreads=5 -jar $PICARDJARPATH/picard.jar MarkDuplicates I=$TMP_DIR/${SAMPLE}_bowtie2.sorted.sam O=$outdir/$SAMPLE/$out_dir/${SAMPLE}_bowtie2.sorted.rmDup.sam REMOVE_DUPLICATES=true METRICS_FILE=$outdir/$SAMPLE/$out_dir/${SAMPLE}_picard.rmDup.txt

	########################################################################################################################
	date
}

function do_calculate_fragment_length () {
	date
	########################################################################################################################

	out_dir="fragment_length_deduplicated"

	##== linux command ==##
	mkdir -p $outdir/$SAMPLE/$out_dir

	## Extract the 9th column from the alignment sam file which is the fragment length
	samtools view -F 0x04 $outdir/$SAMPLE/picard_summary/${SAMPLE}_bowtie2.sorted.rmDup.sam | awk -F'\t' 'function abs(x){return ((x < 0.0) ? -x : x)} {print abs($9)}' | sort | uniq -c | awk -v OFS="\t" '{print $2, $1/2}' > $outdir/$SAMPLE/$out_dir/${SAMPLE}_fragmentLen.txt

	########################################################################################################################
	date

}

function do_file_format_conversion () {

	TMP_DIR=/lscratch/$SLURM_JOBID

	##== linux command ==##
	## Sort the SAM file and convert to BAM
	samtools sort -T $TMP_DIR/${SAMPLE} -o $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.bam $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.sam
	## index BAM 
	samtools index $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.bam
	## Filter and keep the mapped read pairs
	samtools view -bS -F 0x04 $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.bam > $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.bam
	## Sort mapped bam file by read, important for next step
	samtools sort -T $TMP_DIR/${SAMPLE} -n -o $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.sortedbyreads.mapped.bam $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.bam
	## Convert into bed file format
	bedtools bamtobed -i $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.sortedbyreads.mapped.bam -bedpe > $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.bed
	## Keep the read pairs that are on the same chromosome and fragment length less than 1000bp.
	awk '$1==$4 && $6-$2 < 1000 {print $0}' $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.bed > $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.clean.bed
	## Only extract the fragment related columns
	cut -f 1,2,6 $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.clean.bed | sort -k1,1 -k2,2n -k3,3n  > $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.clean.fragments.bed

}

function do_reproducibility () {

	##== linux command ==##
	## We use the mid point of each fragment to infer which 500bp bins does this fragment belong to.
	binLen=500
	awk -v w=$binLen '{print $1, int(($2 + $3)/(2*w))*w + w/2}' $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.clean.fragments.bed | sort -k1,1V -k2,2n | uniq -c | awk -v OFS="\t" '{print $2, $3, $1}' |  sort -k1,1V -k2,2n  > $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.clean.fragmentsCount.bin$binLen.bed

}

function do_spike_in_calibration () {

	out_dir="bedgraph"

	seqDepthDouble=`samtools view -F 0x04 $outdir/$SAMPLE/bowtie2_Ecoli/${SAMPLE}_bowtie2_Ecoli.sam | wc -l`
	seqDepth=$((seqDepthDouble/2))
	echo $seqDepth > $outdir/$SAMPLE/bowtie2_Ecoli/${SAMPLE}_bowtie2_spikeIn.seqDepth
	chromSize="/data/NHLBI_BCB/Tisdale_Lab/08_Bjorg_analysis_FS/hg38_chromosome_sizes.txt"

	##== linux command ==##
#	if [[ "$seqDepth" -gt "1" ]]; then
	    
#	    mkdir -p $outdir/$SAMPLE/$out_dir

	    scale_factor=`echo "10000 / $seqDepth" | bc -l`
	    echo "Scaling factor for $SAMPLE is: $scale_factor!"
	    bedtools genomecov -bg -scale $scale_factor -i $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.clean.fragments.bed -g $chromSize > $outdir/$SAMPLE/$out_dir/${SAMPLE}_bowtie2.fragments.normalized.bedgraph
	
#	fi
	bedtools genomecov -bg -i $outdir/$SAMPLE/bowtie2/${SAMPLE}_bowtie2.mapped.clean.fragments.bed -g $chromSize > $outdir/$SAMPLE/$out_dir/${SAMPLE}_bowtie2.fragments.bedgraph

}

function do_run_SEACR () {

	out_dir="SEACR_peaks"
	input_dir="bedgraph"

	mkdir -p $outdir/$SAMPLE/$out_dir

	bash /data/NHLBI_BCB/bin/SEACR/SEACR_1.3.sh \
		$outdir/$SAMPLE/$input_dir/${SAMPLE}_bowtie2.fragments.bedgraph \
		0.01 \
		"non" \
		"stringent" \
		$outdir/$SAMPLE/$out_dir/$SAMPLE

echo -e "

SEACR: Sparse Enrichment Analysis for CUT&RUN
	
	Usage: bash SEACR_1.3.sh <experimental bedgraph>.bg [<control bedgraph>.bg | <FDR threshold>] [norm | non] [relaxed | stringent] output prefix
	
	Description of input fields:
	
	Field 1: Target data bedgraph file in UCSC bedgraph format (https://genome.ucsc.edu/goldenpath/help/bedgraph.html) that omits regions containing 0 signal.
	
	Field 2: Control (IgG) data bedgraph file to generate an empirical threshold for peak calling. Alternatively, a numeric threshold n between 0 and 1 returns the top n fraction of peaks based on total signal within peaks.
	
	Field 3: “norm” denotes normalization of control to target data, “non” skips this behavior. norm is recommended unless experimental and control data are already rigorously normalized to each other (e.g. via spike-in).
		
	Field 4: “relaxed” uses a total signal threshold between the knee and peak of the total signal curve, and corresponds to the “relaxed” mode described in the text, whereas “stringent” uses the peak of the curve, and corresponds to “stringent” mode.
	
	Field 5: Output prefix
	
	Output file:
	<output prefix>.auc.threshold.merge.bed (Bed file of enriched regions)
	
	Output data structure: 
	
	<chr>	<start>	<end>	<AUC>	<max signal>	<max signal region>
	
	Description of output fields:
	Field 1: Chromosome
	
	Field 2: Start coordinate
	
	Field 3: End coordinate
	
	Field 4: Total signal contained within denoted coordinates
	
	Field 5: Maximum bedgraph signal attained at any base pair within denoted coordinates
	
	Field 6: Region representing the farthest upstream and farthest downstream bases within the denoted coordinates that are represented by the maximum bedgraph signal
	
	Examples:
	bash SEACR_1.3.sh target.bedgraph IgG.bedgraph norm stringent output
	Calls enriched regions in target data using normalized IgG control track with stringent threshold
	
	bash SEACR_1.3.sh target.bedgraph IgG.bedgraph non relaxed output
	Calls enriched regions in target data using non-normalized IgG control track with relaxed threshold
	bash SEACR_1.3.sh target.bedgraph 0.01 non stringent output
	Calls enriched regions in target data by selecting the top 1% of regions by area under the curve (AUC)

" > /dev/null

}

function do_HOMER_annotation () {

	out_dir="MACS2_peaks_minlength25_HOMER"
	input_dir="MACS2_peaks_minlength25"
#	bedgraph_dir="bedgraph"

	mkdir -p $outdir/$SAMPLE/$out_dir

	cat $outdir/$SAMPLE/$input_dir/${SAMPLE}_peaks.narrowPeak | sed '1,1d' | awk '{print $0"\t"$3-$2}' > $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.input.bed

	annotatePeaks.pl $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.input.bed hg38 -bedGraph $outdir/$SAMPLE/$input_dir/${SAMPLE}_treat_pileup.bdg > $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.output.bed

	awk 'NR==FNR{A[$4]=$7"\t"$8"\t"$9"\t"$10"\t"$11;next}{print$0"\t"A[$1]}' $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.input.bed $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.output.bed | sed '1,1d' > $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.output.v2.bed

	cat $outdir/header_HOMER.txt $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.output.v2.bed > $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.output.v2.bed.temp

	mv -f $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.output.v2.bed.temp $outdir/$SAMPLE/$out_dir/${SAMPLE}_peaks.narrowPeak.homer.output.v2.txt

	echo -e "
		Usage: annotatePeaks.pl <peak file | tss> <genome version>  [additional options...]

		Available Genomes (required argument): (name,org,directory,default promoter set)
			rn4	rat	/usr/local/apps/homer/4.11.1/.//data/genomes/rn4/	default
			mm9	mouse	/usr/local/apps/homer/4.11.1/.//data/genomes/mm9/	default
			sacCer2	yeast	/usr/local/apps/homer/4.11.1/.//data/genomes/sacCer2/	default
			galGal4	chicken	/usr/local/apps/homer/4.11.1/.//data/genomes/galGal4/	default
			dm3	fly	/usr/local/apps/homer/4.11.1/.//data/genomes/dm3/	default
			ASM294v1	pombe	/usr/local/apps/homer/4.11.1/.//data/genomes/ASM294v1/	default
			ce6	worm	/usr/local/apps/homer/4.11.1/.//data/genomes/ce6/	default
			rn5	rat	/usr/local/apps/homer/4.11.1/.//data/genomes/rn5/	default
			rheMac8	rhesus	/usr/local/apps/homer/4.11.1/.//data/genomes/rheMac8/	default
			danRer10	zebrafish	/usr/local/apps/homer/4.11.1/.//data/genomes/danRer10/	default
			rn6	rat	/usr/local/apps/homer/4.11.1/.//data/genomes/rn6/	default
			mm8	mouse	/usr/local/apps/homer/4.11.1/.//data/genomes/mm8/	default
			mm10	mouse	/usr/local/apps/homer/4.11.1/.//data/genomes/mm10/	default
			hg38	human	/usr/local/apps/homer/4.11.1/.//data/genomes/hg38/	default
			hg19	human	/usr/local/apps/homer/4.11.1/.//data/genomes/hg19/	default
			hg18	human	/usr/local/apps/homer/4.11.1/.//data/genomes/hg18/	default
			ce10	worm	/usr/local/apps/homer/4.11.1/.//data/genomes/ce10/	default
				-- or --
			Custom: provide the path to genome FASTA files (directory or single file)
			If no genome is available, specify 'none'.
			If using FASTA file or none, may want to specify '-organism <...>'

		User defined annotation files (default is UCSC refGene annotation):
			annotatePeaks.pl accepts GTF (gene transfer formatted) files to annotate positions relative
			to custom annotations, such as those from de novo transcript discovery or Gencode.
			-gtf <gtf format file> (Use -gff and -gff3 if appropriate, but GTF is better)
			-gid (by default the GTF file is processed by transcript_id, use this option for gene_id)
			-ann <custom homer annotation file> (created by assignGenomeAnnotation, see website)

		Peak vs. tss/tts/rna mode (works with custom GTF file):
			If the first argument is tss (i.e. annotatePeaks.pl tss hg18 ...) then a TSS centric
			analysis will be carried out.  Tag counts and motifs will be found relative to the TSS.
			(no position file needed) [tts now works too - e.g. 3' end of gene]
			[rna specifies gene bodies, will automaticall set -size given]
			NOTE: The default TSS peak size is 4000 bp, i.e. +/- 2kb (change with -size option)
			-list <gene id list> (subset of genes to perform analysis [unigene, gene id, accession,
				 probe, etc.], default = all promoters)
			-cTSS <promoter position file i.e. peak file> (should be centered on TSS)

		Primary Annotation Options:
			-mask (Masked repeats, can also add 'r' to end of genome name)
			-m <motif file 1> [motif file 2] ... (list of motifs to find in peaks)
				-mscore (reports the highest log-odds score within the peak)
				-nmotifs (reports the number of motifs per peak)
				-mdist (reports distance to closest motif)
				-mfasta <filename> (reports sites in a fasta file - for building new motifs)
				-fm <motif file 1> [motif file 2] (list of motifs to filter from above)
				-rmrevopp <#> (only count sites found within <#> on both strands once, i.e. palindromic)
				-matrix <prefix> (outputs a motif co-occurrence files:
					prefix.count.matrix.txt - number of peaks with motif co-occurrence
					prefix.ratio.matrix.txt - ratio of observed vs. expected  co-occurrence
					prefix.logPvalue.matrix.txt - co-occurrence enrichment
					prefix.stats.txt - table of pair-wise motif co-occurrence statistics
					additional options:
					-matrixMinDist <#> (minimum distance between motif pairs - to avoid overlap, default: 4)
					-matrixMaxDist <#> (maximum distance between motif pairs)
				-mbed <filename> (Output motif positions to a BED file to load at UCSC (or -mpeak))
				-mlogic <filename> (will output stats on common motif orientations)
			-d <tag directory 1> [tag directory 2] ... (list of experiment directories to show
				tag counts for) NOTE: -dfile <file> where file is a list of directories in first column
			-bedGraph <bedGraph file 1> [bedGraph file 2] ... (read coverage counts from bedGraph files)
			-wig <wiggle file 1> [wiggle file 2] ... (read coverage counts from wiggle files)
			-p <peak file> [peak file 2] ... (to find nearest peaks)
				-pdist to report only distance (-pdist2 gives directional distance)
				-pcount to report number of peaks within region
			-vcf <VCF file> (annotate peaks with genetic variation infomation, one col per individual)
				-editDistance (Computes the # bp changes relative to reference)
				-individuals <name1> [name2] ... (restrict analysis to these individuals)
			-gene <data file> ... (Adds additional data to result based on the closest gene.
				This is useful for adding gene expression data.  The file must have a header,
				and the first column must be a GeneID, Accession number, etc.  If the peak
				cannot be mapped to data in the file then the entry will be left empty.
			-go <output directory> (perform GO analysis using genes near peaks)
			-genomeOntology <output directory> (perform genomeOntology analysis on peaks)
				-gsize <#> (Genome size for genomeOntology analysis, default: 2e9)

		Annotation vs. Histogram mode:
			-hist <bin size in bp> (i.e 1, 2, 5, 10, 20, 50, 100 etc.)
			The -hist option can be used to generate histograms of position dependent features relative
			to the center of peaks.  This is primarily meant to be used with -d and -m options to map
			distribution of motifs and ChIP-Seq tags.  For ChIP-Seq peaks for a Transcription factor
			you might want to use the -center option (below) to center peaks on the known motif
			** If using -size given, histogram will be scaled to each region (i.e. 0-100%), with
			the -hist parameter being the number of bins to divide each region into.
				Histogram Mode specific Options:
				-nuc (calculated mononucleotide frequencies at each position,
					Will report by default if extracting sequence for other purposes like motifs)
				-di (calculated dinucleotide frequencies at each position)
				-histNorm <#> (normalize the total tag count for each region to 1, where <#> is the
					minimum tag total per region - use to avoid tag spikes from low coverage
				-ghist (outputs profiles for each gene, for peak shape clustering)
				-rm <#> (remove occurrences of same motif that occur within # bp)

		Peak Centering: (other options are ignored)
			-center <motif file> (This will re-center peaks on the specified motif, or remove peak
				if there is no motif in the peak.  ONLY recentering will be performed, and all other
				options will be ignored.  This will output a new peak file that can then be reanalyzed
				to reveal fine-grain structure in peaks (It is advised to use -size < 200) with this
				to keep peaks from moving too far (-mirror flips the position)
			-multi (returns genomic positions of all sites instead of just the closest to center)

		Genome comparisons (need genome & liftOver)
			-cmpGenome <genome1> [genome2] (Genomes to compare for sequence/motifs)
			-cmpLiftover <liftover1> [genome2] (Genomes to compare for sequence/motifs)

		Normalization options:
			-fpkm (normalize read counts to million reads or fragments per kilobase mapped)
			-raw (do not adjust the tag counts based on total tags sequenced, -noadj works too)
			-norm <#> (normalize tags to this tag count, default=1e7, 0=average tag count in all directories)
			-normLength <#> (Fragment length to normlize to for experiments with different lens, def: 100)
			-log (output tag counts as log2(x+1+rand) values - for scatter plots)
			-sqrt (output tag counts as sqrt(x+rand) values - for scatter plots)
			-ratio (process tag values as ratios - i.e. chip-seq, or mCpG/CpG)

		Advanced normalization options: (-rlog and -vst require R and DESeq2 to be installed)
			-rlog (quantile/variance normalization on reported genes using DESeq2 rlog funcition, needs R)
			-vst (quantile/variance normalization on reported genes using DESeq2 vst function, needs R)

		Advanced Options:
			-len <#> / -fragLength <#> (Fragment length, default=auto, might want to set to 1 for 5'RNA)
			-size <#> (Peak size[from center of peak], default=inferred from peak file)
				-size #,# (i.e. -size -10,50 count tags from -10 bp to +50 bp from center)
				-size given (count tags etc. using the actual regions - for variable length regions)
			-strand <+|-|both> (Count tags on specific strands relative to peak, default: both)
			-pc <#> (maximum number of tags to count per bp, default=0 [no maximum], -tbp <#> works too)
			-CpG (Calculate CpG/GC content)
			-nfr (report nuclesome free region scores instead of tag counts, also -nfrSize <#>)
			-norevopp (do not search for motifs on the opposite strand [works with -center too])
			-gwasCatalog <gwasCatalog file from UCSC> (list overlapping GWAS risk SNPs)
			-pdist (only report distance to nearest peak using -p, not peak name)
			-map <mapping file> (mapping between peak IDs and promoter IDs, overrides closest assignment)
			-noann, -nogene (skip genome annotation step, skip TSS annotation)
			-homer1/-homer2 (by default, the new version of homer [-homer2] is used for finding motifs)
			-cpu <#> (Number of processors to use when possible - only some parts utilize multiple cores)
			-noblanks (remove peaks/rows with missing data)

	" > /dev/null

}

function do_HOMER_motif_finding	() {

	out_dir="MACS2_peaks_minlength25_HOMER_motif_finding_mask"
	input_dir="MACS2_peaks_minlength25_HOMER"

#	cat $outdir/$SAMPLE/$input_dir/${SAMPLE}_peaks.narrowPeak.homer.input.bed | awk '{if($11<=120){print $0}}' > $outdir/$SAMPLE/$input_dir/${SAMPLE}_peaks.narrowPeak.homer.motif.input.bed

	mkdir -p $outdir/$SAMPLE/$out_dir

	findMotifsGenome.pl $outdir/$SAMPLE/$input_dir/${SAMPLE}_peaks.narrowPeak.homer.input.bed hg38 $outdir/$SAMPLE/$out_dir/ -size -100,100 -mask
#	findMotifsGenome.pl $outdir/$SAMPLE/$input_dir/${SAMPLE}_peaks.narrowPeak.homer.motif.input.bed hg38 $outdir/$SAMPLE/$out_dir/ -size -100,100

	echo -e "

	Program will find de novo and known motifs in regions in the genome

	Usage: findMotifsGenome.pl <pos file> <genome> <output directory> [additional options]
	Example: findMotifsGenome.pl peaks.txt mm8r peakAnalysis -size 200 -len 8

	Possible Genomes:
		sacCer2	yeast
		hg38	human
		ce6	worm
		ce10	worm
		mm10	mouse
		danRer10	zebrafish
		mm9	mouse
		rheMac10	rhesus
		hg19	human
		hg18	human
		rn6	rat
		dm3	fly
		ASM294v1	pombe
		galGal4	chicken
		rheMac8	rhesus
		mm8	mouse
		rn5	rat
		rn4	rat
			-- or --
		Custom: provide the path to genome FASTA files (directory or single file)
			Heads up: will create the directory preparsed/ in same location.

	Basic options:
		-mask (mask repeats/lower case sequence, can also add 'r' to genome, i.e. mm9r)
		-bg <background position file> (genomic positions to be used as background, default=automatic)
			removes background positions overlapping with target positions unless -keepOverlappingBg is used
			-chopify (chop up large background regions to the avg size of target regions)
		-len <#>[,<#>,<#>...] (motif length, default=8,10,12) [NOTE: values greater 12 may cause the program
			to run out of memory - in these cases decrease the number of sequences analyzed (-N),
			or try analyzing shorter sequence regions (i.e. -size 100)]
		-size <#> (fragment size to use for motif finding, default=200)
			-size <#,#> (i.e. -size -100,50 will get sequences from -100 to +50 relative from center)
			-size given (uses the exact regions you give it)
		-S <#> (Number of motifs to optimize, default: 25)
		-mis <#> (global optimization: searches for strings with # mismatches, default: 2)
		-norevopp (don't search reverse strand for motifs)
		-nomotif (don't search for de novo motif enrichment)
		-rna (output RNA motif logos and compare to RNA motif database, automatically sets -norevopp)

	Scanning sequence for motifs
		-find <motif file> (This will cause the program to only scan for motifs)

	Known Motif Options/Visualization
		-mset <vertebrates|insects|worms|plants|yeast|all> (check against motif collects, default: auto)
		-basic (just visualize de novo motifs, don't check similarity with known motifs)
		-bits (scale sequence logos by information content, default: doesn't scale)
		-nocheck (don't search for de novo vs. known motif similarity)
		-mcheck <motif file> (known motifs to check against de novo motifs,
		-float (allow adjustment of the degeneracy threshold for known motifs to improve p-value[dangerous])
		-noknown (don't search for known motif enrichment, default: -known)
		-mknown <motif file> (known motifs to check for enrichment,
		-nofacts (omit humor)
		-seqlogo (use weblogo/seqlogo/ghostscript to generate logos, default uses SVG now)

	Sequence normalization options:
		-gc (use GC% for sequence content normalization, now the default)
		-cpg (use CpG% instead of GC% for sequence content normalization)
		-noweight (no CG correction)
		Also -nlen <#>, -olen <#>, see homer2 section below.

	Advanced options:
		-h (use hypergeometric for p-values, binomial is default)
		-N <#> (Number of sequences to use for motif finding, default=max(50k, 2x input)
		-local <#> (use local background, # of equal size regions around peaks to use i.e. 2)
		-redundant <#> (Remove redundant sequences matching greater than # percent, i.e. -redundant 0.5)
		-maxN <#> (maximum percentage of N's in sequence to consider for motif finding, default: 0.7)
		-maskMotif <motif file1> [motif file 2]... (motifs to mask before motif finding)
		-opt <motif file1> [motif file 2]... (motifs to optimize or change length of)
		-rand (randomize target and background sequences labels)
		-ref <peak file> (use file for target and background - first argument is list of peak ids for targets)
		-oligo (perform analysis of individual oligo enrichment)
		-dumpFasta (Dump fasta files for target and background sequences for use with other programs)
		-preparse (force new background files to be created)
		-preparsedDir <directory> (location to search for preparsed file and/or place new files)
		-keepFiles (keep temporary files)
		-fdr <#> (Calculate empirical FDR for de novo discovery #=number of randomizations)

	homer2 specific options:
		-homer2 (use homer2 instead of original homer, default)
		-nlen <#> (length of lower-order oligos to normalize in background, default: -nlen 3)
			-nmax <#> (Max normalization iterations, default: 160)
			-neutral (weight sequences to neutral frequencies, i.e. 25%, 6.25%, etc.)
		-olen <#> (lower-order oligo normalization for oligo table, use if -nlen isn't working well)
		-p <#> (Number of processors to use, default: 1)
		-e <#> (Maximum expected motif instance per bp in random sequence, default: 0.01)
		-cache <#> (size in MB for statistics cache, default: 500)
		-quickMask (skip full masking after finding motifs, similar to original homer)
		-minlp <#> (stop looking for motifs when seed logp score gets above #, default: -10)

	Original homer specific options:
		-homer1 (to force the use of the original homer)
		-depth [low|med|high|allnight] (time spent on local optimization default: med)

	"> /dev/null

}


function do_function_MACS_peak_calling () {

	out_dir="MACS2_peaks_minlength25"
	input_dir="picard_summary"

	mkdir -p $outdir/$SAMPLE/$out_dir

	macs2 callpeak -t $outdir/$SAMPLE/$input_dir/${SAMPLE}_bowtie2.sorted.rmDup.sam \
		-f SAM \
		--gsize hs \
		--keep-dup all \
		--outdir $outdir/$SAMPLE/$out_dir \
		--name ${SAMPLE} \
		--verbose 2 \
		--bdg \
		--trackline \
		--min-length 25 \
		--tempdir /lscratch/$SLURM_JOBID

echo -e "
	usage: macs2 [-h] [--version] {callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak} ...

	macs2 -- Model-based Analysis for ChIP-Sequencing

	positional arguments:
	  {callpeak,bdgpeakcall,bdgbroadcall,bdgcmp,bdgopt,cmbreps,bdgdiff,filterdup,predictd,pileup,randsample,refinepeak}
	    callpeak            Main MACS2 Function: Call peaks from alignment results.
	    bdgpeakcall         Call peaks from bedGraph output. Note: All regions on the same chromosome in the bedGraph file should be continuous so only bedGraph files
	                        from MACS2 are accpetable.
	    bdgbroadcall        Call broad peaks from bedGraph output. Note: All regions on the same chromosome in the bedGraph file should be continuous so only bedGraph
	                        files from MACS2 are accpetable.
	    bdgcmp              Deduct noise by comparing two signal tracks in bedGraph. Note: All regions on the same chromosome in the bedGraph file should be continuous so
	                        only bedGraph files from MACS2 are accpetable.
	    bdgopt              Operations on score column of bedGraph file. Note: All regions on the same chromosome in the bedGraph file should be continuous so only
	                        bedGraph files from MACS2 are accpetable.
	    cmbreps             Combine BEDGraphs of scores from replicates. Note: All regions on the same chromosome in the bedGraph file should be continuous so only
	                        bedGraph files from MACS2 are accpetable.
	    bdgdiff             Differential peak detection based on paired four bedgraph files. Note: All regions on the same chromosome in the bedGraph file should be
	                        continuous so only bedGraph files from MACS2 are accpetable.
	    filterdup           Remove duplicate reads at the same position, then save the rest alignments to BED or BEDPE file. If you use '--keep-dup all option', this
	                        script can be utilized to convert any acceptable format into BED or BEDPE format.
	    predictd            Predict d or fragment size from alignment results. In case of PE data, report the average insertion/fragment size from all pairs. *Will NOT
	                        filter duplicates*
	    pileup              Pileup aligned reads with a given extension size (fragment size or d in MACS language). Note there will be no step for duplicate reads
	                        filtering or sequencing depth scaling, so you may need to do certain pre/post-processing.
	    randsample          Randomly sample number/percentage of total reads.
	    refinepeak          (Experimental) Take raw reads alignment, refine peak summits and give scores measuring balance of waston/crick tags. Inspired by SPP.

	optional arguments:
	  -h, --help            show this help message and exit
	  --version             show program's version number and exit

	For command line options of each command, type: macs2 COMMAND -h
	(base) [seifuddinft@cn0630 02_cut_and_tag_analysis]$ 
	(base) [seifuddinft@cn0630 02_cut_and_tag_analysis]$ 
	(base) [seifuddinft@cn0630 02_cut_and_tag_analysis]$ macs2 callpeak -h
	usage: macs2 callpeak [-h] -t TFILE [TFILE ...] [-c [CFILE [CFILE ...]]] [-f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}] [-g GSIZE] [-s TSIZE]
	                      [--keep-dup KEEPDUPLICATES] [--outdir OUTDIR] [-n NAME] [-B] [--verbose VERBOSE] [--trackline] [--SPMR] [--nomodel] [--shift SHIFT]
	                      [--extsize EXTSIZE] [--bw BW] [--d-min D_MIN] [-m MFOLD MFOLD] [--fix-bimodal] [-q QVALUE | -p PVALUE] [--scale-to {large,small}]
	                      [--down-sample] [--seed SEED] [--tempdir TEMPDIR] [--nolambda] [--slocal SMALLLOCAL] [--llocal LARGELOCAL] [--max-gap MAXGAP]
	                      [--min-length MINLEN] [--broad] [--broad-cutoff BROADCUTOFF] [--cutoff-analysis] [--call-summits] [--fe-cutoff FECUTOFF]
	                      [--buffer-size BUFFER_SIZE] [--to-large] [--ratio RATIO]

	optional arguments:
	  -h, --help            show this help message and exit

	Input files arguments:
	  -t TFILE [TFILE ...], --treatment TFILE [TFILE ...]
	                        ChIP-seq treatment file. If multiple files are given as '-t A B C', then they will all be read and pooled together. REQUIRED.
	  -c [CFILE [CFILE ...]], --control [CFILE [CFILE ...]]
	                        Control file. If multiple files are given as '-c A B C', they will be pooled to estimate ChIP-seq background noise.
	  -f {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}, --format {AUTO,BAM,SAM,BED,ELAND,ELANDMULTI,ELANDEXPORT,BOWTIE,BAMPE,BEDPE}
	                        Format of tag file, AUTO, BED or ELAND or ELANDMULTI or ELANDEXPORT or SAM or BAM or BOWTIE or BAMPE or BEDPE. The default
	                        AUTO option will let MACS decide which format (except for BAMPE and BEDPE which should be implicitly set) the file is. Please check the
	                        definition in README. Please note that if the format is set as BAMPE or BEDPE, MACS2 will call its special Paired-end mode to call peaks by
	                        piling up the actual ChIPed fragments defined by both aligned ends, instead of predicting the fragment size first and extending reads. Also
	                        please note that the BEDPE only contains three columns, and is NOT the same BEDPE format used by BEDTOOLS. DEFAULT: AUTO
	  -g GSIZE, --gsize GSIZE
	                        Effective genome size. It can be 1.0e+9 or 1000000000, or shortcuts:'hs' for human (2.7e9), 'mm' for mouse (1.87e9), 'ce' for C. elegans (9e7)
	                        and 'dm' for fruitfly (1.2e8), Default:hs
	  -s TSIZE, --tsize TSIZE
	                        Tag size/read length. This will override the auto detected tag size. DEFAULT: Not set
	  --keep-dup KEEPDUPLICATES
	                        It controls the behavior towards duplicate tags at the exact same location -- the same coordination and the same strand. The 'auto' option
	                        makes MACS calculate the maximum tags at the exact same location based on binomal distribution using 1e-5 as pvalue cutoff; and the 'all'
	                        option keeps every tags. If an integer is given, at most this number of tags will be kept at the same location. Note, if you've used samtools
	                        or picard to flag reads as 'PCR/Optical duplicate' in bit 1024, MACS2 will still read them although the reads may be decided by MACS2 as
	                        duplicate later. If you plan to rely on samtools/picard/any other tool to filter duplicates, please remove those duplicate reads and save a
	                        new alignment file then ask MACS2 to keep all by '--keep-dup all'. The default is to keep one tag at the same location. Default: 1

	Output arguments:
	  --outdir OUTDIR       If specified all output files will be written to that directory. Default: the current working directory
	  -n NAME, --name NAME  Experiment name, which will be used to generate output file names. DEFAULT: NA
	  -B, --bdg             Whether or not to save extended fragment pileup, and local lambda tracks (two files) at every bp into a bedGraph file. DEFAULT: False
	  --verbose VERBOSE     Set verbose level of runtime message. 0: only show critical message, 1: show additional warning message, 2: show process information, 3: show
	                        debug messages. DEFAULT:2
	  --trackline           Tells MACS to include trackline with bedGraph files. To include this trackline while displaying bedGraph at UCSC genome browser, can show name
	                        and description of the file as well. However my suggestion is to convert bedGraph to bigWig, then show the smaller and faster binary bigWig
	                        file at UCSC genome browser, as well as downstream analysis. Require -B to be set. Default: Not include trackline.
	  --SPMR                If True, MACS will SAVE signal per million reads for fragment pileup profiles. It won't interfere with computing pvalue/qvalue during peak
	                        calling, since internally MACS2 keeps using the raw pileup and scaling factors between larger and smaller dataset to calculate statistics
	                        measurements. If you plan to use the signal output in bedGraph to call peaks using bdgcmp and bdgpeakcall, you shouldn't use this option
	                        because you will end up with different results. However, this option is recommended for displaying normalized pileup tracks across many
	                        datasets. Require -B to be set. Default: False

	Shifting model arguments:
	  --nomodel             Whether or not to build the shifting model. If True, MACS will not build model. by default it means shifting size = 100, try to set extsize to
	                        change it. It's highly recommended that while you have many datasets to process and you plan to compare different conditions, aka differential
	                        calling, use both 'nomodel' and 'extsize' to make signal files from different datasets comparable. DEFAULT: False
	  --shift SHIFT         (NOT the legacy --shiftsize option!) The arbitrary shift in bp. Use discretion while setting it other than default value. When NOMODEL is set,
	                        MACS will use this value to move cutting ends (5') towards 5'->3' direction then apply EXTSIZE to extend them to fragments. When this value is
	                        negative, ends will be moved toward 3'->5' direction. Recommended to keep it as default 0 for ChIP-Seq datasets, or -1 * half of EXTSIZE
	                        together with EXTSIZE option for detecting enriched cutting loci such as certain DNAseI-Seq datasets. Note, you can't set values other than 0
	                        if format is BAMPE or BEDPE for paired-end data. DEFAULT: 0.
	  --extsize EXTSIZE     The arbitrary extension size in bp. When nomodel is true, MACS will use this value as fragment size to extend each read towards 3' end, then
	                        pile them up. It's exactly twice the number of obsolete SHIFTSIZE. In previous language, each read is moved 5'->3' direction to middle of
	                        fragment by 1/2 d, then extended to both direction with 1/2 d. This is equivalent to say each read is extended towards 5'->3' into a d size
	                        fragment. DEFAULT: 200. EXTSIZE and SHIFT can be combined when necessary. Check SHIFT option.
	  --bw BW               Band width for picking regions to compute fragment size. This value is only used while building the shifting model. Tweaking this is not
	                        recommended. DEFAULT: 300
	  --d-min D_MIN         Minimum fragment size in basepair. Any predicted fragment size less than this will be excluded. DEFAULT: 20
	  -m MFOLD MFOLD, --mfold MFOLD MFOLD
	                        Select the regions within MFOLD range of high-confidence enrichment ratio against background to build model. Fold-enrichment in regions must
	                        be lower than upper limit, and higher than the lower limit. Use as -m 10 30. This setting is only used while building the shifting model.
	                        Tweaking it is not recommended. DEFAULT:5 50
	  --fix-bimodal         Whether turn on the auto pair model process. If set, when MACS failed to build paired model, it will use the nomodel settings, the --exsize
	                        parameter to extend each tags towards 3' direction. Not to use this automate fixation is a default behavior now. DEFAULT: False

	Peak calling arguments:
	  -q QVALUE, --qvalue QVALUE
	                        Minimum FDR (q-value) cutoff for peak detection. DEFAULT: 0.05. -q, and -p are mutually exclusive.
	  -p PVALUE, --pvalue PVALUE
	                        Pvalue cutoff for peak detection. DEFAULT: not set. -q, and -p are mutually exclusive. If pvalue cutoff is set, qvalue will not be calculated
	                        and reported as -1 in the final .xls file.
	  --scale-to {large,small}
	                        When set to 'small', scale the larger sample up to the smaller sample. When set to 'larger', scale the smaller sample up to the bigger sample.
	                        By default, scale to 'small'. This option replaces the obsolete '--to-large' option. The default behavior is recommended since it will lead to
	                        less significant p/q-values in general but more specific results. Keep in mind that scaling down will influence control/input sample more.
	                        DEFAULT: 'small', the choice is either 'small' or 'large'.
	  --down-sample         When set, random sampling method will scale down the bigger sample. By default, MACS uses linear scaling. Warning: This option will make your
	                        result unstable and irreproducible since each time, random reads would be selected. Consider to use 'randsample' script instead. <not
	                        implmented>If used together with --SPMR, 1 million unique reads will be randomly picked.</not implemented> Caution: due to the implementation,
	                        the final number of selected reads may not be as you expected! DEFAULT: False
	  --seed SEED           Set the random seed while down sampling data. Must be a non-negative integer in order to be effective. DEFAULT: not set
	  --tempdir TEMPDIR     Optional directory to store temp files. DEFAULT: /lscratch/6645448
	  --nolambda            If True, MACS will use fixed background lambda as local lambda for every peak region. Normally, MACS calculates a dynamic local lambda to
	                        reflect the local bias due to the potential chromatin accessibility.
	  --slocal SMALLLOCAL   The small nearby region in basepairs to calculate dynamic lambda. This is used to capture the bias near the peak summit region. Invalid if
	                        there is no control data. If you set this to 0, MACS will skip slocal lambda calculation. *Note* that MACS will always perform a d-size local
	                        lambda calculation while the control data is available. The final local bias would be the maximum of the lambda value from d, slocal, and
	                        llocal size windows. While control is not available, d and slocal lambda won't be considered. DEFAULT: 1000
	  --llocal LARGELOCAL   The large nearby region in basepairs to calculate dynamic lambda. This is used to capture the surround bias. If you set this to 0, MACS will
	                        skip llocal lambda calculation. *Note* that MACS will always perform a d-size local lambda calculation while the control data is available.
	                        The final local bias would be the maximum of the lambda value from d, slocal, and llocal size windows. While control is not available, d and
	                        slocal lambda won't be considered. DEFAULT: 10000.
	  --max-gap MAXGAP      Maximum gap between significant sites to cluster them together. The DEFAULT value is the detected read length/tag size.
	  --min-length MINLEN   Minimum length of a peak. The DEFAULT value is the predicted fragment size d. Note, if you set a value smaller than the fragment size, it may
	                        have NO effect on the result. For BROAD peak calling, try to set a large value such as 500bps. You can also use '--cutoff-analysis' option
	                        with default setting, and check the column 'avelpeak' under different cutoff values to decide a reasonable minlen value.
	  --broad               If set, MACS will try to call broad peaks using the --broad-cutoff setting. Please tweak '--broad-cutoff' setting to control the peak calling
	                        behavior. At the meantime, either -q or -p cutoff will be used to define regions with 'stronger enrichment' inside of broad peaks. The maximum
	                        gap is expanded to 4 * MAXGAP (--max-gap parameter). As a result, MACS will output a 'gappedPeak' and a 'broadPeak' file instead of
	                        'narrowPeak' file. Note, a broad peak will be reported even if there is no 'stronger enrichment' inside. DEFAULT: False
	  --broad-cutoff BROADCUTOFF
	                        Cutoff for broad region. This option is not available unless --broad is set. If -p is set, this is a pvalue cutoff, otherwise, it's a qvalue
	                        cutoff. Please note that in broad peakcalling mode, MACS2 uses this setting to control the overall peak calling behavior, then uses -q or -p
	                        setting to define regions inside broad region as 'stronger' enrichment. DEFAULT: 0.1
	  --cutoff-analysis     While set, MACS2 will analyze number or total length of peaks that can be called by different p-value cutoff then output a summary table to
	                        help user decide a better cutoff. The table will be saved in NAME_cutoff_analysis.txt file. Note, minlen and maxgap may affect the results.
	                        WARNING: May take ~30 folds longer time to finish. The result can be useful for users to decide a reasonable cutoff value. DEFAULT: False

	Post-processing options:
	  --call-summits        If set, MACS will use a more sophisticated signal processing approach to find subpeak summits in each enriched peak region. DEFAULT: False
	  --fe-cutoff FECUTOFF  When set, the value will be used to filter out peaks with low fold-enrichment. Note, MACS2 use 1.0 as pseudocount while calculating fold-
	                        enrichment. DEFAULT: 1.0

	Other options:
	  --buffer-size BUFFER_SIZE
	                        Buffer size for incrementally increasing internal array size to store reads alignment information. In most cases, you don't have to change
	                        this parameter. However, if there are large number of chromosomes/contigs/scaffolds in your alignment, it's recommended to specify a smaller
	                        buffer size in order to decrease memory usage (but it will take longer time to read alignment files). Minimum memory requested for reading an
	                        alignment file is about # of CHROMOSOME * BUFFER_SIZE * 8 Bytes. DEFAULT: 100000

	Obsolete options:
	  --to-large            Obsolete option. Please use '--scale-to large' instead.
	  --ratio RATIO         Obsolete option. Originally designed to normalize treatment and control with customized ratio, now it won't have any effect.
"> /dev/null

}

# do_fastqc no
# do_trimmomatic
# do_trimgalore
# do_fastqc yes
# do_get_fastqc_stats no
# do_get_fastqc_stats yes
# do_bowtie2_build_index
# do_bowtie2_alignment no
# do_bowtie2_alignment_to_spikein_genome no
# do_picard_remove_PCR_duplicates
# do_calculate_fragment_length
# do_file_format_conversion
# do_reproducibility
# do_spike_in_calibration
# do_run_SEACR
# do_HOMER_annotation
do_HOMER_motif_finding	
# do_function_MACS_peak_calling
