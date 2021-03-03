# ChIP-Seq_Cut_and_Tag_Cut_and_Run_scripts
 Scripts to analyze ChIP-Seq including Cut&Tag and Cut&Run data
##### cut-and-tag-analysis.sh
###### Functions: 
###### Please modify any of these functions accordingly
```bash
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
# do_HOMER_motif_finding	
# do_function_MACS_peak_calling
```
***
* do_bowtie2_alignment
	+ current parameters are suitable for cut&tag/cut&run data
	+ adjust parameters for ChIP-Seq data accordingly



