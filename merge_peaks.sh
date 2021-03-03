#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

module load bedops

inputdir="/data/NHLBI_BCB/Tisdale_Lab/08_Bjorg_analysis_FS/01-POGZ_cut_and_tag/02_cut_and_tag_analysis"
outputdir="/data/NHLBI_BCB/Tisdale_Lab/08_Bjorg_analysis_FS/01-POGZ_cut_and_tag/02_cut_and_tag_analysis/04-peaks_overlap_between_donors"

SAMPLE1="BGCT1_1"
SAMPLE2="BGCT1_2"
SAMPLE3="BGCT1_3"
# BGCT1_1
# BGCT1_2
# BGCT1_3	

cat $inputdir/$SAMPLE1/MACS2_peaks_minlength25_HOMER/${SAMPLE1}_peaks.narrowPeak.homer.input.bed $inputdir/$SAMPLE2/MACS2_peaks_minlength25_HOMER/${SAMPLE2}_peaks.narrowPeak.homer.input.bed $inputdir/$SAMPLE3/MACS2_peaks_minlength25_HOMER/${SAMPLE3}_peaks.narrowPeak.homer.input.bed | awk '{print $1"\t"$2"\t"$3"\t"$4"\t"$7}' | sort -k1,1 -k2,2n - > "$outputdir/day5.sorted.mapfile.bed"

# while read SAMPLE READ1 READ2; do
#	cat $inputdir/$SAMPLE/${SAMPLE}_peaks.narrowPeak.homer.input.bed
# done < "/data/NHLBI_BCB/Tisdale_Lab/08_Bjorg_analysis_FS/01-POGZ_cut_and_tag/POGZ_cut_and_tag_sampleIDs_read1_read2.txt"

bedops -u $inputdir/$SAMPLE1/MACS2_peaks_minlength25_HOMER/${SAMPLE1}_peaks.narrowPeak.homer.input.bed $inputdir/$SAMPLE2/MACS2_peaks_minlength25_HOMER/${SAMPLE2}_peaks.narrowPeak.homer.input.bed $inputdir/$SAMPLE3/MACS2_peaks_minlength25_HOMER/${SAMPLE3}_peaks.narrowPeak.homer.input.bed \
    | bedmap --echo --echo-map-score --echo-map-id-uniq --echo-overlap-size --echo-map-size --count --echo-ref-name --delim '\t' --mean --bp-ovr 10 - "$outputdir/day5.sorted.mapfile.bed" \
    | awk -F'\t' '(split($13, a, ";") > 1)' \
    > day5.bed

awk 'NR==FNR{A[$4]=$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13"\t"$14"\t"$15"\t"$16"\t"$17"\t"$18;next}{print$0"\t"A[$1]}' day5.bed day5.homer.annotated.bed > day5.merged.txt

echo -e "

# chr1  629136  629846
# 629376 - 629421
# D4: chr1  629136  629846  BGCT1_4_peak_6  172 . 2.14877 22.296  17.2108 196 710
# D5: chr1  629365  629421  BGCT1_5_peak_3  37  . 2.44539 8.8205  3.7247  23  56
# BGCT1_6_peak_4
# D6: chr1  629184  629989  BGCT1_6_peak_4  440 . 2.14259 49.2523 44.0835 230 805

bedmap
  citation: http://bioinformatics.oxfordjournals.org/content/28/14/1919.abstract
            https://doi.org/10.1093/bioinformatics/bts277
  version:  2.4.38 (typical)
  authors:  Shane Neph & Scott Kuehn
                                                                                                    
 USAGE: bedmap [process-flags] [overlap-option] <operation(s)...> <ref-file> [map-file]             
     Any input file must be sorted per the sort-bed utility.                                        
     The program accepts BED and Starch file formats.                                               
     You may use '-' for a BED file to indicate the input comes from stdin.                         
                                                                                                    
     Traverse <ref-file>, while applying <operation(s)> on qualified, overlapping elements from     
       <map-file>.  Output is one line for each line in <ref-file>, sent to standard output.  There 
       is no limit on the number of operations you can specify to compute in one bedmap call.       
     If <map-file> is omitted, the given file is treated as both the <ref-file> and <map-file>.     
       This usage is more efficient than specifying the same file twice.                            
     Arguments may be given in any order before the input file(s).                                  
                                                                                                    
    Process Flags:                                                                                  
     --------                                                                                       
      --chrom <chromosome>  Jump to and process data for given <chromosome> only.                   
      --delim <delim>       Change output delimiter from '|' to <delim> between columns (e.g. '\t').
      --ec                  Error check all input files (slower).                                   
      --faster              (advanced) Strong input assumptions are made.  Compatible with:         
                              --bp-ovr, --range, --fraction-both, and --exact overlap options only. 
      --header              Accept headers (VCF, GFF, SAM, BED, WIG) in any input file.             
      --help                Print this message and exit successfully.                               
      --min-memory          Minimize memory usage (slower).                                         
      --multidelim <delim>  Change delimiter of multi-value output columns from ';' to <delim>.     
      --prec <int>          Change the post-decimal precision of scores to <int>.  0 <= <int>.      
      --sci                 Use scientific notation for score outputs.                              
      --skip-unmapped       Print no output for a row with no mapped elements.                      
      --sweep-all           Ensure <map-file> is read completely (helps to prevent broken pipes).   
      --unmapped-val <val>  Print <val> on unmapped --echo-map* and --min/max-element* operations.  
                              The default is to print nothing.                                      
      --version             Print program information.                                              
                                                                                                    
                                                                                                    
    Overlap Options (At most, one may be selected.  By default, --bp-ovr 1 is used):                
     --------                                                                                       
      --bp-ovr <int>           Require <int> bp overlap between elements of input files.            
      --exact                  First 3 fields from <map-file> must be identical to <ref-file>'s.    
      --fraction-both <val>    Both --fraction-ref <val> and --fraction-map <val> must be true to   
                                 qualify as overlapping.  Expect 0 < val <= 1.                      
      --fraction-either <val>  Either --fraction-ref <val> or --fraction-map <val> must be true to  
                                 qualify as overlapping.  Expect 0 < val <= 1.                      
      --fraction-map <val>     The fraction of the element's size from <map-file> that must overlap 
                                 the element in <ref-file>.  Expect 0 < val <= 1.                   
      --fraction-ref <val>     The fraction of the element's size from <ref-file> that must overlap 
                                 an element in <map-file>.  Expect 0 < val <= 1.                    
      --range <int>            Grab <map-file> elements within <int> bp of <ref-file>'s element,    
                                 where 0 <= int.  --range 0 is an alias for --bp-ovr 1.             
                                                                                                    
                                                                                                    
    Operations:  (Any number of operations may be used any number of times.)                        
     ----------                                                                                     
      SCORE:                                                                                        
       <ref-file> must have at least 3 columns and <map-file> 5 columns.                            
                                                                                                    
      --cv                The result of --stdev divided by the result of --mean.
      --kth <val>         Generalized median. Report the value, x, such that the fraction <val>
                            of overlapping elements' scores from <map-file> is less than x,
                            and the fraction 1-<val> of scores is greater than x.  0 < val <= 1.
      --mad <mult=1>      The median absolute deviation of overlapping elements in <map-file>.
                            Multiply mad score by <mult>.  0 < mult, and mult is 1 by default.
      --max               The highest score from overlapping elements in <map-file>.
      --max-element       A (non-random) highest-scoring and overlapping element in <map-file>.
      --max-element-rand  A random highest-scoring and overlapping element in <map-file>.
      --mean              The average score from overlapping elements in <map-file>.
      --median            The median score from overlapping elements in <map-file>.
      --min               The lowest score from overlapping elements in <map-file>.
      --min-element       A (non-random) lowest-scoring and overlapping element in <map-file>.
      --min-element-rand  A random lowest-scoring and overlapping element in <map-file>.
      --stdev             The square root of the result of --variance.
      --sum               Accumulated scores from overlapping elements in <map-file>.
      --tmean <low> <hi>  The mean score from overlapping elements in <map-file>, after
                            ignoring the bottom <low> and top <hi> fractions of those scores.
                            0 <= low <= 1.  0 <= hi <= 1.  low+hi <= 1.
      --variance          The variance of scores from overlapping elements in <map-file>.
      --wmean             Weighted mean, scaled in proportion to the coverage of the <ref-file>
                            element by each overlapping <map-file> element.
     
     ----------
      NON-SCORE:
       <ref-file> must have at least 3 columns.
       For --echo-map-id/echo-map-id-uniq, <map-file> must have at least 4 columns.
       For --echo-map-score, <map-file> must have at least 5 columns.
       For all others, <map-file> requires at least 3 columns.

      --bases             The total number of overlapping bases from <map-file>.
      --bases-uniq        The number of distinct bases from <ref-file>'s element covered by
                            overlapping elements in <map-file>.
      --bases-uniq-f      The fraction of distinct bases from <ref-file>'s element covered by
                            overlapping elements in <map-file>.
      --count             The number of overlapping elements in <map-file>.
      --echo              Print each line from <ref-file>.
      --echo-map          List all overlapping elements from <map-file>.
      --echo-map-id       List IDs from all overlapping <map-file> elements.
      --echo-map-id-uniq  List unique IDs from overlapping <map-file> elements.
      --echo-map-range    Print genomic range of overlapping elements from <map-file>.
      --echo-map-score    List scores from overlapping <map-file> elements.
      --echo-map-size     List the full length of every overlapping element.
      --echo-overlap-size List lengths of overlaps.
      --echo-ref-name     Print the first 3 fields of <ref-file> using chrom:start-end format.
      --echo-ref-row-id   Print 'id-' followed by the line number of <ref-file>.
      --echo-ref-size     Print the length of each line from <ref-file>.
      --indicator         Print 1 if there exists an overlapping element in <map-file>, 0 otherwise.

" > /dev/null

echo -e "

bedops
  citation: http://bioinformatics.oxfordjournals.org/content/28/14/1919.abstract
            https://doi.org/10.1093/bioinformatics/bts277
  version:  2.4.38 (typical)
  authors:  Shane Neph & Scott Kuehn

      USAGE: bedops [process-flags] <operation> <File(s)>*

          Every input file must be sorted per the sort-bed utility.
          Each operation requires a minimum number of files as shown below.
            There is no fixed maximum number of files that may be used.
          Input files must have at least the first 3 columns of the BED specification.
          The program accepts BED and Starch file formats.
          May use '-' for a file to indicate reading from standard input (BED format only).

      Process Flags:
          --chrom <chromosome> Jump to and process data for given <chromosome> only.
          --ec                 Error check input files (slower).
          --header             Accept headers (VCF, GFF, SAM, BED, WIG) in any input file.
          --help               Print this message and exit successfully.
          --help-<operation>   Detailed help on <operation>.
                                 An example is --help-c or --help-complement
          --range L:R          Add 'L' bp to all start coordinates and 'R' bp to end
                                 coordinates. Either value may be + or - to grow or
                                 shrink regions.  With the -e/-n operations, the first
                                 (reference) file is not padded, unlike all other files.
          --range S            Pad or shrink input file(s) coordinates symmetrically by S.
                                 This is shorthand for: --range -S:S.
          --version            Print program information.

      Operations: (choose one of)
          -c, --complement [-L] File1 [File]*
          -d, --difference ReferenceFile File2 [File]*
          -e, --element-of [bp | percentage] ReferenceFile File2 [File]*
                 by default, -e 100% is used.  'bedops -e 1' is also popular.
          -i, --intersect File1 File2 [File]*
          -m, --merge File1 [File]*
          -n, --not-element-of [bp | percentage] ReferenceFile File2 [File]*
                 by default, -n 100% is used.  'bedops -n 1' is also popular.
          -p, --partition File1 [File]*
          -s, --symmdiff File1 File2 [File]*
          -u, --everything File1 [File]*
          -w, --chop [bp] [--stagger <nt>] [-x] File1 [File]*
                 by default, -w 1 is used with no staggering.
      
Example: bedops --range 10 -u file1.bed
      NOTE: Only operations -e|n|u preserve all columns (no flattening)

" > /dev/null