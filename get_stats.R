
library(tidyverse)
library(ggplot2)
library(corrplot)
library(scales)

# 1. Sequencing depth
##=== R command ===## 
## Path to the project and histone list
projPath = "/data/NHLBI_BCB/Tisdale_Lab/08_Bjorg_analysis_FS/01-POGZ_cut_and_tag/02_cut_and_tag_analysis/"
sampleList = c("BGCT1_1_day5", "BGCT1_2_day5", "BGCT1_3_day5", "BGCT1_4_day9", "BGCT1_5_day9", "BGCT1_6_day9")
histList = c("day5", "day9")

## Collect the alignment results from the bowtie2 alignment summary files
alignResult = c()
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  hist = paste(histInfo[1],histInfo[2],sep="_")
  alignRes = read.table(paste0(projPath, "/00-log_files_bowtie2/", hist, ".POGZ.cut.and.tag.err.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(alignRes$V1[9], 1, nchar(as.character(alignRes$V1[9]))-1)
#  histInfo = strsplit(hist, "_")[[1]]
  alignResult = data.frame(Histone = histInfo[3], Replicate = histInfo[2], 
                           SequencingDepth = alignRes$V1[3] %>% as.character %>% as.numeric, 
                           MappedFragNum_hg38 = alignRes$V1[6] %>% as.character %>% as.numeric + alignRes$V1[8] %>% as.character %>% as.numeric, 
                           AlignmentRate_hg38 = alignRate %>% as.numeric)  %>% rbind(alignResult, .)
}
alignResult$Histone = factor(alignResult$Histone, levels = histList)
alignResult %>% mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"))


# 2. Spike-in alignment 
##=== R command ===## 
spikeAlign = c()
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  hist = paste(histInfo[1],histInfo[2],sep="_")
  spikeRes = read.table(paste0(projPath, "/00-log_files_bowtie2_Ecoli/", hist, ".POGZ.cut.and.tag.err.txt"), header = FALSE, fill = TRUE)
  alignRate = substr(spikeRes$V1[9], 1, nchar(as.character(spikeRes$V1[9]))-1)
#  histInfo = strsplit(hist, "_")[[1]]
  spikeAlign = data.frame(Histone = histInfo[3], Replicate = histInfo[2], 
                          SequencingDepth = spikeRes$V1[3] %>% as.character %>% as.numeric, 
                          MappedFragNum_spikeIn = spikeRes$V1[6] %>% as.character %>% as.numeric + spikeRes$V1[8] %>% as.character %>% as.numeric, 
                          AlignmentRate_spikeIn = alignRate %>% as.numeric)  %>% rbind(spikeAlign, .)
}
spikeAlign$Histone = factor(spikeAlign$Histone, levels = histList)
spikeAlign %>% mutate(AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))


# 3. Summarize the alignment to hg38 and E.coli
##=== R command ===## 
alignSummary = left_join(alignResult, spikeAlign, by = c("Histone", "Replicate", "SequencingDepth")) %>%
  mutate(AlignmentRate_hg38 = paste0(AlignmentRate_hg38, "%"), 
         AlignmentRate_spikeIn = paste0(AlignmentRate_spikeIn, "%"))
alignSummary

# 4. Visualizing the sequencing depth and alignment results.
##=== R command ===## 
## Generate sequencing depth boxplot
fig3A = alignResult %>% ggplot(aes(x = Histone, y = SequencingDepth/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Sequencing Depth per Million") +
    xlab("") + 
    ggtitle("A. Sequencing Depth")

fig3B = alignResult %>% ggplot(aes(x = Histone, y = MappedFragNum_hg38/1000000, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Mapped Fragments per Million") +
    xlab("") +
    ggtitle("B. Alignable Fragment (hg38)")

fig3C = alignResult %>% ggplot(aes(x = Histone, y = AlignmentRate_hg38, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("% of Mapped Fragments") +
    xlab("") +
    ggtitle("C. Alignment Rate (hg38)")

fig3D = spikeAlign %>% ggplot(aes(x = Histone, y = AlignmentRate_spikeIn, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Spike-in Alignment Rate") +
    xlab("") +
    ggtitle("D. Alignment Rate (E.coli)")

# ggarrange(fig3A, fig3B, fig3C, fig3D, ncol = 2, nrow=2, common.legend = TRUE, legend="bottom")
print(fig3A)
print(fig3B)
print(fig3C)
print(fig3D)

##=== R command ===## 
## Summarize the duplication information from the picard summary outputs.
dupResult = c()
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  hist = paste(histInfo[1],histInfo[2],sep="_")
  dupRes = read.table(paste0(projPath, hist,"/picard_summary/", hist, "_picard.rmDup.txt"), header = TRUE, fill = TRUE)
# histInfo = strsplit(hist, "_")[[1]]
  dupResult = data.frame(Histone = histInfo[3], Replicate = histInfo[2], MappedFragNum_hg38 = dupRes$READ_PAIRS_EXAMINED[1] %>% as.character %>% as.numeric, DuplicationRate = dupRes$PERCENT_DUPLICATION[1] %>% as.character %>% as.numeric * 100, EstimatedLibrarySize = dupRes$ESTIMATED_LIBRARY_SIZE[1] %>% as.character %>% as.numeric) %>% mutate(UniqueFragNum = MappedFragNum_hg38 * (1-DuplicationRate/100))  %>% rbind(dupResult, .)
}
dupResult$Histone = factor(dupResult$Histone, levels = histList)
alignDupSummary = left_join(alignSummary, dupResult, by = c("Histone", "Replicate", "MappedFragNum_hg38")) %>% mutate(DuplicationRate = paste0(DuplicationRate, "%"))
alignDupSummary

##=== R command ===## 
## generate boxplot figure for the  duplication rate
fig4A = dupResult %>% ggplot(aes(x = Histone, y = DuplicationRate, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Duplication Rate (*100%)") +
    xlab("") 

fig4B = dupResult %>% ggplot(aes(x = Histone, y = EstimatedLibrarySize, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("Estimated Library Size") +
    xlab("") 

fig4C = dupResult %>% ggplot(aes(x = Histone, y = UniqueFragNum, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 18) +
    ylab("# of Unique Fragments") +
    xlab("")

# ggarrange(fig4A, fig4B, fig4C, ncol = 3, common.legend = TRUE, legend="bottom")
print(fig4A)
print(fig4B)
print(fig4C)

##=== R command ===## 
## Collect the fragment size information
fragLen = c()
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  hist = paste(histInfo[1],histInfo[2],sep="_")
#  histInfo = strsplit(hist, "_")[[1]]
  fragLen = read.table(paste0(projPath, hist, "/fragment_length/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[3], Replicate = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
  fragLendup = read.table(paste0(projPath, hist, "/fragment_length_deduplicated/", hist, "_fragmentLen.txt"), header = FALSE) %>% mutate(fragLen = V1 %>% as.numeric, fragCount = V2 %>% as.numeric, Weight = as.numeric(V2)/sum(as.numeric(V2)), Histone = histInfo[3], Replicate = histInfo[2], sampleInfo = hist) %>% rbind(fragLen, .) 
}
# fragLen$sampleInfo = factor(fragLen$sampleInfo, levels = sampleList)
fragLen$Histone = factor(fragLen$Histone, levels = histList)
as_tibble(fragLen)

fragLendup$Histone = factor(fragLendup$Histone, levels = histList)
as_tibble(fragLendup)

## Generate the fragment size density plot (violin plot)
fig5A = fragLen %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
    geom_violin(bw = 5) +
    scale_y_continuous(breaks = seq(0, 800, 50)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ggpubr::rotate_x_text(angle = 20) +
    ylab("Fragment Length") +
    xlab("")

fig5B = fragLen %>% ggplot(aes(x = fragLen, y = fragCount, group = sampleInfo, linetype = Replicate, color = Histone)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

## Generate the fragment size density plot (violin plot) - deduplicated BAM file
fig6A = fragLendup %>% ggplot(aes(x = sampleInfo, y = fragLen, weight = Weight, fill = Histone)) +
    geom_violin(bw = 5) +
    scale_y_continuous(breaks = seq(0, 800, 50)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ggpubr::rotate_x_text(angle = 20) +
    ylab("Fragment Length") +
    xlab("")

fig6B = fragLendup %>% ggplot(aes(x = fragLen, y = fragCount, group = sampleInfo, linetype = Replicate, color = Histone)) +
  geom_line(size = 1) +
  scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma") +
  theme_bw(base_size = 20) +
  xlab("Fragment Length") +
  ylab("Count") +
  coord_cartesian(xlim = c(0, 500))

# ggarrange(fig5A, fig5B, ncol = 2)
print(fig5A)
print(fig5B)
print(fig6A)
print(fig6B)

##== R command ==##
reprod = c()
fragCount = NULL
for(hist in sampleList){
  
  histInfo = strsplit(hist, "_")[[1]]
  hist = paste(histInfo[1],histInfo[2],sep="_")
  if(is.null(fragCount)){
    
    fragCount = read.table(paste0(projPath, hist, "/bowtie2/", hist, "_bowtie2.mapped.clean.fragmentsCount.bin500.bed"), header = FALSE) 
    colnames(fragCount) = c("chrom", "bin", hist)
  
  }else{
    
    fragCountTmp = read.table(paste0(projPath, hist, "/bowtie2/", hist, "_bowtie2.mapped.clean.fragmentsCount.bin500.bed"), header = FALSE)
    colnames(fragCountTmp) = c("chrom", "bin", hist)
    fragCount = full_join(fragCount, fragCountTmp, by = c("chrom", "bin"))
    
  }
}

M = cor(fragCount %>% select(-c("chrom", "bin")) %>% log2(), use = "complete.obs") 

corrplot(M, method = "color", outline = T, addgrid.col = "darkgray", order="hclust", cl.pos = "b", tl.col = "indianred4", tl.cex = 1, cl.cex = 1, addCoef.col = "black", number.digits = 2, number.cex = 1, col = colorRampPalette(c("midnightblue","white","darkred"))(100))

##=== R command ===## 
scaleFactor = c()
multiplier = 10000
for(hist in sampleList){

  histInfo = strsplit(hist, "_")[[1]]
  hist = paste(histInfo[1],histInfo[2],sep="_")

  spikeDepth = read.table(paste0(projPath, hist, "/bowtie2_Ecoli/", hist, "_bowtie2_spikeIn.seqDepth"), header = FALSE, fill = TRUE)$V1[1]
  
#  histInfo = strsplit(hist, "_")[[1]]
  scaleFactor = data.frame(scaleFactor = multiplier/spikeDepth, Histone = histInfo[3], Replicate = histInfo[2])  %>% rbind(scaleFactor, .)
}
# scaleFactor$Histone = factor(scaleFactor$Histone, levels = histList)
left_join(alignDupSummary, scaleFactor, by = c("Histone", "Replicate"))

##=== R command ===##
## Generate sequencing depth boxplot
fig7A = scaleFactor %>% ggplot(aes(x = Histone, y = scaleFactor, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ylab("Spike-in Scalling Factor") +
    xlab("")

normDepth = inner_join(scaleFactor, alignResult, by = c("Histone", "Replicate")) %>% mutate(normDepth = MappedFragNum_hg38 * scaleFactor)

fig7B = normDepth %>% ggplot(aes(x = Histone, y = normDepth, fill = Histone)) +
    geom_boxplot() +
    geom_jitter(aes(color = Replicate), position = position_jitter(0.15)) +
    scale_fill_viridis(discrete = TRUE, begin = 0.1, end = 0.9, option = "magma", alpha = 0.8) +
    scale_color_viridis(discrete = TRUE, begin = 0.1, end = 0.9) +
    theme_bw(base_size = 20) +
    ylab("Normalization Fragment Count") +
    xlab("") + 
    coord_cartesian(ylim = c(1000000, 130000000))

#ggarrange(fig6A, fig6B, ncol = 2, common.legend = TRUE, legend="bottom")
print(fig7A)
print(fig7B)

##=== R command ===## 
peakN = c()
peakWidth = c()
peakType = "top0.01"
for(hist in sampleList){
  histInfo = strsplit(hist, "_")[[1]]
  hist = paste(histInfo[1],histInfo[2],sep="_")
  if(histInfo[1] != "IgG"){
    for(type in peakType){
      peakInfo = read.table(paste0(projPath, hist, "/SEACR_peaks/", hist,".stringent.bed"), header = FALSE, fill = TRUE)  %>% mutate(width = abs(V3-V2))
      peakN = data.frame(peakN = nrow(peakInfo), peakType = type, Histone = histInfo[3], Replicate = histInfo[2]) %>% rbind(peakN, .)
      peakWidth = data.frame(width = peakInfo$width, peakType = type, Histone = histInfo[3], Replicate = histInfo[2])  %>% rbind(peakWidth, .)
    }
  }
}
peakN %>% select(Histone, Replicate, peakType, peakN)



dev.off()
