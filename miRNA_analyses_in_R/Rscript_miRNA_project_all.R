rm(list = ls())
library(openxlsx)

draw = read.xlsx('directory/read counts.xlsx')
dquality = read.xlsx('directory/quality unique read counts.xlsx')
druns = read.csv('directory/information on runs for unique samples.csv')

## Barplot of all read counts (matched and unmatched) ordered by matched read count
quality_plot = dquality[order(dquality$Construct_Matched_Without_Mito_and_Y),] # order according to matched reads
quality_plot = subset(quality_plot, Condition != "spike_in_plasmid_DNA_from_CP1608") # exclude library
quality_plot = cbind(quality_plot, Number_of_Screens = c(1:nrow(quality_plot)),
                     Unmatched = (quality_plot$Sample_Matched - quality_plot$Construct_Matched_Without_Mito_and_Y))
pdf('Matched_and_Unmatched_Read_Counts_of_all_Unique_Replicates.pdf', height = 10, width = 15)
par(mar = c(5.1, 4.6, 4.1, 2.1))
barplot(rbind(quality_plot$Construct_Matched_Without_Mito_and_Y, quality_plot$Unmatched), width = 1,
        space = 0, names.arg = quality_plot$Number_of_Screens, beside = FALSE, col =
          rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549))), cex.main = 2,
        border = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549))), xlab = 'Number of Replicates', ylab = 'Number of Reads', cex.axis = 1.5, 
        main = 'Matched and Unmatched Reads', cex.lab = 1.5, cex.names = 1.5)
legend(-3,12000000, legend = c('Reads matched to CP1608', 'Total number of reads'),
       fill = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549))), box.lty = 0, cex = 1.5 )
abline(h = 1000000, col = rgb(0.8, 0, 0), lwd = 3) # After 3rd screen: 3 samples are below, only one with high enough reads to be visible
dev.off()

## Same plot for screens not meeting 500x coverage
coverage_plot = subset(quality_plot, Construct_Matched_Without_Mito_and_Y < (500*9434))
pdf('Matched_and_Unmatched_Read_Counts_of_Replicates_Below_500x_Coverage.pdf', width = 12, height = 10)
par(mar = c(7.1, 6.6, 4.1, 2.1))
barplot(rbind(coverage_plot$Construct_Matched_Without_Mito_and_Y, coverage_plot$Unmatched), width = 1, space = 0,
        names.arg = sub('_REP_', ' ', coverage_plot$Condition), las = 2, beside = FALSE, log = 'y',
        col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549))),
        border = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549))),
        cex.axis = 1.5, cex.main = 2, cex.lab = 1.5,
        main = 'Matched and Unmatched Reads of Replicates below 500x Coverage')
title(ylab = 'Number of Reads', mgp = c(5,1,0), cex.lab = 1.5)
legend(26.2,4700, legend = c('Reads matched to CP1608', 'Total number of reads'), cex = 1.5,
       fill = rgb(rbind(c(0.28, 0.52, 0.72), c(0.7291, 0.6588, 0.4549))), box.lty = 0, bg = 'white')
abline(h = 1000000, col = rgb(0.8, 0, 0), lwd = 3)
dev.off()


## Remove failed conditions
removed = rbind.data.frame(cbind(draw$Construct_Barcode,raw$Construct_IDs),
                           c('Construct_Matched_With_Mito','Construct_Matched_With_Mito'),
                           c('Sample_Matched','Sample_Matched'),
                           c('Fraction_Matched','Fraction_Matched'),
                           c('Mito_and_Y_counts', 'Mito_and_Y_counts'),
                           c('Construct_Matched_Without_Mito_and_Y', 'Construct_Matched_Without_Mito_and_Y'),
                           c('Fraction_Matched_Without_Mito_and_Y', 'Fraction_Matched_Without_Mito_and_Y'),
                           c('Sequencing_Run', 'Sequencing_Run'),
                           c('Plate', 'Plate'),
                           stringsAsFactors = FALSE)
colnames(removed) = c('Construct_Barcode','Construct_ID')

set_coverage = 1000000 ## Set cutoff for coverage, all read counts below will be deleted
set_match = 0.6 ## Set cutoff for fraction of read counts matched to library

cquality = subset(dquality, Fraction_Matched >= set_match & Construct_Matched_Without_Mito_and_Y >= set_coverage)
fails_quality = subset(dquality, Fraction_Matched < set_match | Construct_Matched_Without_Mito_and_Y < set_coverage)
craw = draw
cruns = druns
for (ii in 1:nrow(fails_quality)){
  temp_pos = which(colnames(craw) == fails_quality$Condition[ii]) # number of col in which the current sample is located (i.e. which col has to be removed)
  removed = cbind(removed,
                  rbind(as.data.frame(as.numeric(craw[,temp_pos])),
                        fails_quality$Construct_Matched_With_Mito[ii], fails_quality$Sample_Matched[ii],
                        fails_quality$Fraction_Matched[ii], fails_quality$Mito_and_Y_counts[ii],
                        fails_quality$Construct_Matched_Without_Mito_and_Y[ii],
                        fails_quality$Fraction_Matched_Without_Mito_and_Y[ii],
                        cruns[1,(temp_pos-2)],cruns[2,(temp_pos-2)]))
  colnames(removed)[length(colnames(removed))] = fails_quality$Condition[ii]
  craw = craw[,-temp_pos]
  cruns = cruns[,-(temp_pos-2)] # same with information on runs (shifted by 2)
}

## Do any samples fail the match cutoff while meeting the coverage?
as.logical(nrow(subset(dquality, Fraction_Matched < set_match & Construct_Matched_Without_Mito_and_Y >= set_coverage)))

## List with cell lines that fail 500x coverage with more than 1 replicate
multiple_below = data.frame( # remove replicate name and include all conditions that have multiple entries
  Cell_Line = sub('_REP_.', '', coverage_plot$Condition[which(duplicated(sub('_REP_.', '', coverage_plot$Condition)))]),
  Number_of_Replicates = rep(2, times = length(coverage_plot$Condition[which(duplicated(sub(
    '_REP_.', '', coverage_plot$Condition)))]))) 
while (anyDuplicated(multiple_below$Cell_Line)) { # redo that until there are no duplicates present 
  for (ii in 1:length(anyDuplicated(multiple_below$Cell_Line))) { # ATTENTION: Potential to be endless loop
    multiple_below$Number_of_Replicates[ # be cautious!
      which(multiple_below$Cell_Line == multiple_below$Cell_Line[anyDuplicated(multiple_below$Cell_Line)[ii]])[1]] =
      multiple_below$Number_of_Replicates[
        which(multiple_below$Cell_Line == multiple_below$Cell_Line[anyDuplicated(multiple_below$Cell_Line)[ii]])[1]] + 1
    multiple_below = multiple_below[
      -(which(multiple_below$Cell_Line == multiple_below$Cell_Line[anyDuplicated(multiple_below$Cell_Line)[ii]])[2]),]
  }
  
}

## Sanity Checks
SanityCheck2_1 = (ncol(draw) == (ncol(craw) +  ncol(removed) - 2))
SanityCheck2_2 = all(cquality$Condition == colnames(craw[,-(1:2)]))
SanityCheck2_3 = any(removed[nrow(removed),3:ncol(removed)] >= set_match)
SanityCheck2_4 = all(apply(craw[,3:ncol(craw)], 2, sum) >= set_match)
SanityCheck2_5 = (nrow(dquality) == (nrow(cquality) + ncol(removed) - 2))
SanityCheck2_6 = all(removed[nrow(removed)-3,3:ncol(removed)] ==
                       apply(type.convert(removed[1:(nrow(removed)-8),3:ncol(removed)], as.is = FALSE), 2, sum))
SanityCheck2_7 = all(colnames(cruns) == colnames(craw)[-(1:2)])
# SanityCheck2_1 # Do the divided tables together have the same number of conditions as the original table? (read counts)
# SanityCheck2_2 # Are the names of the conditions meeting the cutoffs the same for the quality and read count table?
# SanityCheck2_3 # Do the fractions of all excluded conditions fail to meet the cutoff?
# SanityCheck2_4 # Do the matched read counts of all included conditions meet the cutoff?
# SanityCheck2_5 # Do the divided tables together have the same number of conditions as the original table? (quality)
# SanityCheck2_6 # Do the the read counts of the removed samples add up to the read counts extracted from the quality table?
# SanityCheck2_7 # Are all the names of the conditions equal between the file with information on sequencing run and the read counts? 
all(SanityCheck2_1,SanityCheck2_2,SanityCheck2_3,SanityCheck2_4,SanityCheck2_5,SanityCheck2_6,SanityCheck2_7)

write.xlsx(craw, 'read counts above thresholds.xlsx')
write.xlsx(removed, 'unique samples beneath threshold.xlsx')
write.xlsx(cquality, 'quality above threshold.xlsx')
write.csv2(cruns, 'information on runs for conditions above threshold.csv')
write.table(multiple_below, 'Cell lines with multiple replicates below 500 x coverage.txt', row.names = FALSE)


reads = read.xlsx('directory/read counts above thresholds.xlsx')
quality_complete = read.xlsx('directory/quality unique read counts.xlsx')
quality =  read.xlsx('directory/quality above threshold.xlsx')
## Descriptives for conditions (cell line replicates)
cell_lines = data.frame(Condition = colnames(reads[,-(1:2)]), Mean = colMeans(reads[,-(1:2)]),
                        Median = apply(reads[,-(1:2)], 2, median), Max_Reads = apply(reads[,-(1:2)], 2, max),
                        Min_Reads = apply(reads[,-(1:2)], 2, min), Missing_guides = apply(reads[,-(1:2)]==0, 2, sum),
                        Reads_Without_Mito_and_Y = quality$Construct_Matched_Without_Mito_and_Y)

## Descriptives for guides
reads_guide = reads[,-(which(colnames(reads) == "spike_in_plasmid_DNA_from_CP1608"))] ## Remove library
guide_des = data.frame(reads_guide[,1:2], Mean = rowMeans((reads_guide[,-(1:2)])),
                       Median = apply(reads_guide[,-(1:2)], 1, median),
                       Max_Reads = apply(reads_guide[,-(1:2)], 1, max), Min_Reads = apply(reads_guide[,-(1:2)], 1, min),
                       Missing_in_lines = apply(reads_guide[,-(1:2)]==0, 1, sum))
SanityCheck1 = all(quality$Condition == colnames(reads[,-(1:2)]))
SanityCheck1

## build the different intervals
sub_pie_0 = subset(cell_lines, Missing_guides == 0) 
sub_pie_1 = subset(cell_lines, Missing_guides > 0 & Missing_guides <= 10)
sub_pie_2 = subset(cell_lines, Missing_guides > 10 & Missing_guides <= 20)
sub_pie_3 = subset(cell_lines, Missing_guides > 20 & Missing_guides <= 50)
sub_pie_4 = subset(cell_lines, Missing_guides > 50)

## bind information on intervals and number of replicates in the corresponding interval
inter_pie = data.frame(Intervals = c('0', '1-10', '11-20', '21-50', '> 50'),
                       Number_of_Conditions = c(nrow(sub_pie_0), nrow(sub_pie_1), nrow(sub_pie_2), nrow(sub_pie_3),
                                                nrow(sub_pie_4)))

pdf('Fractions_of_Replicates_with_Missing_sgRNAs.pdf')
pie(inter_pie$Number_of_Conditions, inter_pie$Intervals, main ='Fractions of Replicates with Missing sgRNAs',
    col = rgb(rbind(c(0.9411765, 0.9725490, 1), c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549),
                    c(0.9333333, 0.3607843, 0.2588235), c(0.8, 0, 0)), alpha = c(1, 1, 0.7, 1, 1)))
dev.off()

ordered = cell_lines[order(cell_lines$Missing_guides),]
plot(ordered$Missing_guides, pch = 20, col = rgb(0.2117647, 0.3921569, 0.5450980), bty = 'n', xlab = 'Replicates',
     ylab = 'Number of missing guides')
hist(cell_lines$Missing_guides, breaks = 100, col = rgb(0.2117647, 0.3921569, 0.5450980), bty = 'n',
     xlab = 'Number of missing guides', main = '')
write.xlsx(cell_lines, 'Descriptives for cell lines.xlsx')
write.xlsx(guide_des, 'Descriptives for guides.xlsx')

## test the linear correlation of number of matched read counts with number of missing guides in the corresponding
## replicates
cor.test(cell_lines$Missing_guides, cell_lines$Reads_Without_Mito_and_Y, 'two.sided', 'spearman')
pdf('Correlation_Read_Counts_Missing_Guides.pdf')
plot(cell_lines$Missing_guides, cell_lines$Reads_Without_Mito_and_Y, xlab = 'Number of Missing Guides',
     ylab = 'Read Counts', pch = 20)
abline(lm(Reads_Without_Mito_and_Y ~ Missing_guides, data = cell_lines), col = rgb(0.8, 0, 0))
dev.off()

library(CRISPRcleanR)

##load annotation file for lentiG-miR (here called CP1608)

## Read csv file annotation with all guides (n=9435)
CP1608.anno = read.csv(file="CP1608_anno.csv", row.names = 1) # first column will be rownames
## Load the clean annotation file without NO_SITE guides (n=8931) in the same way
CP1608.anno.clean = read.csv(file="CP1608_anno_clean.csv", row.names = 1)
## Load essential genes to be disregarded in segmentation
control_essentials = read.csv(file="Control_essentials.csv", header = FALSE)
essential_string = as.vector(control_essentials[,1])
##load read counts file with all unique replicates good enough for analysis
counts = read.xlsx('directory/read counts above thresholds.xlsx')

## Data frame with sgRNA, gene and read counts of the lentiG-miR plasmid to add last column with read counts of a replicate in a loop
## Data frame needed for normalisation and calculation of log fold changes
single_counts =  data.frame(sgRNA = rownames(CP1608.anno), gene = CP1608.anno$GENES, CP1608_plasmid = counts$spike_in_plasmid_DNA_from_CP1608,
                            Replicate = rep(0, nrow(counts)))
## Generate shortened data frame with read counts to run a loop for all replicates
counts_short = counts[,-(1:2)] # remove information on constructs
## Prepare data frame to take up normalised counts of each condition
counts_norm = data.frame(sgRNA = rownames(CP1608.anno), gene = CP1608.anno$GENES)
## Prepare data frame to take up LFCs of each condition
lfc = data.frame(sgRNA = rownames(CP1608.anno), gene = CP1608.anno$GENES)
## Prepare the same data frames for corrected counts and LFCs
counts_cor = data.frame(sgRNA = rownames(CP1608.anno), gene = CP1608.anno$GENES)
lfc_cor = data.frame(sgRNA = rownames(CP1608.anno.clean), gene = CP1608.anno.clean$GENES)

for (ii in 1:ncol(counts_short)) {
  single_counts$Replicate = counts_short[,ii] # exchange the column for the sample to the current replicate
  normANDfcs = ccr.NormfoldChanges(Dframe = single_counts, # calculate normalised read counts and log fold changes for
                                   libraryAnnotation = CP1608.anno, # current replicate
                                   saveToFig = TRUE,
                                   EXPname = colnames(counts_short)[ii],
                                   outdir = 'directory/normANDfcs/')
  counts_norm = cbind(counts_norm, normANDfcs$norm_counts$Replicate) # take the generated normalised counts of the condition and add it to the data frame
  colnames(counts_norm)[ii + 2] = colnames(counts_short)[ii] # name the added column like the replicate
  lfc = cbind(lfc, normANDfcs$logFCs$Replicate) # same for log fold change
  colnames(lfc)[ii + 2] = colnames(counts_short)[ii] # same for log fold change
  gwSortedFCs = ccr.logFCs2chromPos(normANDfcs$logFCs,CP1608.anno.clean) # sort sgRNA by genome coordinates
  correctedFCs = ccr.GWclean(gwSortedFCs, # correct LFCs
                             display = TRUE,
                             label = colnames(counts_short)[ii],
                             ignoredGenes = essential_string,
                             saveTO = 'directory/correctedFCs')
  lfc_cor = cbind(lfc_cor, correctedFCs$corrected_logFCs$correctedFC) # add corrected LFCs for condition
  colnames(lfc_cor)[ii + 2] = colnames(counts_short)[ii]
  correctedCounts = ccr.correctCounts(colnames(counts_short)[ii], # corrected read counts from corrected LFCs 
                                      normANDfcs$norm_counts, correctedFCs, CP1608.anno.clean, # leaves other guides uncorrected but still in the data
                                      OutDir = 'directory/correctedCounts')
  counts_cor = cbind(counts_cor, correctedCounts$Replicate) # add corrected read counts for condition
  colnames(counts_cor)[ii + 2] = colnames(counts_short)[ii]
}
lfc_cor$gene = correctedFCs$corrected_logFCs$genes # Change gene names because order changed according to position in genome
lfc_cor$sgRNA = rownames(correctedFCs$corrected_logFCs) # Change guide names because order changed according to position in genome
lfc_library = lfc[,(which(colnames(lfc) == 'spike_in_plasmid_DNA_from_CP1608'))] # save LFCs of library for sanity check
lfc = lfc[,-(which(colnames(lfc) == 'spike_in_plasmid_DNA_from_CP1608'))] # delete log fold changes for the library
lfc_cor_library = lfc_cor[,(which(colnames(lfc_cor) == 'spike_in_plasmid_DNA_from_CP1608'))] # see above
lfc_cor = lfc_cor[,-(which(colnames(lfc_cor) == 'spike_in_plasmid_DNA_from_CP1608'))] # see above
nosites = subset(lfc, substr(gene, 1, 8) == 'NO_SITE_') 

## Sanity Checks for normalisation and log fold change
SanityCheck1_1 = ncol(counts) == ncol(counts_norm)
SanityCheck1_2 = ncol(counts) == ncol(lfc) + 1
SanityCheck1_3 = nrow(counts) == nrow(counts_norm)
SanityCheck1_4 = nrow(counts) == nrow(lfc)
SanityCheck1_5 = all(colnames(counts)[-(1:2)] == colnames(counts_norm)[-(1:2)])
SanityCheck1_6 = all(colnames(counts)[-(c((1:2),which(colnames(counts) == 'spike_in_plasmid_DNA_from_CP1608')))]
                     == colnames(lfc)[-(1:2)])
SanityCheck1_7 = all(lfc_library == 0)
SanityCheck1_8 = all(normANDfcs$norm_counts$sgRNA == counts_norm$sgRNA)
SanityCheck1_9 = all(normANDfcs$logFCs$sgRNA == lfc$sgRNA)

# SanityCheck1_1 # Does the data frame with normalized counts contain the right amount of conditions?
# SanityCheck1_2 # Does the data frame with log fold changes contain the right amount of conditions?
# SanityCheck1_3 # Does the data frame with normalized counts contain the right amount of sgRNAs?
# SanityCheck1_4 # Does the data frame with log fold changes contain the right amount of sgRNAs?
# SanityCheck1_5 # Does the data frame with normalized counts contain all conditions?
# SanityCheck1_6 # Does the data frame with log fold changes contain all conditions?
# SanityCheck1_7 # Are the calculated log fold changes for the library all 0?
# SanityCheck1_8 # Do the generated normalised read counts correspond to the predefined sgRNAs?
# SanityCheck1_9 # Do the generated LFCs correspond to the predefinded sgRNAs?
all(SanityCheck1_1,SanityCheck1_2,SanityCheck1_3,SanityCheck1_4,SanityCheck1_5,SanityCheck1_6, SanityCheck1_7,
    SanityCheck1_8,SanityCheck1_9)

## Write files with normalised read counts and log fold changes
write.xlsx(lfc, 'log fold changes.xlsx')
write.xlsx(counts_norm, 'normalized read counts.xlsx')

## Sanity Checks for correction
SanityCheck2_1 = ncol(counts) == ncol(counts_cor)
SanityCheck2_2 = ncol(counts) == ncol(lfc_cor) + 1
SanityCheck2_3 = nrow(counts) == nrow(counts_cor)
SanityCheck2_4 = nrow(counts) == (nrow(lfc_cor) + nrow(nosites))
SanityCheck2_5 = all(colnames(counts)[-(1:2)] == colnames(counts_cor)[-(1:2)])
SanityCheck2_6 = all(colnames(counts)[-(c((1:2),which(colnames(counts) == 'spike_in_plasmid_DNA_from_CP1608')))]
                     == colnames(lfc_cor)[-(1:2)])
SanityCheck2_7 = all(lfc_cor_library == 0)
SanityCheck2_8 = all(correctedCounts$sgRNA == counts_cor$sgRNA)
SanityCheck2_9 = all(rownames(correctedFCs) == lfc_cor$sgRNA)

# SanityCheck2_1 # Does the data frame with corrected read counts contain the right amount of conditions?
# SanityCheck2_2 # Does the data frame with corrected log fold changes contain the right amount of conditions?
# SanityCheck2_3 # Does the data frame with corrected read counts contain the right amount of sgRNAs?
# SanityCheck2_4 # Does the data frame with corrected log fold changes contain the right amount of sgRNAs?
# SanityCheck2_5 # Does the data frame with corrected read counts contain all conditions?
# SanityCheck2_6 # Does the data frame with corrected log fold changes contain all conditions?
# SanityCheck2_7 # Are the calculated corrected log fold changes for the library all 0?
# SanityCheck2_8 # Do the generated corrected read counts correspond to the predefined sgRNAs?
# SanityCheck2_9 # Do the generated corrected log fold changes correspond to the predefined sgRNAS?
all(SanityCheck2_1,SanityCheck2_2,SanityCheck2_3,SanityCheck2_4,SanityCheck2_5,SanityCheck2_6, SanityCheck2_7,
    SanityCheck2_8,SanityCheck2_9)

## Write files with corrected read counts and log fold changes
write.xlsx(lfc_cor, 'corrected log fold changes.xlsx')
write.xlsx(counts_cor, 'corrected read counts.xlsx')

## Alternatively read in finished files for downstream manipulaiton
## Read files at the start, skip code manipulating data, read these files
# lfc = read.xlsx('log fold changes.xlsx')
# counts_norm = read.xlsx('normalized read counts.xlsx')
# lfc_cor = read.xlsx('corrected log fold changes.xlsx')
# counts_cor = read.xlsx('corrected read counts.xlsx')
# nosites = subset(lfc, substr(gene, 1, 8) == 'NO_SITE_')

## LFCs on gene level

lfc_gene = aggregate(lfc[,-(1:2)], list(lfc$gene), mean)
colnames(lfc_gene)[1] = 'gene'
lfc_cor_gene = aggregate(lfc_cor[,-(1:2)], list(lfc_cor$gene), mean)
colnames(lfc_cor_gene)[1] = 'gene'

## Sanity Checks
SanityCheck3_1 = ncol(lfc_gene) == (ncol(lfc) - 1)
SanityCheck3_2 = ncol(lfc_cor_gene) == (ncol(lfc_cor) - 1)
SanityCheck3_3 = nrow(lfc_gene) == (nrow(lfc_cor_gene) + length(unique(nosites$gene)))
SanityCheck3_4 = all(lfc_gene$gene[order(lfc_gene$gene)] == c(lfc_cor_gene$gene, unique(nosites$gene))
                     [order(c(lfc_cor_gene$gene, unique(nosites$gene)))])
SanityCheck3_5 = all(colnames(lfc_gene) == colnames(lfc)[-1])
SanityCheck3_6 = all(colnames(lfc_cor_gene) == colnames(lfc_cor)[-1])
SanityCheck3_7 = all(lfc_cor_gene$gene == unique(CP1608.anno.clean$GENES[order(CP1608.anno.clean$GENES)])) &&
  all(lfc_gene$gene == unique(CP1608.anno$GENES[order(CP1608.anno$GENES)]))

# SanityCheck3_1 # Does the data frame with uncorrected LFCs on gene level contain the right amount of conditions?
# SanityCheck3_2 # Does the data frame with corrected LFCs on gene level contain the right amount of conditions?
# SanityCheck3_3 # Are the numbers of genes equal for corrected and uncorrected LFCs on gene level?
# SanityCheck3_4 # Are the gene names the same for corrected and uncorrected LFCs on gene level?
# SanityCheck3_5 # Are the names of the conditions the same for the uncorrected LFCs on gene and on guide level?
# SanityCheck3_6 # Are the names of the conditions the same for the corrected LFC on gene and on guide level?
# SanityCheck3_7 # Are all genes present in both data frames for LFCs on gene level?
all(SanityCheck3_1,SanityCheck3_2,SanityCheck3_3,SanityCheck3_4,SanityCheck3_5,SanityCheck3_6,SanityCheck3_7)

## Write files with LFCs on gene level
write.xlsx(lfc_gene, 'uncorrected log fold changes on gene level.xlsx')
write.xlsx(lfc_cor_gene, 'corrected log fold changes on gene level.xlsx')

library(psych)

lfc_corrected = read.xlsx('directory/corrected log fold changes.xlsx')
lfc_uncorrected = read.xlsx('directory/log fold changes.xlsx')
control_essentials = read.csv('directory/Control_essentials.csv', header = FALSE)
genelfc_cor = read.xlsx('directory/corrected log fold changes on gene level.xlsx')
genelfc_uncor = read.xlsx('directory/uncorrected log fold changes on gene level.xlsx')
reads_corrected = read.xlsx('directory/corrected read counts.xlsx')
reads_uncorrected = read.xlsx('directory/read counts above thresholds.xlsx')
runs = read.csv2('directory/information on runs for conditions above threshold.csv', header = TRUE)
rownames(runs) = runs[,1]
runs = runs[,-1]
colnames(runs) = sub('\\.', '-', colnames(runs))
colnames(runs) = sub('\\.', '-', colnames(runs))
runs = runs[,-(which(colnames(runs) == "spike_in_plasmid_DNA_from_CP1608"))]
## Make subsets for essentials and nonessentials, nosites and one intergenic as subsets of nonessentials
## All subsets from corrected as well as uncorrected LFCs
essentials_cor = subset(lfc_corrected, gene %in% control_essentials$V1)
essentials_uncor = subset(lfc_uncorrected, gene %in% control_essentials$V1)
nosites = subset(lfc_uncorrected, substr(gene, 1, 8) == 'NO_SITE_')
oneinter_cor = subset(lfc_corrected, substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
oneinter_uncor = subset(lfc_uncorrected, substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
nonessentials_cor = subset(lfc_corrected, substr(gene, 1, 8) == 'NO_SITE_'
                           | substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
nonessentials_uncor = subset(lfc_uncorrected, substr(gene, 1, 8) == 'NO_SITE_'
                             | substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
genes_cor = subset(lfc_corrected, substr(gene, 1, 3) == 'MIR')
genes_uncor = subset(lfc_uncorrected, substr(gene, 1, 3) == 'MIR')

## Plots of distributions of LFC
## Compute means of LFC for each sg RNA
## For corrected LFCs
messentials_cor = data.frame(sgRNA = essentials_cor$sgRNA, type = rep('Essential genes', nrow(essentials_cor)),
                             mean_lfc = apply(essentials_cor[,-(1:2)], 1, mean))
moneinter_cor = data.frame(sgRNA = oneinter_cor$sgRNA, type = rep('One intergenic site sgRNAs', nrow(oneinter_cor)),
                           mean_lfc = apply(oneinter_cor[,-(1:2)], 1, mean))
# And for uncorrected LFCs
messentials_uncor = data.frame(sgRNA = essentials_uncor$sgRNA, type = rep('Essential genes', nrow(essentials_uncor)),
                               mean_lfc = apply(essentials_uncor[,-(1:2)], 1, mean))
mnosites = data.frame(sgRNA = nosites$sgRNA, type = rep('No site sgRNAs', nrow(nosites)),
                      mean_lfc = apply(nosites[,-(1:2)], 1, mean))
moneinter_uncor = data.frame(sgRNA = oneinter_uncor$sgRNA,
                             type = rep('One intergenic site sgRNAs', nrow(oneinter_uncor)),
                             mean_lfc = apply(oneinter_uncor[,-(1:2)], 1, mean))
## Create the plots
pdf('Distributions_Corrected_LFC.pdf')
plot(density(messentials_cor$mean_lfc), xlim = c(-6,1.2), ylim = c(0, 2), xlab = 'logfold change', col = 'transparent',
     main = 'Distributions of Corrected log FC for Control sgRNAs', bty = 'n')
polygon(density(messentials_cor$mean_lfc), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(moneinter_cor$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(messentials_cor$mean_lfc), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(moneinter_cor$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 2, c('Essential genes', 'One intergenic'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0)
dev.off()
pdf('Distributions_Uncorrected_LFC.pdf')
plot(density(messentials_uncor$mean_lfc), xlim = c(-5.6,1.3), ylim = c(0, 2.2), xlab = 'logfold change',
     main = 'Distributions of Uncorrected log FC for Control sgRNAs', bty = 'n', col = 'transparent')
polygon(density(messentials_uncor$mean_lfc), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(mnosites$mean_lfc), col = rgb(0.7291, 0.6588, 0.4549, alpha = .7), border = 'transparent')
polygon(density(moneinter_uncor$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(messentials_uncor$mean_lfc), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(mnosites$mean_lfc), col = rgb(0.7291, 0.6588, 0.4549), lty = 2, lw = 3)
abline(v = median(moneinter_uncor$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-5.6, 2.1, c('Essential genes', 'One intergenic', 'Nosites'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549)),
                 alpha = c(.6,.6,.7)), box.lty = 0)
dev.off()

## Calculate NNMD as measure for distinction between positive and negative controls
NNMD_cor = ((apply(essentials_cor[,-(1:2)], 2, mean) - apply(oneinter_cor[,-(1:2)], 2, mean))
            /apply(oneinter_cor[,-(1:2)], 2, sd))
NNMD_uncor = ((apply(essentials_uncor[,-(1:2)], 2, mean) - apply(oneinter_uncor[,-(1:2)], 2, mean))
              /apply(oneinter_uncor[,-(1:2)], 2, sd))

## Plot distribution of NNMD
pdf('Distribution_of_NNMD_all_replicates_corrected.pdf')
hist(NNMD_cor, freq = FALSE, 40,xlim = c(-12,2), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .4),
     main = '', xlab = 'NNMD')
lines(density(NNMD_cor), lwd = 2, col = rgb(0.2117647, 0.3921569, 0.5450980))
rug(NNMD_cor)
# abline(v = -1, col = rgb(0.8, 0, 0), lw = 3)
dev.off()
pdf('Distribution_of_NNMD_all_replicates_uncorrected.pdf')
hist(NNMD_uncor, freq = FALSE, 40, xlim = c(-12,1), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .4),
     main = '', xlab = 'NNMD')
lines(density(NNMD_uncor), lwd = 2, col = rgb(0.2117647, 0.3921569, 0.5450980))
rug(NNMD_uncor)
# abline(v = -1, col = rgb(0.8, 0, 0), lw = 3)
dev.off()
pdf('Comparison_of_distributions_NNMD_corrected_vs_uncorrected.pdf')
plot(density(NNMD_cor), xlim = c(-12,2), main = '', xlab = 'NNMD', bty = 'n', col = 'transparent')
polygon(density(NNMD_cor), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(NNMD_uncor), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
legend(-12.5, 0.155, c('Before correction', 'After correction'),
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = 0.6), lty = 1, lwd = 5, box.lty = 0)
dev.off()

## Identify bad replicates based on NNMD
## Set threshold
threshold_NNMD = -1
## Split NNMD into subsets above and below threshold, not final decision, only to illustrate
good_NNMD_c = subset(NNMD_cor, NNMD_cor <= threshold_NNMD)
bad_NNMD_c = subset(NNMD_cor, NNMD_cor > threshold_NNMD)
good_NNMD_uc = subset(NNMD_uncor, NNMD_uncor <= threshold_NNMD)
bad_NNMD_uc = subset(NNMD_uncor, NNMD_uncor > threshold_NNMD)

## Calculate Cohen's D
## Transform LFCs in FCs
D_cor = (apply((2^oneinter_cor[,-(1:2)]), 2, mean) - apply((2^essentials_cor[,-(1:2)]), 2, mean)) / 
  sqrt((((nrow(oneinter_cor)-1)*apply((2^oneinter_cor[,-(1:2)]), 2, var))+
          (nrow(essentials_cor)-1)*apply((2^essentials_cor[,-(1:2)]), 2, var)) /
         (nrow(oneinter_cor)+nrow(essentials_cor)-2))
D_uncor = (apply((2^oneinter_uncor[,-(1:2)]), 2, mean) - apply((2^essentials_uncor[,-(1:2)]), 2, mean)) / 
  sqrt((((nrow(oneinter_uncor)-1)*apply((2^oneinter_uncor[,-(1:2)]), 2, var))+
          (nrow(essentials_uncor)-1)*apply((2^essentials_uncor[,-(1:2)]), 2, var)) /
         (nrow(oneinter_uncor)+nrow(essentials_uncor)-2))
## Remoce replicates based on Choen's D
##Set threshold
threshold_D = 1
good_D_c = subset(D_cor, D_cor > threshold_D)
bad_D_c = subset(D_cor, D_cor <= threshold_D)

pdf('Distribution_of_Cohens_D_all_replicates_corrected.pdf')
hist(D_cor, 45, freq = FALSE, xlim = c(-0.3, 4), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .4),
     main = '', xlab = "Cohen's D")
lines(density(D_cor), lwd = 2, col = rgb(0.2117647, 0.3921569, 0.5450980))
rug(D_cor)
# abline(v = 1, col = rgb(0.8, 0, 0), lw = 3)
dev.off()

pdf('Distribution_of_Cohens_D_all_replicates_uncorrected.pdf')
hist(D_uncor, freq = FALSE, 45, xlim = c(-0.3,4), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .4),
     main = '', xlab = "Cohen's D")
lines(density(D_uncor), lwd = 2, col = rgb(0.2117647, 0.3921569, 0.5450980))
rug(D_uncor)
# abline(v = 1, col = rgb(0.8, 0, 0), lw = 3)
dev.off()

pdf('Comparison_of_distributions_Cohens_D_corrected_vs_uncorrected.pdf')
plot(density(D_uncor), xlim = c(-0.5,4.5), main = '', xlab = "Cohen's D", bty = 'n', col = 'transparent')
polygon(density(D_uncor), col = rgb(0.8, 0, 0, alpha = 0.6), border = 'transparent')
polygon(density(D_cor), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = 0.6), border = 'transparent')
legend(-0.5, 0.61, c('Before correction', 'After correction'),
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), lty = 1, lwd = 5,
       box.lty = 0, bty = 'n')
dev.off()

## Sanity Checks for corrected LFCs
SanityCheck1_1 = nrow(lfc_corrected) == (nrow(essentials_cor) + nrow(nonessentials_cor) + nrow(genes_cor))
SanityCheck1_2 = all(lfc_corrected$gene[order(lfc_corrected$gene)] ==
                       c(essentials_cor$gene, nonessentials_cor$gene, genes_cor$gene)[order(
                         c(essentials_cor$gene,nonessentials_cor$gene, genes_cor$gene))])
SanityCheck1_3 = nrow(nonessentials_cor) == nrow(oneinter_cor)
SanityCheck1_4 = all(nonessentials_cor$gene[order(nonessentials_cor$gene)] == oneinter_cor$gene
                     [order(oneinter_cor$gene)])
SanityCheck1_5 = all.equal.character(colnames(lfc_corrected), colnames(essentials_cor), colnames(nosites_cor),
                                     colnames(oneinter_cor), colnames(nonessentials_cor), colnames(genes_cor))
SanityCheck1_6 = all(colnames(lfc_corrected)[-(1:2)] == names(NNMD_cor))
SanityCheck1_7 = all(names(NNMD_cor) == names(D_cor))

# SanityCheck1_1 # Do the subsets together have the same amount of rows as the original data frame?
# SanityCheck1_2 # Do the subsets together contain all sgRNAs?
# SanityCheck1_3 # Do the two subsets with different negative controls (nosites and oneinter) together have the same amount of rows as the subset with nonessential genes?
# SanityCheck1_4 # Do the two subsets with different negative controls (nosites and oneinter) contain all sgRNAs also contained in the subset with nonessentials?
# SanityCheck1_5 # Do all subsets and the original data frame have the same conditions in the same order?
# SanityCheck1_6 # Does NNMD contain all replicates?
# SanityCheck1_7 # Does Cohen's D contain all replicates?
all(SanityCheck1_1,SanityCheck1_2,SanityCheck1_3,SanityCheck1_4,SanityCheck1_5,SanityCheck1_6,
    SanityCheck1_7)

## Sanity checks for uncorrected LFCs
SanityCheck2_1 = nrow(lfc_uncorrected) == (nrow(essentials_uncor) + nrow(nonessentials_uncor) + nrow(genes_uncor))
SanityCheck2_2 = all(lfc_uncorrected$gene[order(lfc_uncorrected$gene)] ==
                       c(essentials_uncor$gene, nonessentials_uncor$gene, genes_uncor$gene)[order(
                         c(essentials_uncor$gene,nonessentials_uncor$gene, genes_uncor$gene))])
SanityCheck2_3 = nrow(nonessentials_uncor) == (nrow(nosites) + nrow(oneinter_uncor))
SanityCheck2_4 = all(nonessentials_uncor$gene[order(nonessentials_uncor$gene)] == c(nosites$gene, oneinter_uncor$gene)
                     [order(c(nosites$gene, oneinter_uncor$gene))])
SanityCheck2_5 = all.equal.character(colnames(lfc_uncorrected), colnames(essentials_uncor), colnames(nosites_uncor),
                                     colnames(oneinter_uncor), colnames(nonessentials_uncor), colnames(genes_uncor))
SanityCheck2_6 = all(colnames(lfc_uncorrected)[-(1:2)] == names(NNMD_uncor))
SanityCheck2_7 = all(names(NNMD_uncor) == names(D_uncor))
SanityCheck2_8 = all(names(NNMD_cor) == names(NNMD_uncor))

# SanityCheck2_1 # Do the subsets together have the same amount of rows as the original data frame?
# SanityCheck2_2 # Do the subsets together contain all sgRNAs?
# SanityCheck2_3 # Do the two subsets with different negative controls (nosites and oneinter) together have the same amount of rows as the subset with nonessential genes?
# SanityCheck2_4 # Do the two subsets with different negative controls (nosites and oneinter) contain all sgRNAs also contained in the subset with nonessentials?
# SanityCheck2_5 # Do all subsets and the original data frame have the same conditions in the same order?
# SanityCheck2_6 # Does NNMD contain all replicates?
# SanityCheck2_7 # Does Cohen's D contain all replicates?
# SanityCheck2_8 # Are the replicate names the same for corrected and uncorrected measures?
all(SanityCheck2_1,SanityCheck2_2,SanityCheck2_3,SanityCheck2_4,SanityCheck2_5,SanityCheck2_6,
    SanityCheck2_7,SanityCheck2_8)

## F measures
files = list.files() # all files in the folder
F_mes = data.frame() # prepare data frame to hold F measures
for (ii in 1:length(files)) { # for all files in the folder CAUTION: ii does not count txt files but all files! 
  if (grepl('.txt', files[ii], fixed = TRUE)) { # only the txt files
    temp = read.table(files[ii], header = TRUE) # read file
    tempBF = subset(temp, BF >= 5) # only BF over 5
    if (tempBF$Recall[nrow(tempBF)] == 0 && tempBF$Precision[nrow(tempBF)] == 0){ # harmonic mean defined 0 if one input is 0
      F_mes[nrow(F_mes)+1,1] = 0
    }
    else {
      F_mes[nrow(F_mes)+1,1] = 2*(tempBF$Recall[nrow(tempBF)]*tempBF$Precision[nrow(tempBF)])/
        (tempBF$Recall[nrow(tempBF)]+tempBF$Precision[nrow(tempBF)])
      
    }
    rownames(F_mes)[length(rownames(F_mes))] = sub('_pr\\.txt', '', files[ii])
  }
}


pdf('Distribution_of_F_Measure.pdf')
hist(F_mes$V1, freq = FALSE, 60, xlim = c(0,1.1), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .4),
     main = '', xlab = 'F measure')
lines(density(F_mes$V1), col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 2)
rug(F_mes$V1)
dev.off()

## Build data frame with measures of separation of gene controls
separation = cbind(NNMD_uncor, NNMD_cor, D_uncor, D_cor)
colnames(separation) = c('NNMD before correction', 'NNMD after correction',
                         "Cohen's D before correction", "Cohen's D after correction")
write.xlsx(as.data.frame(separation), 'Measures of separation of gene controls replicates.xlsx', rowNames = TRUE)
cl_separation = cbind.data.frame(separation, sub('_REP_.', '', rownames(separation)))
colnames(cl_separation)[5] = 'Cell_Line'
cl_separation = aggregate(cl_separation[,-5], list(cl_separation$Cell_Line), mean)
colnames(cl_separation)[1] = 'Cell_Line'
cl_separation = cbind(cl_separation, F_mes)
colnames(cl_separation)[6] = 'F'
bad_sep = subset(cl_separation, Cell_Line == 'Jurkat' | Cell_Line == 'MDA-MB453')
good_sep = subset(cl_separation, Cell_Line != 'Jurkat' & Cell_Line != 'MDA-MB453')

pdf('Correlation_Cohen_D_NNMD.pdf')
plot(good_sep$`NNMD after correction`, good_sep$`Cohen's D after correction`, xlim = c(-9,0), ylim = c(0,3.5),
     pch = 16, xlab = 'NNMD', ylab = "Cohen's D", bty = 'n')
points(bad_sep$`NNMD after correction`, bad_sep$`Cohen's D after correction`, col = rgb(0.8, 0, 0),
       pch = 16)
# lines(lowess(cl_separation$`NNMD after correction`, cl_separation$`Cohen's D after correction`, f = 4/5),
# col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 3)
abline(lm(cl_separation$`Cohen's D after correction` ~ cl_separation$`NNMD after correction`),
       col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 3)
text(bad_sep$`NNMD after correction`[1], bad_sep$`Cohen's D after correction`[1]*0.82, 'Jurkat', adj = 0.5)
text(bad_sep$`NNMD after correction`[2]*1.15, bad_sep$`Cohen's D after correction`[2], 'MDA-MB453', adj = 1)
dev.off()

pdf('Correlation_Cohen_D_F_Measure.pdf')
plot(good_sep$`Cohen's D after correction`, good_sep$F, xlim = c(0,3), ylim = c(0,1), pch = 16, xlab = "Cohen's D",
     ylab = 'F Measure', bty = 'n')
points(bad_sep$`Cohen's D after correction`, bad_sep$F, col = rgb(0.8, 0, 0), pch = 16)
# lines(lowess(cl_separation$`Cohen's D after correction`, cl_separation$F, f = 2/5),
#       col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 3)
abline(lm(cl_separation$F ~ cl_separation$`Cohen's D after correction`),
       col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 3)
text(bad_sep$`Cohen's D after correction`[1]*1.065, bad_sep$F[1], 'Jurkat', adj = 0)
text(bad_sep$`Cohen's D after correction`[2]*1.065, bad_sep$F[2], 'MDA-MB453', adj = 0)
dev.off()

pdf('Correlation_NNMD_F_Measure.pdf')
plot(good_sep$`NNMD after correction`, good_sep$F, xlim = c(-9,0), ylim = c(0,1), pch = 16, xlab = 'NNMD',
     ylab = 'F Measure', bty = 'n')
points(bad_sep$`NNMD after correction`, bad_sep$F, pch = 16, col = rgb(0.8, 0, 0))
# lines(lowess(cl_separation$`NNMD after correction`, cl_separation$F, f = 3/5),
# col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 3)
abline(lm(cl_separation$F ~ cl_separation$`NNMD after correction`),
       col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 3)
text(bad_sep$`NNMD after correction`[1]*1.2, bad_sep$F[1], 'Jurkat', adj = 1)
text(bad_sep$`NNMD after correction`[2]*1.15, bad_sep$F[2], 'MDA-MB453', adj = 1)
dev.off()
## Remove Jurkat and Ben_MEN
lfc_corrected_rem = lfc_corrected[,-c(which(sub('_REP_.', '', colnames(lfc_corrected)) == 'Jurkat'),
                                      which(sub('_REP_.', '', colnames(lfc_corrected)) == 'MDA-MB453'))]
lfc_uncorrected_rem = lfc_uncorrected[,-c(which(sub('_REP_.', '', colnames(lfc_uncorrected)) == 'Jurkat'),
                                          which(sub('_REP_.', '', colnames(lfc_uncorrected)) == 'MDA-MB453'))]
separation_rem = separation[-c(which(sub('_REP_.', '', rownames(separation)) == 'Jurkat'),
                               which(sub('_REP_.', '', rownames(separation)) == 'MDA-MB453')),]
separation_bad = separation[c(which(sub('_REP_.', '', rownames(separation)) == 'Jurkat'),
                              which(sub('_REP_.', '', rownames(separation)) == 'MDA-MB453')),]
runs_rem = runs[,-c(which(sub('_REP_.', '', colnames(runs)) == 'Jurkat'),
                    which(sub('_REP_.', '', colnames(runs)) == 'MDA-MB453'))]
genelfc_cor_rem = genelfc_cor[,-c(which(sub('_REP_.', '', colnames(genelfc_cor)) == 'Jurkat'),
                                  which(sub('_REP_.', '', colnames(genelfc_cor)) == 'MDA-MB453'))]
genelfc_uncor_rem = genelfc_uncor[,-c(which(sub('_REP_.', '', colnames(genelfc_uncor)) == 'Jurkat'),
                                      which(sub('_REP_.', '', colnames(genelfc_uncor)) == 'MDA-MB453'))]
reads_cor_rem = reads_corrected[,-c(which(sub('_REP_.', '', colnames(reads_corrected)) == 'Jurkat'),
                                    which(sub('_REP_.', '', colnames(reads_corrected)) == 'MDA-MB453'))]
reads_uncor_rem = reads_uncorrected[,-c(which(sub('_REP_.', '', colnames(reads_uncorrected)) == 'Jurkat'),
                                        which(sub('_REP_.', '', colnames(reads_uncorrected)) == 'MDA-MB453'))]

## Sanity Checks for removal
SanityCheck3_1 = (ncol(lfc_corrected_rem) == ncol(lfc_corrected) - 5) && 
  (ncol(lfc_uncorrected_rem) == ncol(lfc_uncorrected) - 5) &&
  (nrow(separation_rem) == nrow(separation) - 5) && (ncol(runs_rem) == ncol(runs) - 5) &&
  (ncol(genelfc_cor_rem) == ncol(genelfc_cor) - 5) && (ncol(genelfc_uncor_rem) == ncol(genelfc_uncor) - 5) &&
  (ncol(reads_cor_rem) == ncol(reads_corrected) - 5) && (ncol(reads_uncor_rem) == ncol(reads_uncorrected) - 5)
SanityCheck3_2 = all(all(lfc_corrected_rem != 'Jurkat_REP_A') && all(lfc_corrected_rem != 'Jurkat_REP_B') &&
                       all(lfc_corrected_rem != 'Jurkat_REP_C') && all(lfc_corrected_rem != 'MDA-MB453_REP_A') &&
                       all(lfc_corrected_rem != 'MDA-MB453_REP_B')) &&
  all(all(lfc_uncorrected_rem != 'Jurkat_REP_A') && all(lfc_uncorrected_rem != 'Jurkat_REP_B') &&
        all(lfc_uncorrected_rem != 'Jurkat_REP_C') && all(lfc_uncorrected_rem != 'MDA-MB453_REP_A') &&
        all(lfc_uncorrected_rem != 'MDA-MB453_REP_B')) &&
  all(all(separation_rem != 'Jurkat_REP_A') && all(separation_rem != 'Jurkat_REP_B') &&
        all(separation_rem != 'Jurkat_REP_C') && all(separation_rem != 'MDA-MB453_REP_A') &&
        all(separation_rem != 'MDA-MB453_REP_B')) &&
  all(all(runs_rem != 'Jurkat_REP_A') && all(runs_rem != 'Jurkat_REP_B') &&
        all(runs_rem != 'Jurkat_REP_C') && all(runs_rem != 'MDA-MB453_REP_A') &&
        all(runs_rem != 'MDA-MB453_REP_B')) &&
  all(all(genelfc_cor_rem != 'Jurkat_REP_A') && all(genelfc_cor_rem != 'Jurkat_REP_B') &&
        all(genelfc_cor_rem != 'Jurkat_REP_C') && all(genelfc_cor_rem != 'MDA-MB453_REP_A') &&
        all(genelfc_cor_rem != 'MDA-MB453_REP_B')) &&
  all(all(genelfc_uncor_rem != 'Jurkat_REP_A') && all(genelfc_uncor_rem != 'Jurkat_REP_B') &&
        all(genelfc_uncor_rem != 'Jurkat_REP_C') && all(genelfc_uncor_rem != 'MDA-MB453_REP_A') &&
        all(genelfc_uncor_rem != 'MDA-MB453_REP_B')) &&
  all(all(reads_cor_rem != 'Jurkat_REP_A') && all(reads_cor_rem != 'Jurkat_REP_B') &&
        all(reads_cor_rem != 'Jurkat_REP_C') && all(reads_cor_rem != 'MDA-MB453_REP_A') &&
        all(reads_cor_rem != 'MDA-MB453_REP_B')) &&
  all(all(reads_uncor_rem != 'Jurkat_REP_A') && all(reads_uncor_rem != 'Jurkat_REP_B') &&
        all(reads_uncor_rem != 'Jurkat_REP_C') && all(reads_uncor_rem != 'MDA-MB453_REP_A') &&
        all(reads_uncor_rem != 'MDA-MB453_REP_B'))
# SanityCheck3_1 # Is the right amount of replicates removed from all data frames?
# SanityCheck3_2 # Are all bad replicates removed?
all(SanityCheck3_1,SanityCheck3_2)

write.xlsx(lfc_corrected_rem, 'Corrected LFC without bad cell lines.xlsx')
write.xlsx(lfc_uncorrected_rem, 'Uncorrected LFC without bad cell lines.xlsx')
write.xlsx(good_sep, 'Measures of separation of gene controls good cell lines.xlsx')
write.xlsx(bad_sep, 'Measures of separation of gene controls excluded cell lines.xlsx')
write.xlsx(runs_rem, 'Information on runs good cell lines.xlsx', rowNames = TRUE)
write.xlsx(genelfc_cor_rem, 'Corrected LFC on gene level without bad cell lines.xlsx')
write.xlsx(genelfc_uncor_rem, 'Uncorrected LFC on gene level without bad cell lines.xlsx')
write.xlsx(reads_cor_rem, 'Corrected read counts without bad cell lines.xlsx')
write.xlsx(reads_uncor_rem, 'Uncorrected read counts without bad cell lines.xlsx')

## Gene-level read counts
reads_c_gene = aggregate(reads_cor_rem[,-(1:2)], list(reads_cor_rem$gene), mean)
colnames(reads_c_gene)[1] = 'gene'
reads_uc_gene = aggregate(reads_uncor_rem[,-(1:2)], list(reads_cor_rem$gene), mean)
colnames(reads_uc_gene)[1] = 'gene'
reads_c_gene_wo_nosites = reads_c_gene[-which(grepl('NO_SITE', reads_c_gene$gene)),]

SanityCheck4_1 = all(substr(reads_cor_rem$sgRNA, 1, 20) == reads_uncor_rem$Construct_Barcode)
SanityCheck4_2 = all(genelfc_uncor_rem$gene == reads_c_gene$gene) && all(genelfc_uncor_rem$gene == reads_uc_gene$gene)
SanityCheck4_3 = all(genelfc_cor_rem$gene == reads_c_gene_wo_nosites$gene)

# SanityCheck4_1 # Are the guides in the same order for both read count data frames? If yes, gene names can be assigned from the other data frame.
# SanityCheck4_2 # Are the genes correct?
# SanityCheck4_3 # Are the genes correct for the corrected reads without nosites?
all(SanityCheck4_1,SanityCheck4_2,SanityCheck4_3)

write.xlsx(reads_c_gene, 'Corrected read counts on gene level with nosites.xlsx')
write.xlsx(reads_uc_gene, 'Unorrected read counts on gene level.xlsx')
write.xlsx(reads_c_gene_wo_nosites, 'Corrected read counts on gene level without nosites.xlsx')

## Distributions without excluded replicates
essentials_cor_rem = subset(lfc_corrected_rem, gene %in% control_essentials$V1)
essentials_uncor_rem = subset(lfc_uncorrected_rem, gene %in% control_essentials$V1)
nosites_rem = subset(lfc_uncorrected_rem, substr(gene, 1, 8) == 'NO_SITE_')
oneinter_cor_rem = subset(lfc_corrected_rem, substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
oneinter_uncor_rem = subset(lfc_uncorrected_rem, substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
nonessentials_cor_rem = subset(lfc_corrected_rem, substr(gene, 1, 8) == 'NO_SITE_'
                               | substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
nonessentials_uncor_rem = subset(lfc_uncorrected_rem, substr(gene, 1, 8) == 'NO_SITE_'
                                 | substr(gene, 1, 20) == 'ONE_INTERGENIC_SITE_')
genes_cor_rem = subset(lfc_corrected_rem, substr(gene, 1, 3) == 'MIR')
genes_uncor_rem = subset(lfc_uncorrected_rem, substr(gene, 1, 3) == 'MIR')

## Plots of distributions of LFC
## Compute means of LFC for each sg RNA
## For corrected LFCs
messentials_cor_rem = data.frame(sgRNA = essentials_cor_rem$sgRNA, type = rep('Essential genes',
                                                                              nrow(essentials_cor_rem)),
                                 mean_lfc = apply(essentials_cor_rem[,-(1:2)], 1, mean))
moneinter_cor_rem = data.frame(sgRNA = oneinter_cor_rem$sgRNA, type = rep('One intergenic site sgRNAs',
                                                                          nrow(oneinter_cor_rem)),
                               mean_lfc = apply(oneinter_cor_rem[,-(1:2)], 1, mean))
# And for uncorrected LFCs
messentials_uncor_rem = data.frame(sgRNA = essentials_uncor_rem$sgRNA, type = rep('Essential genes',
                                                                                  nrow(essentials_uncor_rem)),
                                   mean_lfc = apply(essentials_uncor_rem[,-(1:2)], 1, mean))
mnosites_rem = data.frame(sgRNA = nosites_rem$sgRNA, type = rep('No site sgRNAs', nrow(nosites_rem)),
                          mean_lfc = apply(nosites_rem[,-(1:2)], 1, mean))
moneinter_uncor_rem = data.frame(sgRNA = oneinter_uncor_rem$sgRNA,
                                 type = rep('One intergenic site sgRNAs', nrow(oneinter_uncor_rem)),
                                 mean_lfc = apply(oneinter_uncor_rem[,-(1:2)], 1, mean))
## Create the plots
pdf('Distributions_Corrected_LFC_wo_Removed_Replicates.pdf')
plot(density(messentials_cor_rem$mean_lfc), xlim = c(-6,1.2), ylim = c(0, 2), xlab = 'logfold change',
     col = 'transparent', main = 'Distributions of Corrected log FC for Control sgRNAs', bty = 'n')
polygon(density(messentials_cor_rem$mean_lfc), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(moneinter_cor_rem$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(messentials_cor_rem$mean_lfc), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(moneinter_cor_rem$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-6, 2, c('Essential genes', 'One intergenic'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980)), alpha = .6), box.lty = 0)
dev.off()
pdf('Distributions_Uncorrected_LFC_wo_Removed_Replicates.pdf')
plot(density(messentials_uncor_rem$mean_lfc), xlim = c(-5.6,1.3), ylim = c(0, 2.2), xlab = 'logfold change',
     main = 'Distributions of Uncorrected log FC for Control sgRNAs', bty = 'n', col = 'transparent')
polygon(density(messentials_uncor_rem$mean_lfc), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(mnosites_rem$mean_lfc), col = rgb(0.7291, 0.6588, 0.4549, alpha = .7), border = 'transparent')
polygon(density(moneinter_uncor_rem$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(messentials_uncor_rem$mean_lfc), col = rgb(0.8, 0, 0), lty = 2, lw = 3)
abline(v = median(mnosites_rem$mean_lfc), col = rgb(0.7291, 0.6588, 0.4549), lty = 2, lw = 3)
abline(v = median(moneinter_uncor_rem$mean_lfc), col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
legend(-5.6, 2.1, c('Essential genes', 'One intergenic', 'Nosites'), lty = 1, lw = 5,
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549)),
                 alpha = c(.6,.6,.7)), box.lty = 0)
dev.off()

library(pheatmap)
library(dendextend)
library(umap)
library(RColorBrewer)
library(stringr)
library(ggplot2)
library(webr)
library(dplyr)
library(tidyr)

genelfc_cor = read.xlsx('directory/Corrected LFC on gene level without bad cell lines.xlsx')
genelfc_uncor = read.xlsx('directory/Uncorrected LFC on gene level without bad cell lines.xlsx')
runs = read.xlsx('directory/Information on runs good cell lines.xlsx', rowNames = TRUE)
lfc_cor = read.xlsx('directory/Corrected LFC without bad cell lines.xlsx')
lfc_uncor = read.xlsx('directory/Uncorrected LFC without bad cell lines.xlsx')
rc_cor = read.xlsx('directory/Corrected read counts without bad cell lines.xlsx')
rc_uncor = read.xlsx('directory/Uncorrected read counts without bad cell lines.xlsx')
lineage = read.xlsx('directory/B41_miRNA_screen_cell_line_overview_stripped.xlsx')

## Correlation matrix and clustering
## Sort Columns

sortedFCs_c = cbind(gene = genelfc_cor$gene, genelfc_cor[,(order(colnames(genelfc_cor)[-1]) + 1)])
sortedFCs_uc = cbind(gene = genelfc_uncor$gene, genelfc_uncor[,(order(colnames(genelfc_uncor)[-1]) + 1)])

## Identify the top 3 % most variable scores
## Calculate variance for all genes
variance_c = cbind(sortedFCs_c, variance = apply(sortedFCs_c[,-1], 1, var))
variance_uc = cbind(sortedFCs_uc, variance = apply(sortedFCs_uc[,-1], 1, var))

## Sort genes according to variance
sort_variance_c = variance_c[order(variance_c$variance),]
sort_variance_uc = variance_uc[order(variance_uc$variance),]

set_fraction = 0.03
high_variance_c = sort_variance_c[(nrow(sort_variance_c):(set_fraction*nrow(sort_variance_c))),]
high_variance_uc = sort_variance_uc[(nrow(sort_variance_uc):(set_fraction*nrow(sort_variance_uc))),]
table_c = cor(high_variance_c[,-c(1,ncol(high_variance_c))])
table_uc = cor(high_variance_uc[,-c(1,ncol(high_variance_uc))])
## Corrected LFCs
sorted_runs = runs[,order(colnames(runs))]
pheatmap(table_c, treeheight_row =  0, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(t(sorted_runs[1,])))
sorted_runs_cl = rbind(sorted_runs, sub('_REP_.', '', colnames(sorted_runs)))
rownames(sorted_runs_cl)[3] = 'Cell_Line'
pheatmap(table_uc, treeheight_row =  0, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(t(sorted_runs_cl[c(1,3),])), )
hclust_c = hclust(dist(table_c), method = "average")
pdf('Dendrogram_of_Replicates.pdf', height = 10, width = 20)
par(mar = c(7.1, 5.1, 4.1, 2.1))
as.dendrogram(hclust_c) %>%
  plot(horiz = F)
dev.off()

pdf('Heatmap_Replicates_after_Correction_Variable_Genes.pdf') # generate figures with heat maps
pheatmap(table_c, treeheight_row = 0, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(t(sorted_runs[1,])), clustering_method = 'average')
dev.off()
## Uncorrected LFCs
pdf('Heatmap_Replicates_before_Correction_Variable_Genes.pdf')
pheatmap(table_uc, treeheight_row = 0, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(t(sorted_runs[1,])), clustering_method = 'average')
dev.off()
## Use all genes instead of most variable
table_all_c = cor(sortedFCs_c[,-1])
table_all_uc = cor(sortedFCs_uc[,-1])
## Corrected LFCs
pdf('Heatmap_Replicates_after_Correction_all_Genes.pdf')
pheatmap(table_all_c, treeheight_row = 0, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(t(sorted_runs[1,])), clustering_method = 'average')
dev.off()
## Uncorrected LFCs
pdf('Heatmap_Replicates_before_Correction_all_Genes.pdf')
pheatmap(table_all_uc, treeheight_row = 0, show_rownames = FALSE, show_colnames = FALSE,
         annotation_col = as.data.frame(t(sorted_runs[1,])), clustering_method = 'average')
dev.off()
## Compare distributions
pdf('Distribution_of_Correlation_Coefficients_between_Replicates_of_all_Cell_Lines_Variable_genes.pdf')
plot(density(table_c), col = 'transparent', ylim = c(0, 5), main = 'Distributions of Correlation Coefficients', xlab = '', bty = 'n')
polygon(density(table_c), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(table_uc), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
legend(0.4, 5.2, c('Corrected LFC', 'Uncorrected LFC'), lty = 1, lwd = 5, box.lty = 0, bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980), alpha = 0.6)))
dev.off()

## Build data frame with correlation within cell lines only
## Highly variable genes, corrected LFC
lines_c = unique(sub('_REP_.', '', colnames(table_c)))
line_correlation_c = data.frame(Cell_line = c(), Cor_AB = c(), Cor_AC = c(), Cor_BC = c())
for (ii in 1:length(lines_c)) {
  if(any(colnames(table_c) == paste(lines_c[ii], 'REP_A', sep = '_'))){
    if(any(rownames(table_c) == paste(lines_c[ii], 'REP_B', sep = '_'))){
      AB = table_c[which(rownames(table_c) == paste(lines_c[ii], 'REP_B', sep = '_')),
                   which(colnames(table_c) == paste(lines_c[ii], 'REP_A', sep = '_'))]
    }
    else{
      AB = NA
    }
    if(any(rownames(table_c) == paste(lines_c[ii], 'REP_C', sep = '_'))){
      AC = table_c[which(rownames(table_c) == paste(lines_c[ii], 'REP_C', sep = '_')),
                   which(colnames(table_c) == paste(lines_c[ii], 'REP_A', sep = '_'))]
    }
    else{
      AC = NA
    }
  }
  else{
    AB = NA
    AC = NA
  }
  if(any(colnames(table_c) == paste(lines_c[ii], 'REP_B', sep = '_'))){
    if(any(rownames(table_c) == paste(lines_c[ii], 'REP_C', sep = '_'))){
      BC = table_c[which(rownames(table_c) == paste(lines_c[ii], 'REP_C', sep = '_')),
                   which(colnames(table_c) == paste(lines_c[ii], 'REP_B', sep = '_'))]
    }
    else{
      BC = NA
    }
  }
  else{
    BC = NA
  }
  line_correlation_c = rbind(line_correlation_c, cbind.data.frame(AB, AC, BC))
  rownames(line_correlation_c)[nrow(line_correlation_c)] = lines_c[ii]
}
colnames(line_correlation_c) = c('cor(AB)', 'cor(AC)', 'cor(BC)')
## Highly variable genes, uncorrected LFCs
lines_uc = unique(sub('_REP_.', '', colnames(table_uc)))
line_correlation_uc = data.frame(Cell_line = c(), Cor_AB = c(), Cor_AC = c(), Cor_BC = c())
for (ii in 1:length(lines_uc)) {
  if(any(colnames(table_uc) == paste(lines_uc[ii], 'REP_A', sep = '_'))){
    if(any(rownames(table_uc) == paste(lines_uc[ii], 'REP_B', sep = '_'))){
      AB = table_uc[which(rownames(table_uc) == paste(lines_uc[ii], 'REP_B', sep = '_')),
                    which(colnames(table_uc) == paste(lines_uc[ii], 'REP_A', sep = '_'))]
    }
    else{
      AB = NA
    }
    if(any(rownames(table_uc) == paste(lines_uc[ii], 'REP_C', sep = '_'))){
      AC = table_uc[which(rownames(table_uc) == paste(lines_uc[ii], 'REP_C', sep = '_')),
                    which(colnames(table_uc) == paste(lines_uc[ii], 'REP_A', sep = '_'))]
    }
    else{
      AC = NA
    }
  }
  else{
    AB = NA
    AC = NA
  }
  if(any(colnames(table_uc) == paste(lines_uc[ii], 'REP_B', sep = '_'))){
    if(any(rownames(table_uc) == paste(lines_uc[ii], 'REP_C', sep = '_'))){
      BC = table_uc[which(rownames(table_uc) == paste(lines_uc[ii], 'REP_C', sep = '_')),
                    which(colnames(table_uc) == paste(lines_uc[ii], 'REP_B', sep = '_'))]
    }
    else{
      BC = NA
    }
  }
  else{
    BC = NA
  }
  line_correlation_uc = rbind(line_correlation_uc, cbind.data.frame(AB, AC, BC))
  rownames(line_correlation_uc)[nrow(line_correlation_uc)] = lines_uc[ii]
}
colnames(line_correlation_uc) = c('cor(AB)', 'cor(AC)', 'cor(BC)')
## All genes, corrected LFC
lines_all_c = unique(sub('_REP_.', '', colnames(table_all_c)))
line_correlation_all_c = data.frame(Cell_line = c(), Cor_AB = c(), Cor_AC = c(), Cor_BC = c())
for (ii in 1:length(lines_all_c)) {
  if(any(colnames(table_all_c) == paste(lines_all_c[ii], 'REP_A', sep = '_'))){
    if(any(rownames(table_all_c) == paste(lines_all_c[ii], 'REP_B', sep = '_'))){
      AB = table_all_c[which(rownames(table_all_c) == paste(lines_all_c[ii], 'REP_B', sep = '_')),
                       which(colnames(table_all_c) == paste(lines_all_c[ii], 'REP_A', sep = '_'))]
    }
    else{
      AB = NA
    }
    if(any(rownames(table_all_c) == paste(lines_all_c[ii], 'REP_C', sep = '_'))){
      AC = table_all_c[which(rownames(table_all_c) == paste(lines_all_c[ii], 'REP_C', sep = '_')),
                       which(colnames(table_all_c) == paste(lines_all_c[ii], 'REP_A', sep = '_'))]
    }
    else{
      AC = NA
    }
  }
  else{
    AB = NA
    AC = NA
  }
  if(any(colnames(table_all_c) == paste(lines_all_c[ii], 'REP_B', sep = '_'))){
    if(any(rownames(table_all_c) == paste(lines_all_c[ii], 'REP_C', sep = '_'))){
      BC = table_all_c[which(rownames(table_all_c) == paste(lines_all_c[ii], 'REP_C', sep = '_')),
                       which(colnames(table_all_c) == paste(lines_all_c[ii], 'REP_B', sep = '_'))]
    }
    else{
      BC = NA
    }
  }
  else{
    BC = NA
  }
  line_correlation_all_c = rbind(line_correlation_all_c, cbind.data.frame(AB, AC, BC))
  rownames(line_correlation_all_c)[nrow(line_correlation_all_c)] = lines_all_c[ii]
}
colnames(line_correlation_all_c) = c('cor(AB)', 'cor(AC)', 'cor(BC)')
## All genes, uncorrected LFCs
lines_all_uc = unique(sub('_REP_.', '', colnames(table_all_uc)))
line_correlation_all_uc = data.frame(Cell_line = c(), Cor_AB = c(), Cor_AC = c(), Cor_BC = c())
for (ii in 1:length(lines_all_uc)) {
  if(any(colnames(table_all_uc) == paste(lines_all_uc[ii], 'REP_A', sep = '_'))){
    if(any(rownames(table_all_uc) == paste(lines_all_uc[ii], 'REP_B', sep = '_'))){
      AB = table_all_uc[which(rownames(table_all_uc) == paste(lines_all_uc[ii], 'REP_B', sep = '_')),
                        which(colnames(table_all_uc) == paste(lines_all_uc[ii], 'REP_A', sep = '_'))]
    }
    else{
      AB = NA
    }
    if(any(rownames(table_all_uc) == paste(lines_all_uc[ii], 'REP_C', sep = '_'))){
      AC = table_all_uc[which(rownames(table_all_uc) == paste(lines_all_uc[ii], 'REP_C', sep = '_')),
                        which(colnames(table_all_uc) == paste(lines_all_uc[ii], 'REP_A', sep = '_'))]
    }
    else{
      AC = NA
    }
  }
  else{
    AB = NA
    AC = NA
  }
  if(any(colnames(table_all_uc) == paste(lines_all_uc[ii], 'REP_B', sep = '_'))){
    if(any(rownames(table_all_uc) == paste(lines_all_uc[ii], 'REP_C', sep = '_'))){
      BC = table_all_uc[which(rownames(table_all_uc) == paste(lines_all_uc[ii], 'REP_C', sep = '_')),
                        which(colnames(table_all_uc) == paste(lines_all_uc[ii], 'REP_B', sep = '_'))]
    }
    else{
      BC = NA
    }
  }
  else{
    BC = NA
  }
  line_correlation_all_uc = rbind(line_correlation_all_uc, cbind.data.frame(AB, AC, BC))
  rownames(line_correlation_all_uc)[nrow(line_correlation_all_uc)] = lines_all_uc[ii]
}
colnames(line_correlation_all_uc) = c('cor(AB)', 'cor(AC)', 'cor(BC)')

## Distribution of all correlations between replicates of different cell lines
## corrected, highly variable genes
between_correlation_c = table_c
for (ii in 1:nrow(table_c)) {
  for (jj in 1:ncol(table_c)) {
    if (sub('_REP_.', '', rownames(between_correlation_c)[ii]) ==
        sub('_REP_.', '', colnames(between_correlation_c)[jj])){
      between_correlation_c[ii,jj] = NA
    }
  }
}
pdf('Correlation_Coefficients_Replicates_Within_Between_Cell_Lines_After_Correction_Highly_Variable_Genes.pdf')
plot(density(c(line_correlation_c$`cor(AB)`,line_correlation_c$`cor(AC)`,line_correlation_c$`cor(BC)`),na.rm = TRUE),
     xlim = c(0.1, 1.05), col = 'transparent', xlab = 'Correlation Coefficient', main = '', bty = 'n')
polygon(density(between_correlation_c, na.rm = TRUE), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(c(line_correlation_c$`cor(AB)`,line_correlation_c$`cor(AC)`,line_correlation_c$`cor(BC)`),
                na.rm = TRUE), col = rgb (0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(c(line_correlation_c$`cor(AB)`,line_correlation_c$`cor(AC)`,line_correlation_c$`cor(BC)`), na.rm = TRUE),
       col = rgb (0.2117647, 0.3921569, 0.5450980), lty = 2, lw = 3)
abline(v = median(between_correlation_c, na.rm = TRUE), col = rgb(0.8, 0, 0), lty = 2, lwd = 3)
legend(0.08, 5, c('Within Cell Lines', 'Between Cell Lines'), lty = 1, lwd = 5, box.lty = 0,
       col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.8, 0, 0)), alpha = 0.6))
dev.off()

## uncorrected, highly variable genes
between_correlation_uc = table_uc
for (ii in 1:nrow(table_uc)) {
  for (jj in 1:ncol(table_uc)) {
    if (sub('_REP_.', '', rownames(between_correlation_uc)[ii]) ==
        sub('_REP_.', '', colnames(between_correlation_uc)[jj])){
      between_correlation_uc[ii,jj] = NA
    }
  }
}
pdf('Correlation_Coefficients_Replicates_Within_Between_Cell_Lines_Before_Correction_Highly_Variable_Genes.pdf')
plot(density(c(line_correlation_uc$`cor(AB)`,line_correlation_uc$`cor(AC)`,line_correlation_uc$`cor(BC)`),na.rm = TRUE),
     xlim = c(0.1, 1.05), col = 'transparent', xlab = 'Correlation Coefficient', main = '', bty = 'n')
polygon(density(between_correlation_uc, na.rm = TRUE), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(c(line_correlation_uc$`cor(AB)`,line_correlation_uc$`cor(AC)`,line_correlation_uc$`cor(BC)`),
                na.rm = TRUE), col = rgb (0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(c(line_correlation_uc$`cor(AB)`,line_correlation_uc$`cor(AC)`,line_correlation_uc$`cor(BC)`),na.rm = TRUE),
       col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lwd = 3)
abline(v = median(between_correlation_uc, na.rm = TRUE), col = rgb(0.8, 0, 0), lty = 2, lwd = 3)
legend(0.08, 6, c('Within Cell Lines', 'Between Cell Lines'), lty = 1, lwd = 5, box.lty = 0,
       col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.8, 0, 0)), alpha = .6))
dev.off()

## corrected, all genes
between_correlation_all_c = table_all_c
for (ii in 1:nrow(table_all_c)) {
  for (jj in 1:ncol(table_all_c)) {
    if (sub('_REP_.', '', rownames(between_correlation_all_c)[ii]) ==
        sub('_REP_.', '', colnames(between_correlation_all_c)[jj])){
      between_correlation_all_c[ii,jj] = NA
    }
  }
}
pdf('Correlation_Coefficients_Replicates_Within_Between_Cell_Lines_After_Correction_All_Genes.pdf')
plot(density(c(line_correlation_all_c$`cor(AB)`,line_correlation_all_c$`cor(AC)`,line_correlation_all_c$`cor(BC)`),na.rm = TRUE),
     xlim = c(0.1, 1.05), col = 'transparent', xlab = 'Correlation Coefficient', main = '', bty = 'n')
polygon(density(between_correlation_all_c, na.rm = TRUE), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(c(line_correlation_all_c$`cor(AB)`,line_correlation_all_c$`cor(AC)`,line_correlation_all_c$`cor(BC)`),
                na.rm = TRUE), col = rgb (0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(c(line_correlation_all_c$`cor(AB)`,line_correlation_all_c$`cor(AC)`,line_correlation_all_c$`cor(BC)`), na.rm = TRUE),
       col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lwd = 3)
abline(v = median(between_correlation_all_c, na.rm = TRUE), col = rgb(0.8, 0, 0), lty = 2, lwd = 3)
legend(0.08, 5, c('Within Cell Lines', 'Between Cell Lines'), lty = 1, lwd = 5, box.lty = 0,
       col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.8, 0, 0)), alpha = .6))
dev.off()

## Uncorrected, all genes
between_correlation_all_uc = table_all_uc
for (ii in 1:nrow(table_all_uc)) {
  for (jj in 1:ncol(table_all_uc)) {
    if (sub('_REP_.', '', rownames(between_correlation_all_uc)[ii]) ==
        sub('_REP_.', '', colnames(between_correlation_all_uc)[jj])){
      between_correlation_all_uc[ii,jj] = NA
    }
  }
}
pdf('Correlation_Coefficients_Replicates_Within_Between_Cell_Lines_Before_Correction_All_Genes.pdf')
plot(density(c(line_correlation_all_uc$`cor(AB)`,line_correlation_all_uc$`cor(AC)`,line_correlation_all_uc$`cor(BC)`),na.rm = TRUE),
     xlim = c(0.1, 1.05), col = 'transparent', xlab = 'Correlation Coefficient', main = '', bty = 'n')
polygon(density(between_correlation_all_uc, na.rm = TRUE), col = rgb(0.8, 0, 0, alpha = .6), border = 'transparent')
polygon(density(c(line_correlation_all_uc$`cor(AB)`,line_correlation_all_uc$`cor(AC)`,line_correlation_all_uc$`cor(BC)`),
                na.rm = TRUE), col = rgb (0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')
abline(v = median(c(line_correlation_all_uc$`cor(AB)`,line_correlation_all_uc$`cor(AC)`,line_correlation_all_uc$`cor(BC)`),na.rm = TRUE),
       col = rgb(0.2117647, 0.3921569, 0.5450980), lty = 2, lwd = 3)
abline(v = median(between_correlation_all_uc, na.rm = TRUE), col = rgb(0.8, 0, 0), lty = 2, lwd = 3)
legend(0.08, 6, c('Within Cell Lines', 'Between Cell Lines'), lty = 1, lwd = 5, box.lty = 0,
       col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.8, 0, 0)), alpha = .6))
dev.off()

## Sanity checks for correlations
SanityCheck3_1 = isSymmetric(table_c) && isSymmetric(table_uc) && isSymmetric(table_all_c) && isSymmetric(table_all_uc)
SanityCheck3_2 = all(colnames(line_correlation_c) == colnames(line_correlation_uc)) &&
  all(colnames(line_correlation_all_c) == colnames(line_correlation_all_uc)) &&
  all(colnames(line_correlation_c) == colnames(line_correlation_all_c))
SanityCheck3_3 = all(rownames(line_correlation_c) == rownames(line_correlation_uc)) &&
  all(rownames(line_correlation_all_c) == rownames(line_correlation_all_uc)) &&
  all(rownames(line_correlation_c) == rownames(line_correlation_all_c))
SanityCheck3_4 = all(colnames(sorted_runs) == colnames(table_c))

# SanityCheck3_1 # Are all correlation matrices symmetrical?
# SanityCheck3_2 # Do all data frames with intraline correlations have the same colnames?
# SanityCheck3_3 # Do all data frames with intraline correlations have the same cell lines?
# SanityCheck3_4 # Do the names of the condition line up for the data frame with correlations and the one with information on sequencing run?
all(SanityCheck3_1,SanityCheck3_2,SanityCheck3_3,SanityCheck3_4)

write.xlsx(as.data.frame(table_c), 'correlation matrix after correction most variable genes.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(table_uc), 'correlation matrix before correction most variable genes.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(table_all_c), 'correlation matrix after correction all genes.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(table_all_uc), 'correlation matrix before correction all genes.xlsx', rowNames = TRUE)

write.xlsx(line_correlation_c, 'Intra-line correlations top variable genes corrected LFC.xlsx', rowNames = TRUE)
write.xlsx(line_correlation_uc, 'Intra-line correlations top variable genes uncorrected LFC.xlsx', rowNames = TRUE)
write.xlsx(line_correlation_all_c, 'Intra-line correlations all genes corrected LFC.xlsx', rowNames = TRUE)
write.xlsx(line_correlation_all_uc, 'Intra-line correlations all genes uncorrected LFC.xlsx', rowNames = TRUE)

plot_line_correlation_c = cbind.data.frame(apply(line_correlation_c, 1, median, na.rm = TRUE),
                                           apply(line_correlation_c, 1, min, na.rm = TRUE),
                                           apply(line_correlation_c, 1, max, na.rm = TRUE),
                                           apply(line_correlation_c, 1, mean, na.rm = TRUE))
colnames(plot_line_correlation_c) = c('median', 'min', 'max', 'mean')
grand_median = median(plot_line_correlation_c$median)
plot_line_correlation_c = plot_line_correlation_c[order(plot_line_correlation_c$median),]
pdf('Intraline_Correlations.pdf', 12, 10)
par(mar = c(6.1, 4.1, 4.1, 2.1))
barplot(plot_line_correlation_c$median, width = 1, space = 0, col = 'grey', 
        border = rgb(0.7291, 0.6588, 0.4549), las = 2, ylim = c(0,1),
        ylab = "Pearson's Correlation Coefficient", names.arg = rownames(plot_line_correlation_c))
abline(h = grand_median, col = rgb(0.8, 0, 0), lwd = 3)
segments((c(1:nrow(plot_line_correlation_c))-1), plot_line_correlation_c$min, c(1:nrow(plot_line_correlation_c)),
         plot_line_correlation_c$min, col = 'black', lwd = 2)
segments((c(1:nrow(plot_line_correlation_c))-1), plot_line_correlation_c$max, c(1:nrow(plot_line_correlation_c)),
         plot_line_correlation_c$max, col = 'black', lwd = 2)
# points((c(1:nrow(plot_line_correlation_c))-0.5), plot_line_correlation_c$mean, pch = 20)
# abline(h = median(c(line_correlation_c$`cor(AB)`,line_correlation_c$`cor(AC)`,line_correlation_c$`cor(BC)`),
#                   na.rm = TRUE), lwd = 3)
segments((c(1:nrow(plot_line_correlation_c))-0.5), plot_line_correlation_c$max,
         (c(1:nrow(plot_line_correlation_c))-0.5), plot_line_correlation_c$min, col = 'black')
legend(-1.3, 1.005, c('Median', 'Max/Min', 'Grand median'), lty = 1, lwd = 5, box.lty = 0,
       col = rgb(rbind((col2rgb('grey')[,1]/255), c(0, 0, 0), c(0.8, 0, 0))))
dev.off()
# plot(density(between_correlation_c[which(sub('_REP_.', '', rownames(between_correlation_c)) == 'SY5Y'),],
# na.rm = TRUE))
# abline(v = line_correlation_c[which(sub('_REP_.', '', rownames(line_correlation_c)) == 'SY5Y'),],
# col = rgb(0.8, 0, 0))
# null_distribution_wor = c()
# for (ii in 1:1000000) {
# col_per = sample(1:ncol(table_c), 1)
# row_per = sample(1:nrow(table_c), 1)
# if (row_per == col_per) {
# while (row_per == col_per) {
# row_per = sample(1:nrow(table_c), 1)
# }
# }
# null_distribution_wor[ii] = table_c[row_per, col_per]
# }
# null_distribution_wol = c()
# col_per = sample(1:ncol(table_c), 1000000, replace = TRUE)
# row_per = sample(1:nrow(table_c), 1000000, replace = TRUE)
# for (ii in 1:1000000) {
# col_per = sample(1:ncol(table_c), 1)
# row_per = sample(1:nrow(table_c), 1)
# if (sub('_REP_.', '', colnames(table_c)[col_per]) == sub('_REP_.', '', rownames(table_c)[row_per])) {
# while (sub('_REP_.', '', colnames(table_c)[col_per]) == sub('_REP_.', '', rownames(table_c)[row_per])) {
# row_per = sample(1:nrow(table_c), 1)
# }
# }
# null_distribution_wol[ii] = table_c[row_per, col_per]
# }
# plot(density(between_correlation_c, na.rm = TRUE))
# polygon(density(null_distribution_wol), col = rgb(0.8,0,0, alpha = .6), border = 'transparent')
# polygon(density(null_distribution_wor), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = .6), border = 'transparent')

## Umap

## Names of cell lines in alphabetical order
cl = unique(sub('_REP_.', '', colnames(lfc_cor[,-(1:2)])))
cl = cl[order(cl)]

## Umap with all sgRNAs
umap_c = umap(t(lfc_cor[,-(1:2)]))
umap_c_plot = as.data.frame(umap_c$layout)
colnames(umap_c_plot) = c('x', 'y')
umap_c_plot = umap_c_plot[order(rownames(umap_c_plot)),]
umap_c_plot = read.xlsx('umap_c_plot.xlsx', rowNames = TRUE) # to have comparable axis every run


## Define colors
## Cell lines
palette_cl = brewer.pal(8, 'Dark2')
palette_cl = palette_cl[-5]
numbers_cl = cbind(c(rep(1:6, each = 7), rep(7, 3)),c(rep(1:7,6),1:3))
colors_cl = data.frame()
for (ii in 1:nrow(numbers_cl)){
  colors_cl[ii,1] = as.character(palette_cl[numbers_cl[ii,1]])
  colors_cl[ii,2] = as.character(palette_cl[numbers_cl[ii,2]])
}
rownames(colors_cl) = cl
assigned_colors_cl = data.frame()
for (ii in 1:nrow(umap_c_plot)) {
  assigned_colors_cl[ii,(1:2)] = colors_cl[which(rownames(colors_cl) ==
                                                   sub('_REP_.', '', rownames(umap_c_plot)[ii])),]
}
names_plot1 = rownames(colors_cl)[1:23]
col_plot1 = colors_cl[1:23,]
names_plot2 = rownames(colors_cl)[24:nrow(colors_cl)]
col_plot2 = colors_cl[24:nrow(colors_cl),]

## Sequencing run
palette_sr = brewer.pal(4, 'Dark2')
colors_sr = c()
for (ii in 1:ncol(runs)) {
  colors_sr = c(colors_sr, palette_sr[as.numeric(runs[1,ii])])
  names(colors_sr)[ii] = colnames(runs)[ii]
}

## DNA extraction
palette_DNA = brewer.pal(7, 'Dark2')
numbers_DNA = cbind(rep(1:7, each = 2), rep(1:2, 7))
colors_DNA = data.frame()
for (ii in 1:nrow(numbers_DNA)) { # replace numbers by colors
  colors_DNA[ii,1] = palette_DNA[numbers_DNA[ii,1]] # first column with palette
  if (numbers_DNA[ii,2] == 1) { # second colum for circle
    colors_DNA[ii,2] = palette_DNA[numbers_DNA[ii,1]] # 1 same color
  }
  else { # 2 black
    colors_DNA[ii,2] = '#000000'
  }
}
assigned_colors_DNA = data.frame()
for (ii in 1:ncol(runs)) {
  assigned_colors_DNA[ii,(1:2)] = colors_DNA[as.numeric(runs[3,ii]),]
  rownames(assigned_colors_DNA)[ii] = colnames(runs)[ii]
}

## Lineage of cell lines
categories_lineage = unique(lineage$Lineage)
palette_lin = brewer.pal(8, 'Dark2')
numbers_lin = cbind(rep(1:8, 2), rep(c(1,2), each = 8))
colors_lin = data.frame()
for (ii in 1:length(categories_lineage)) { # replace numbers by colors
  colors_lin[ii,1] = palette_lin[numbers_lin[ii,1]]
  if (numbers_lin[ii,2] == 1) {
    colors_lin[ii,2] = palette_lin[numbers_lin[ii,1]]
  }
  else {
    colors_lin[ii,2] = '#000000'
  }
}
rownames(colors_lin) = categories_lineage
lines_colors_lin = data.frame()
for (ii in 1:nrow(lineage)) { # match colors to cell lines
  lines_colors_lin[ii,(1:2)] = colors_lin[which(rownames(colors_lin) == lineage$Lineage[ii]),]
  lines_colors_lin[ii,3] = lineage$Lineage[ii]
  rownames(lines_colors_lin)[ii] = lineage$Cell.line[ii]
}
colnames(lines_colors_lin)[3] = 'Lineage'
assigned_colors_lin = data.frame()
for (ii in 1:nrow(umap_c_plot)) { # match colors to replicates
  assigned_colors_lin[ii,(1:2)] = lines_colors_lin[which(rownames(lines_colors_lin) ==
                                                           sub('_REP_.', '', rownames(umap_c_plot)[ii])),(1:2)]
}
## Change rownames to 'pretty' names for legend
rownames(colors_lin) = str_to_title(rownames(colors_lin))
rownames(colors_lin) = gsub('_', ' ', rownames(colors_lin))
rownames(colors_lin) = gsub('Central nervous system', 'CNS', rownames(colors_lin))
rownames(colors_lin) = gsub('Peripheral nervous system', 'PNS', rownames(colors_lin))
## Change colnames to understandable names to use color coding in other scripts
colnames(colors_lin) = c('inside', 'outside')
## Same for color code of lineages to cell lines
colors_lin_cl = assigned_colors_lin[which(rownames(assigned_colors_lin) %in% lineage$Cell.line),]
colnames(colors_lin_cl) = c('inside', 'outside')
## Save color coding for other use
write.xlsx(colors_lin, 'Color coding lineage.xlsx', rowNames = TRUE)
write.xlsx(colors_lin_cl, 'Color coding lineage for cell lines.xlsx', rowNames = TRUE)

## Make plot for umap with all sgRNAs
## Cell lines
pdf('Umap_of_Replicates_Cell_Lines_Corrected_log_fold_changes.pdf', 13.3, 10)
par(mar = c(5.1, 4.1, 4.1, 21.5))
par(xpd = TRUE)
plot(umap_c_plot, col = assigned_colors_cl$V1[order(rownames(assigned_colors_DNA))],
     bg = assigned_colors_cl$V2[order(rownames(assigned_colors_DNA))], pch = 21,
     lwd = 2.5, cex = 1.7)
legend(1.55, 4.1, names_plot1, pch = 21, col = col_plot1$V1, pt.bg = col_plot1$V2,
       pt.lwd =2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
legend(2.2, 4.1, names_plot2, pch = 21, col = col_plot2$V1, pt.bg = col_plot2$V2, pt.lwd = 2.5,
       pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Sequencing Runs
pdf('Umap_of_Replicates_Sequencing_Runs_Corrected_log_fold_changes.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_c_plot, col = colors_sr, pch = 19,
     lwd = 2.5, cex = 1.7)
legend(1.55, -2.8, c('Sequencing run 1', 'Sequencing run 2', 'Sequencing run 3', 'Sequencing run 4'),
       pch = 19, col = palette_sr, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## DNA extraction
pdf('Umap_of_Replicates_DNA_Extraction_Corrected_log_fold_changes.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_c_plot, col = assigned_colors_DNA$V2, bg = assigned_colors_DNA$V1, pch = 21,
     lwd = 2.5, cex = 1.7)
legend(1.55, 0.2, c('DNA extraction #1', 'DNA extraction #2', 'DNA extraction #3', 'DNA extraction #4',
                    'DNA extraction #5', 'DNA extraction #6', 'DNA extraction #7', 'DNA extraction #8',
                    'DNA extraction #9', 'DNA extraction #10', 'DNA extraction #11', 'DNA extraction #12',
                    'DNA extraction #13', 'DNA extraction #14'),
       col = colors_DNA$V2, pt.bg = colors_DNA$V1, pch = 21, pt.cex = 1.7, pt.lwd = 2.5, bty = 'n', cex = 1.5)
dev.off()

## Lineage
pdf('Umap_of_Replicates_Lineage_Corrected_log_fold_changes.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_c_plot, col = assigned_colors_lin$V2, bg = assigned_colors_lin$V1, pch = 21,
     lwd = 2.5, cex = 1.7)
legend(1.55, 0.8, rownames(colors_lin), pch = 21, col = colors_lin$outside, pt.bg = colors_lin$inside,
       pt.lwd = 2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Lineage subtypes
umap_c_plot_sub = data.frame()
for (ii in 1:nrow(umap_c_plot)) {
  umap_c_plot_sub[ii,(1:2)] = umap_c_plot[ii,]
  umap_c_plot_sub[ii,3] = str_replace(lineage$Lineage.subtype[
    which(lineage$Cell.line == sub('_REP_.', '', rownames(umap_c_plot)[ii]))], '^\\w{1}', toupper)
  rownames(umap_c_plot_sub)[ii] = rownames(umap_c_plot)[ii]
}
colnames(umap_c_plot_sub)[3] = 'Lineage.subtype'
umap_c_plot_sub$Lineage.subtype = gsub('_', ' ', umap_c_plot_sub$Lineage.subtype)

## ATRT
umap_c_plot_sub_nATRT = subset(umap_c_plot_sub, Lineage.subtype != 'ATRT')
umap_c_plot_sub_nATRT = umap_c_plot_sub_nATRT[,-3]
umap_c_plot_sub_ATRT = subset(umap_c_plot_sub, Lineage.subtype == 'ATRT')
umap_c_plot_sub_ATRT = umap_c_plot_sub_ATRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_ATRT_Corrected_log_fold_changes.pdf')
plot(umap_c_plot, col = 'transparent', main = 'ATRT')
points(umap_c_plot_sub_nATRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Colorectal adenocarcinoma
umap_c_plot_sub_nCA = subset(umap_c_plot_sub, Lineage.subtype != 'Colorectal adenocarcinoma')
umap_c_plot_sub_nCA = umap_c_plot_sub_nCA[,-3]
umap_c_plot_sub_CA = subset(umap_c_plot_sub, Lineage.subtype == 'Colorectal adenocarcinoma')
umap_c_plot_sub_CA = umap_c_plot_sub_CA[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Colorectal_adenocarcinoma_Corrected_log_fold_changes.pdf')
plot(umap_c_plot, col = 'transparent', main = 'Colorectal Adenocarcinoma')
points(umap_c_plot_sub_nCA, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_CA, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Glioma
umap_c_plot_sub_nGlioma = subset(umap_c_plot_sub, Lineage.subtype != 'Glioma')
umap_c_plot_sub_nGlioma = umap_c_plot_sub_nGlioma[,-3]
umap_c_plot_sub_Glioma = subset(umap_c_plot_sub, Lineage.subtype == 'Glioma')
umap_c_plot_sub_Glioma = umap_c_plot_sub_Glioma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Glioma_Corrected_log_fold_changes.pdf')
plot(umap_c_plot, col = 'transparent', main = 'Glioma')
points(umap_c_plot_sub_nGlioma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_Glioma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Malignant rhabdoid tumor
umap_c_plot_sub_nMRT = subset(umap_c_plot_sub, Lineage.subtype != 'Malignant rhabdoid tumor')
umap_c_plot_sub_nMRT = umap_c_plot_sub_nMRT[,-3]
umap_c_plot_sub_MRT = subset(umap_c_plot_sub, Lineage.subtype == 'Malignant rhabdoid tumor')
umap_c_plot_sub_MRT = umap_c_plot_sub_MRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Malignant_rhabdoid_tumor_Corrected_log_fold_changes.pdf')
plot(umap_c_plot, col = 'transparent', main = 'Malignant rhabdoid tumor')
points(umap_c_plot_sub_nMRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_MRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Medulloblastoma
umap_c_plot_sub_nMedulloblastoma = subset(umap_c_plot_sub, Lineage.subtype != 'Medulloblastoma')
umap_c_plot_sub_nMedulloblastoma = umap_c_plot_sub_nMedulloblastoma[,-3]
umap_c_plot_sub_Medulloblastoma = subset(umap_c_plot_sub, Lineage.subtype == 'Medulloblastoma')
umap_c_plot_sub_Medulloblastoma = umap_c_plot_sub_Medulloblastoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Medulloblastoma_Corrected_log_fold_changes.pdf')
plot(umap_c_plot, col = 'transparent', main = 'Medulloblastoma')
points(umap_c_plot_sub_nMedulloblastoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_Medulloblastoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Melanoma
umap_c_plot_sub_nMelanoma = subset(umap_c_plot_sub, Lineage.subtype != 'Melanoma')
umap_c_plot_sub_nMelanoma = umap_c_plot_sub_nMelanoma[,-3]
umap_c_plot_sub_Melanoma = subset(umap_c_plot_sub, Lineage.subtype == 'Melanoma')
umap_c_plot_sub_Melanoma = umap_c_plot_sub_Melanoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Melanoma_Corrected_log_fold_changes.pdf')
plot(umap_c_plot, col = 'transparent', main = 'Melanoma')
points(umap_c_plot_sub_nMelanoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_Melanoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## CNS
CNS = subset(lineage, Lineage == 'central_nervous_system')
CNS$Lineage.subtype = str_replace(CNS$Lineage.subtype, '\\w{1}', toupper)
umap_c_plot_sub_nCNS = subset(umap_c_plot_sub, !(Lineage.subtype %in% CNS$Lineage.subtype))
umap_c_plot_sub_Meningioma = subset(umap_c_plot_sub, Lineage.subtype == 'Meningioma')
pdf('Umap_of_Replicates_Lineage_CNS_Subtypes_Corrected_log_fold_change.pdf')
plot(umap_c_plot, col = 'transparent', main = 'CNS')
points(umap_c_plot_sub_nCNS, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_c_plot_sub_Glioma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
points(umap_c_plot_sub_Medulloblastoma, col = rgb(0.7377778, 0.5244444, 0.1288889), pch = 19, cex = 1.7)
points(umap_c_plot_sub_Meningioma, col = rgb(0.3, 0.4, 0.35), pch = 19, cex = 1.7)
legend('bottomright', c('ATRT', 'Glioma', 'Medulloblastoma', 'Meningioma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980), c(0.7377778, 0.5244444, 0.1288889),
                       c(0.3, 0.4, 0.35))), pch = 19, pt.cex = 1.7)
dev.off()

## Embryonal brain tumors
umap_c_plot_sub_nEBT = subset(umap_c_plot_sub, Lineage.subtype != 'ATRT' | Lineage.subtype != 'Medulloblastoma')
umap_c_plot_sub_nEBT = umap_c_plot_sub_nEBT[,-3]
pdf('Umap_of_Replicates_Lineage_EBT_Subtypes_Corrected_log_fold_change.pdf')
plot(umap_c_plot, col = 'transparent', main = 'Embryonal Brain Tumors')
points(umap_c_plot_sub_nEBT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_c_plot_sub_Medulloblastoma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('bottomright', c('ATRT', 'Medulloblastoma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

## Hematopoietic cancers
hem = subset(lineage, Lineage == 'blood' | Lineage == 'lymphocyte')
hem$Lineage.subtype = str_replace(hem$Lineage.subtype, '\\w{1}', toupper)
umap_c_plot_sub_nhem = subset(umap_c_plot_sub, !(Lineage.subtype %in% hem$Lineage.subtype))
umap_c_plot_sub_blood = subset(umap_c_plot_sub, Lineage.subtype == 'AML' | Lineage.subtype == 'CLL')
umap_c_plot_sub_lymphocyte = subset(umap_c_plot_sub, Lineage.subtype == 'Hodgkin lymphoma' |
                                      Lineage.subtype == 'Non hodgkin lymphoma')
pdf('Umap_of_Replicates_Lineage_Hematopoietic_Corrected_log_fold_change.pdf')
plot(umap_c_plot, col = 'transparent', main = 'Hematopoietic Cancers')
points(umap_c_plot_sub_nhem, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_c_plot_sub_blood, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_c_plot_sub_lymphocyte, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('bottomright', c('Blood', 'Lymphocyte'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

## 500 most variable sgRNAs
var_c = cbind(lfc_cor, variance = apply(lfc_cor[,-(1:2)], 1, var))
ordered_var_c = var_c[order(var_c$variance, decreasing = TRUE),]
high_var_c = ordered_var_c[1:500,]

## Umap of most variable sgRNAs
high_var_c_t = t(high_var_c[,-c(1,2,ncol(high_var_c))])
umap_hv_c = umap(high_var_c_t)
umap_hv_c_plot = as.data.frame(umap_hv_c$layout)
colnames(umap_hv_c_plot) = c('x', 'y')
umap_hv_c_plot = umap_hv_c_plot[order(rownames(umap_hv_c_plot)),]
umap_hv_c_plot = read.xlsx('umap_hv_c_plot.xlsx', rowNames = TRUE) # to have comparable axis


## Make plot
## Cell lines
pdf('Umap_of_Replicates_Cell_Lines_Highly_Variable_Guides_Corrected_log_fold_changes.pdf', 13.3, 10)
par(mar = c(5.1, 4.1, 4.1, 21.5))
par(xpd = TRUE)
plot(umap_hv_c_plot, col = assigned_colors_cl$V1, bg = assigned_colors_cl$V2, pch = 21,
     lwd = 2, cex = 1.7)
legend(5, 3.15, names_plot1, pch = 21, col = col_plot1$V1, pt.bg = col_plot1$V2,
       pt.lwd = 2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
legend(7, 3.15, names_plot2, pch = 21, col = col_plot2$V1, pt.bg = col_plot2$V2, pt.lwd = 2.5,
       pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Sequencing Runs
pdf('Umap_of_Replicates_Sequencing_Runs_Highly_Variable_Guides_Corrected_log_fold_changes.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_hv_c_plot, col = colors_sr, pch = 19, lwd = 2, cex = 1.7)
legend(5, -1.05, c('Sequencing run 1', 'Sequencing run 2', 'Sequencing run 3', 'Sequencing run 4'),
       pch = 19, col = palette_sr, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## DNA extraction
pdf('Umap_of_Replicates_DNA_Extraction_Highly_Variable_Guides_Corrected_log_fold_changes.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_hv_c_plot, col = assigned_colors_DNA$V2[order(rownames(assigned_colors_DNA))],
     bg = assigned_colors_DNA$V1[order(rownames(assigned_colors_DNA))], pch = 21,
     lwd = 2.5, cex = 1.7)
legend(5, 0.42, c('DNA extraction #1', 'DNA extraction #2', 'DNA extraction #3', 'DNA extraction #4',
                  'DNA extraction #5', 'DNA extraction #6', 'DNA extraction #7', 'DNA extraction #8',
                  'DNA extraction #9', 'DNA extraction #10', 'DNA extraction #11', 'DNA extraction #12',
                  'DNA extraction #13', 'DNA extraction #14'),
       col = colors_DNA$V2, pt.bg = colors_DNA$V1, pch = 21, pt.cex = 1.7, pt.lwd = 2.5, bty = 'n', cex = 1.5)
dev.off()


## Lineage
pdf('Umap_of_Replicates_Lineage_Highly_Variable_Guides_Corrected_log_fold_changes.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_hv_c_plot, col = assigned_colors_lin$V2, bg = assigned_colors_lin$V1, pch = 21,
     lwd = 2, cex = 1.7)
legend(5, 1.15, rownames(colors_lin), pch = 21, col = colors_lin$V2, pt.bg = colors_lin$V1,
       pt.lwd = 2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Lineage subtypes
umap_hv_c_plot_sub = data.frame()
for (ii in 1:nrow(umap_hv_c_plot)) {
  umap_hv_c_plot_sub[ii,(1:2)] = umap_hv_c_plot[ii,]
  umap_hv_c_plot_sub[ii,3] = str_replace(lineage$Lineage.subtype[
    which(lineage$Cell.line == sub('_REP_.', '', rownames(umap_hv_c_plot)[ii]))], '^\\w{1}', toupper)
  rownames(umap_hv_c_plot_sub)[ii] = rownames(umap_hv_c_plot)[ii]
}
colnames(umap_hv_c_plot_sub)[3] = 'Lineage.subtype'
umap_hv_c_plot_sub$Lineage.subtype = gsub('_', ' ', umap_hv_c_plot_sub$Lineage.subtype)

## ATRT
umap_hv_c_plot_sub_nATRT = subset(umap_hv_c_plot_sub, Lineage.subtype != 'ATRT')
umap_hv_c_plot_sub_nATRT = umap_hv_c_plot_sub_nATRT[,-3]
umap_hv_c_plot_sub_ATRT = subset(umap_hv_c_plot_sub, Lineage.subtype == 'ATRT')
umap_hv_c_plot_sub_ATRT = umap_hv_c_plot_sub_ATRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_ATRT_Highly_Variable_Guides_Corrected_log_fold_changes.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'ATRT')
points(umap_hv_c_plot_sub_nATRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Colorectal adenocarcinoma
umap_hv_c_plot_sub_nCA = subset(umap_hv_c_plot_sub, Lineage.subtype != 'Colorectal adenocarcinoma')
umap_hv_c_plot_sub_nCA = umap_hv_c_plot_sub_nCA[,-3]
umap_hv_c_plot_sub_CA = subset(umap_hv_c_plot_sub, Lineage.subtype == 'Colorectal adenocarcinoma')
umap_hv_c_plot_sub_CA = umap_hv_c_plot_sub_CA[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Colorectal_adenocarcinoma_Highly_Variable_Guides_Corrected_log_fold_changes.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'Colorectal adenocarcinoma')
points(umap_hv_c_plot_sub_nCA, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_CA, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Glioma
umap_hv_c_plot_sub_nGlioma = subset(umap_hv_c_plot_sub, Lineage.subtype != 'Glioma')
umap_hv_c_plot_sub_nGlioma = umap_hv_c_plot_sub_nGlioma[,-3]
umap_hv_c_plot_sub_Glioma = subset(umap_hv_c_plot_sub, Lineage.subtype == 'Glioma')
umap_hv_c_plot_sub_Glioma = umap_hv_c_plot_sub_Glioma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Glioma_Highly_Variable_Guides_Corrected_log_fold_changes.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'Glioma')
points(umap_hv_c_plot_sub_nGlioma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_Glioma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Malignant rhabdoid tumor
umap_hv_c_plot_sub_nMRT = subset(umap_hv_c_plot_sub, Lineage.subtype != 'Malignant rhabdoid tumor')
umap_hv_c_plot_sub_nMRT = umap_hv_c_plot_sub_nMRT[,-3]
umap_hv_c_plot_sub_MRT = subset(umap_hv_c_plot_sub, Lineage.subtype == 'Malignant rhabdoid tumor')
umap_hv_c_plot_sub_MRT = umap_hv_c_plot_sub_MRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Malignant_rhabdoid_tumor_Highly_Variable_Guides_Corrected_log_fold_changes.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'Malignant rhabdoid tumor')
points(umap_hv_c_plot_sub_nMRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_MRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Medulloblastoma
umap_hv_c_plot_sub_nMedulloblastoma = subset(umap_hv_c_plot_sub, Lineage.subtype != 'Medulloblastoma')
umap_hv_c_plot_sub_nMedulloblastoma = umap_hv_c_plot_sub_nMedulloblastoma[,-3]
umap_hv_c_plot_sub_Medulloblastoma = subset(umap_hv_c_plot_sub, Lineage.subtype == 'Medulloblastoma')
umap_hv_c_plot_sub_Medulloblastoma = umap_hv_c_plot_sub_Medulloblastoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Medulloblastoma_Highly_Variable_Guides_Corrected_log_fold_changes.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'Medulloblastoma')
points(umap_hv_c_plot_sub_nMedulloblastoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_Medulloblastoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Melanoma
umap_hv_c_plot_sub_nMelanoma = subset(umap_hv_c_plot_sub, Lineage.subtype != 'Melanoma')
umap_hv_c_plot_sub_nMelanoma = umap_hv_c_plot_sub_nMelanoma[,-3]
umap_hv_c_plot_sub_Melanoma = subset(umap_hv_c_plot_sub, Lineage.subtype == 'Melanoma')
umap_hv_c_plot_sub_Melanoma = umap_hv_c_plot_sub_Melanoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Melanoma_Highly_Variable_Guides_Corrected_log_fold_changes.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'Melanoma')
points(umap_hv_c_plot_sub_nMelanoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_Melanoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## CNS
umap_hv_c_plot_sub_nCNS = subset(umap_hv_c_plot_sub, !(Lineage.subtype %in% CNS$Lineage.subtype))
umap_hv_c_plot_sub_Meningioma = subset(umap_hv_c_plot_sub, Lineage.subtype == 'Meningioma')
pdf('Umap_of_Replicates_Lineage_CNS_Subtypes_Highly_Variable_Guides_Corrected_log_fold_change.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'CNS')
points(umap_hv_c_plot_sub_nCNS, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_Glioma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_Medulloblastoma, col = rgb(0.7377778, 0.5244444, 0.1288889), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_Meningioma, col = rgb(0.3, 0.4, 0.35), pch = 19, cex = 1.7)
legend('topright', c('ATRT', 'Glioma', 'Medulloblastoma', 'Meningioma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980), c(0.7377778, 0.5244444, 0.1288889),
                       c(0.3, 0.4, 0.35))), pch = 19, pt.cex = 1.7)
dev.off()

## Embryonal brain tumors
umap_hv_c_plot_sub_nEBT = subset(umap_hv_c_plot_sub, Lineage.subtype != 'ATRT' | Lineage.subtype != 'Medulloblastoma')
umap_hv_c_plot_sub_nEBT = umap_hv_c_plot_sub_nEBT[,-3]
pdf('Umap_of_Replicates_Lineage_EBT_Subtypes_Highly_Variable_Guides_Corrected_log_fold_change.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'Embryonal Brain Tumors')
points(umap_hv_c_plot_sub_nEBT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_Medulloblastoma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('topright', c('ATRT', 'Medulloblastoma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

## Hematopoietic cancers
umap_hv_c_plot_sub_nhem = subset(umap_hv_c_plot_sub, !(Lineage.subtype %in% hem$Lineage.subtype))
umap_hv_c_plot_sub_blood = subset(umap_hv_c_plot_sub, Lineage.subtype == 'AML' | Lineage.subtype == 'CLL')
umap_hv_c_plot_sub_lymphocyte = subset(umap_hv_c_plot_sub, Lineage.subtype == 'Hodgkin lymphoma' |
                                         Lineage.subtype == 'Non hodgkin lymphoma')
pdf('Umap_of_Replicates_Lineage_Hematopoietic_Highly_Variable_Guides_Corrected_log_fold_change.pdf')
plot(umap_hv_c_plot, col = 'transparent', main = 'Hematopoietic Cancers')
points(umap_hv_c_plot_sub_nhem, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_blood, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_hv_c_plot_sub_lymphocyte, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('topright', c('Blood', 'Lymphocyte'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

## Umaps using corrected read counts
umap_rcc = umap(t(rc_cor[,-c((1:2), which(colnames(rc_cor) == 'spike_in_plasmid_DNA_from_CP1608'))]))
umap_rcc_plot = as.data.frame(umap_rcc$layout)
colnames(umap_rcc_plot) = c('x', 'y')
umap_rcc_plot = umap_rcc_plot[order(rownames(umap_rcc_plot)),]
umap_rcc_plot = read.xlsx('umap_rcc_plot.xlsx', rowNames = TRUE)

pdf('Umap_of_Replicates_Cell_Lines_Corrected_read_counts.pdf', 13.3, 10)
par(mar = c(5.1, 4.1, 4.1, 21.5))
par(xpd = TRUE)
plot(umap_rcc_plot, col = assigned_colors_cl$V1, bg = assigned_colors_cl$V2, pch = 21,
     lwd = 2.5, cex = 1.7)
legend(2.2, 3.75, names_plot1, pch = 21, col = col_plot1$V1, pt.bg = col_plot1$V2,
       pt.lwd =2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
legend(3.4, 3.75, names_plot2, pch = 21, col = col_plot2$V1, pt.bg = col_plot2$V2, pt.lwd = 2.5,
       pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Sequencing Runs
pdf('Umap_of_Replicates_Sequencing_Runs_Corrected_read_counts.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_rcc_plot, col = colors_sr, pch = 19,
     lwd = 2.5, cex = 1.7)
legend(2.2, -2.6, c('Sequencing run 1', 'Sequencing run 2', 'Sequencing run 3', 'Sequencing run 4'),
       pch = 19, col = palette_sr, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## DNA extraction
pdf('Umap_of_Replicates_DNA_Extraction_Corrected_read_counts.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_rcc_plot, col = assigned_colors_DNA$V2, bg = assigned_colors_DNA$V1, pch = 21,
     lwd = 2.5, cex = 1.7)
legend(2.2, 0.23, c('DNA extraction #1', 'DNA extraction #2', 'DNA extraction #3', 'DNA extraction #4',
                    'DNA extraction #5', 'DNA extraction #6', 'DNA extraction #7', 'DNA extraction #8',
                    'DNA extraction #9', 'DNA extraction #10', 'DNA extraction #11', 'DNA extraction #12',
                    'DNA extraction #13', 'DNA extraction #14'),
       col = colors_DNA$V2, pt.bg = colors_DNA$V1, pch = 21, pt.cex = 1.7, pt.lwd = 2.5, bty = 'n', cex = 1.5)
dev.off()

## Lineage
pdf('Umap_of_Replicates_Lineage_Corrected_read_counts.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_rcc_plot, col = assigned_colors_lin$V2, bg = assigned_colors_lin$V1, pch = 21,
     lwd = 2.5, cex = 1.7)
legend(2.2, 0.72, rownames(colors_lin), pch = 21, col = colors_lin$V2, pt.bg = colors_lin$V1,
       pt.lwd = 2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Lineage subtypes
umap_rcc_plot_sub = data.frame()
for (ii in 1:nrow(umap_rcc_plot)) {
  umap_rcc_plot_sub[ii,(1:2)] = umap_rcc_plot[ii,]
  umap_rcc_plot_sub[ii,3] = str_replace(lineage$Lineage.subtype[
    which(lineage$Cell.line == sub('_REP_.', '', rownames(umap_rcc_plot)[ii]))], '^\\w{1}', toupper)
  rownames(umap_rcc_plot_sub)[ii] = rownames(umap_rcc_plot)[ii]
}
colnames(umap_rcc_plot_sub)[3] = 'Lineage.subtype'
umap_rcc_plot_sub$Lineage.subtype = gsub('_', ' ', umap_rcc_plot_sub$Lineage.subtype)

## ATRT
umap_rcc_plot_sub_nATRT = subset(umap_rcc_plot_sub, Lineage.subtype != 'ATRT')
umap_rcc_plot_sub_nATRT = umap_rcc_plot_sub_nATRT[,-3]
umap_rcc_plot_sub_ATRT = subset(umap_rcc_plot_sub, Lineage.subtype == 'ATRT')
umap_rcc_plot_sub_ATRT = umap_rcc_plot_sub_ATRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_ATRT_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'ATRT')
points(umap_rcc_plot_sub_nATRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Colorectal adenocarcinoma
umap_rcc_plot_sub_nCA = subset(umap_rcc_plot_sub, Lineage.subtype != 'Colorectal adenocarcinoma')
umap_rcc_plot_sub_nCA = umap_rcc_plot_sub_nCA[,-3]
umap_rcc_plot_sub_CA = subset(umap_rcc_plot_sub, Lineage.subtype == 'Colorectal adenocarcinoma')
umap_rcc_plot_sub_CA = umap_rcc_plot_sub_CA[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Colorectal_adenocarcinoma_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'CA')
points(umap_rcc_plot_sub_nCA, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_CA, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Glioma
umap_rcc_plot_sub_nGlioma = subset(umap_rcc_plot_sub, Lineage.subtype != 'Glioma')
umap_rcc_plot_sub_nGlioma = umap_rcc_plot_sub_nGlioma[,-3]
umap_rcc_plot_sub_Glioma = subset(umap_rcc_plot_sub, Lineage.subtype == 'Glioma')
umap_rcc_plot_sub_Glioma = umap_rcc_plot_sub_Glioma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Glioma_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'Glioma')
points(umap_rcc_plot_sub_nGlioma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_Glioma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Malignant rhabdoid tumor
umap_rcc_plot_sub_nMRT = subset(umap_rcc_plot_sub, Lineage.subtype != 'Malignant rhabdoid tumor')
umap_rcc_plot_sub_nMRT = umap_rcc_plot_sub_nMRT[,-3]
umap_rcc_plot_sub_MRT = subset(umap_rcc_plot_sub, Lineage.subtype == 'Malignant rhabdoid tumor')
umap_rcc_plot_sub_MRT = umap_rcc_plot_sub_MRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Malignant_rhabdoid_tumor_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'Malignant rhabdoid tumor')
points(umap_rcc_plot_sub_nMRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_MRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Medulloblastoma
umap_rcc_plot_sub_nMedulloblastoma = subset(umap_rcc_plot_sub, Lineage.subtype != 'Medulloblastoma')
umap_rcc_plot_sub_nMedulloblastoma = umap_rcc_plot_sub_nMedulloblastoma[,-3]
umap_rcc_plot_sub_Medulloblastoma = subset(umap_rcc_plot_sub, Lineage.subtype == 'Medulloblastoma')
umap_rcc_plot_sub_Medulloblastoma = umap_rcc_plot_sub_Medulloblastoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Medulloblastoma_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'Medulloblastoma')
points(umap_rcc_plot_sub_nMedulloblastoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_Medulloblastoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Melanoma
umap_rcc_plot_sub_nMelanoma = subset(umap_rcc_plot_sub, Lineage.subtype != 'Melanoma')
umap_rcc_plot_sub_nMelanoma = umap_rcc_plot_sub_nMelanoma[,-3]
umap_rcc_plot_sub_Melanoma = subset(umap_rcc_plot_sub, Lineage.subtype == 'Melanoma')
umap_rcc_plot_sub_Melanoma = umap_rcc_plot_sub_Melanoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Melanoma_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'Melanoma')
points(umap_rcc_plot_sub_nMelanoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_Melanoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## CNS
CNS = subset(lineage, Lineage == 'central_nervous_system')
CNS$Lineage.subtype = str_replace(CNS$Lineage.subtype, '\\w{1}', toupper)
umap_rcc_plot_sub_nCNS = subset(umap_rcc_plot_sub, !(Lineage.subtype %in% CNS$Lineage.subtype))
umap_rcc_plot_sub_Meningioma = subset(umap_rcc_plot_sub, Lineage.subtype == 'Meningioma')
pdf('Umap_of_Replicates_Lineage_CNS_Subtypes_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'CNS')
points(umap_rcc_plot_sub_nCNS, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_Glioma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_Medulloblastoma, col = rgb(0.7377778, 0.5244444, 0.1288889), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_Meningioma, col = rgb(0.3, 0.4, 0.35), pch = 19, cex = 1.7)
legend('bottomright', c('ATRT', 'Glioma', 'Medulloblastoma', 'Meningioma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980), c(0.7377778, 0.5244444, 0.1288889),
                       c(0.3, 0.4, 0.35))), pch = 19, pt.cex = 1.7)
dev.off()

## Embryonal brain tumors
umap_rcc_plot_sub_nEBT = subset(umap_rcc_plot_sub, Lineage.subtype != 'ATRT' | Lineage.subtype != 'Medulloblastoma')
umap_rcc_plot_sub_nEBT = umap_rcc_plot_sub_nEBT[,-3]
pdf('Umap_of_Replicates_Lineage_EBT_Subtypes_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'Embryonal Brain Tumors')
points(umap_rcc_plot_sub_nEBT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_Medulloblastoma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('bottomright', c('ATRT', 'Medulloblastoma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

## Hematopoietic cancers
umap_rcc_plot_sub_nhem = subset(umap_rcc_plot_sub, !(Lineage.subtype %in% hem$Lineage.subtype))
umap_rcc_plot_sub_blood = subset(umap_rcc_plot_sub, Lineage.subtype == 'AML' | Lineage.subtype == 'CLL')
umap_rcc_plot_sub_lymphocyte = subset(umap_rcc_plot_sub, Lineage.subtype == 'Hodgkin lymphoma' |
                                        Lineage.subtype == 'Non hodgkin lymphoma')
pdf('Umap_of_Replicates_Lineage_Hematopoietic_Corrected_read_counts.pdf')
plot(umap_rcc_plot, col = 'transparent', main = 'Hematopoietic Cancers')
points(umap_rcc_plot_sub_nhem, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_blood, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_rcc_plot_sub_lymphocyte, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('bottomright', c('Blood', 'Lymphocyte'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

## 500 most variable sgRNAs
var_rcc = cbind(rc_cor, variance = apply(rc_cor[,-(1:2)], 1, var))
ordered_var_rcc = var_rcc[order(var_rcc$variance, decreasing = TRUE),]
high_var_rcc = ordered_var_rcc[1:500,]

## Umap of most variable sgRNAs
high_var_rcc_t = t(high_var_rcc[,-c(1,2,ncol(high_var_rcc),
                                    which(colnames(high_var_rcc) == "spike_in_plasmid_DNA_from_CP1608"))])
umap_hv_rcc = umap(high_var_rcc_t)
umap_hv_rcc_plot = as.data.frame(umap_hv_rcc$layout)
colnames(umap_hv_rcc_plot) = c('x', 'y')
umap_hv_rcc_plot = umap_hv_rcc_plot[order(rownames(umap_hv_rcc_plot)),]
umap_hv_rcc_plot = read.xlsx('umap_hv_rcc_plot.xlsx', rowNames = TRUE) # to have comparable axis

## Make plot
## Cell lines
pdf('Umap_of_Replicates_Cell_Lines_Highly_Variable_Guides_Corrected_read_counts.pdf', 13.3, 10)
par(mar = c(5.1, 4.1, 4.1, 21.5))
par(xpd = TRUE)
plot(umap_hv_rcc_plot, col = assigned_colors_cl$V1, bg = assigned_colors_cl$V2, pch = 21,
     lwd = 2, cex = 1.7)
legend(1.95, 6.1, names_plot1, pch = 21, col = col_plot1$V1, pt.bg = col_plot1$V2,
       pt.lwd = 2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
legend(2.85, 6.1, names_plot2, pch = 21, col = col_plot2$V1, pt.bg = col_plot2$V2, pt.lwd = 2.5,
       pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Sequencing Runs
pdf('Umap_of_Replicates_Sequencing_Runs_Highly_Variable_Guides_Corrected_read_counts.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_hv_rcc_plot, col = colors_sr, pch = 19, lwd = 2, cex = 1.7)
legend(1.95, -4.1, c('Sequencing run 1', 'Sequencing run 2', 'Sequencing run 3', 'Sequencing run 4'),
       pch = 19, col = palette_sr, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## DNA extraction
pdf('Umap_of_Replicates_DNA_Extraction_Highly_Variable_Guides_Corrected_read_counts.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_hv_rcc_plot, col = assigned_colors_DNA$V2, bg = assigned_colors_DNA$V1, pch = 21,
     lwd = 2.5, cex = 1.7)
legend(1.95, 0.4, c('DNA extraction #1', 'DNA extraction #2', 'DNA extraction #3', 'DNA extraction #4',
                    'DNA extraction #5', 'DNA extraction #6', 'DNA extraction #7', 'DNA extraction #8',
                    'DNA extraction #9', 'DNA extraction #10', 'DNA extraction #11', 'DNA extraction #12',
                    'DNA extraction #13', 'DNA extraction #14'),
       col = colors_DNA$V2, pt.bg = colors_DNA$V1, pch = 21, pt.cex = 1.7, pt.lwd = 2.5, bty = 'n', cex = 1.5)
dev.off()

## Lineage
pdf('Umap_of_Replicates_Lineage_Highly_Variable_Guides_Corrected_read_counts.pdf', 12.05, 10)
par(mar = c(5.1, 4.1, 4.1, 15.5))
par(xpd = TRUE)
plot(umap_hv_rcc_plot, col = assigned_colors_lin$V2, bg = assigned_colors_lin$V1, pch = 21, lwd = 2, cex = 1.7)
legend(1.95, 1.25, rownames(colors_lin), pch = 21, col = colors_lin$V2, pt.bg = colors_lin$V1,
       pt.lwd = 2.5, pt.cex = 1.7, bty = 'n', cex = 1.5)
dev.off()

## Lineage subtypes
umap_hv_rcc_plot_sub = data.frame()
for (ii in 1:nrow(umap_hv_rcc_plot)) {
  umap_hv_rcc_plot_sub[ii,(1:2)] = umap_hv_rcc_plot[ii,]
  umap_hv_rcc_plot_sub[ii,3] = str_replace(lineage$Lineage.subtype[
    which(lineage$Cell.line == sub('_REP_.', '', rownames(umap_hv_rcc_plot)[ii]))], '^\\w{1}', toupper)
  rownames(umap_hv_rcc_plot_sub)[ii] = rownames(umap_hv_rcc_plot)[ii]
}
colnames(umap_hv_rcc_plot_sub)[3] = 'Lineage.subtype'
umap_hv_rcc_plot_sub$Lineage.subtype = gsub('_', ' ', umap_hv_rcc_plot_sub$Lineage.subtype)

## ATRT
umap_hv_rcc_plot_sub_nATRT = subset(umap_hv_rcc_plot_sub, Lineage.subtype != 'ATRT')
umap_hv_rcc_plot_sub_nATRT = umap_hv_rcc_plot_sub_nATRT[,-3]
umap_hv_rcc_plot_sub_ATRT = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'ATRT')
umap_hv_rcc_plot_sub_ATRT = umap_hv_rcc_plot_sub_ATRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_ATRT_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'ATRT')
points(umap_hv_rcc_plot_sub_nATRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Colorectal adenocarcinoma
umap_hv_rcc_plot_sub_nCA = subset(umap_hv_rcc_plot_sub, Lineage.subtype != 'Colorectal adenocarcinoma')
umap_hv_rcc_plot_sub_nCA = umap_hv_rcc_plot_sub_nCA[,-3]
umap_hv_rcc_plot_sub_CA = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'Colorectal adenocarcinoma')
umap_hv_rcc_plot_sub_CA = umap_hv_rcc_plot_sub_CA[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Colorectal_adenocarcinoma_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'Colorectal Adenocarcinoma')
points(umap_hv_rcc_plot_sub_nCA, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_CA, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Glioma
umap_hv_rcc_plot_sub_nGlioma = subset(umap_hv_rcc_plot_sub, Lineage.subtype != 'Glioma')
umap_hv_rcc_plot_sub_nGlioma = umap_hv_rcc_plot_sub_nGlioma[,-3]
umap_hv_rcc_plot_sub_Glioma = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'Glioma')
umap_hv_rcc_plot_sub_Glioma = umap_hv_rcc_plot_sub_Glioma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Glioma_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'Glioma')
points(umap_hv_rcc_plot_sub_nGlioma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_Glioma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Malignant rhabdoid tumor
umap_hv_rcc_plot_sub_nMRT = subset(umap_hv_rcc_plot_sub, Lineage.subtype != 'Malignant rhabdoid tumor')
umap_hv_rcc_plot_sub_nMRT = umap_hv_rcc_plot_sub_nMRT[,-3]
umap_hv_rcc_plot_sub_MRT = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'Malignant rhabdoid tumor')
umap_hv_rcc_plot_sub_MRT = umap_hv_rcc_plot_sub_MRT[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Malignant_rhabdoid_tumor_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'Malignant rhabdoid tumor')
points(umap_hv_rcc_plot_sub_nMRT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_MRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Medulloblastoma
umap_hv_rcc_plot_sub_nMedulloblastoma = subset(umap_hv_rcc_plot_sub, Lineage.subtype != 'Medulloblastoma')
umap_hv_rcc_plot_sub_nMedulloblastoma = umap_hv_rcc_plot_sub_nMedulloblastoma[,-3]
umap_hv_rcc_plot_sub_Medulloblastoma = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'Medulloblastoma')
umap_hv_rcc_plot_sub_Medulloblastoma = umap_hv_rcc_plot_sub_Medulloblastoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Medulloblastoma_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'Medulloblastoma')
points(umap_hv_rcc_plot_sub_nMedulloblastoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_Medulloblastoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## Melanoma
umap_hv_rcc_plot_sub_nMelanoma = subset(umap_hv_rcc_plot_sub, Lineage.subtype != 'Melanoma')
umap_hv_rcc_plot_sub_nMelanoma = umap_hv_rcc_plot_sub_nMelanoma[,-3]
umap_hv_rcc_plot_sub_Melanoma = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'Melanoma')
umap_hv_rcc_plot_sub_Melanoma = umap_hv_rcc_plot_sub_Melanoma[,-3]
pdf('Umap_of_Replicates_Lineage_Subtype_Melanoma_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'Melanoma')
points(umap_hv_rcc_plot_sub_nMelanoma, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_Melanoma, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
dev.off()

## CNS
CNS = subset(lineage, Lineage == 'central_nervous_system')
CNS$Lineage.subtype = str_replace(CNS$Lineage.subtype, '\\w{1}', toupper)
umap_hv_rcc_plot_sub_nCNS = subset(umap_hv_rcc_plot_sub, !(Lineage.subtype %in% CNS$Lineage.subtype))
umap_hv_rcc_plot_sub_Meningioma = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'Meningioma')
pdf('Umap_of_Replicates_Lineage_CNS_Subtypes_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'CNS')
points(umap_hv_rcc_plot_sub_nCNS, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_Glioma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_Medulloblastoma, col = rgb(0.7377778, 0.5244444, 0.1288889), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_Meningioma, col = rgb(0.3, 0.4, 0.35), pch = 19, cex = 1.7)
legend('bottomright', c('ATRT', 'Glioma', 'Medulloblastoma', 'Meningioma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980), c(0.7377778, 0.5244444, 0.1288889),
                       c(0.3, 0.4, 0.35))), pch = 19, pt.cex = 1.7)
dev.off()

## Embryonal brain tumors
umap_hv_rcc_plot_sub_nEBT = subset(umap_hv_rcc_plot_sub, Lineage.subtype != 'ATRT' | Lineage.subtype != 'Medulloblastoma')
umap_hv_rcc_plot_sub_nEBT = umap_hv_rcc_plot_sub_nEBT[,-3]
pdf('Umap_of_Replicates_Lineage_EBT_Subtypes_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'Embryonal Brain Tumors')
points(umap_hv_rcc_plot_sub_nEBT, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_ATRT, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_Medulloblastoma, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('bottomright', c('ATRT', 'Medulloblastoma'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

## Hematopoietic cancers
umap_hv_rcc_plot_sub_nhem = subset(umap_hv_rcc_plot_sub, !(Lineage.subtype %in% hem$Lineage.subtype))
umap_hv_rcc_plot_sub_blood = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'AML' | Lineage.subtype == 'CLL')
umap_hv_rcc_plot_sub_lymphocyte = subset(umap_hv_rcc_plot_sub, Lineage.subtype == 'Hodgkin lymphoma' |
                                           Lineage.subtype == 'Non hodgkin lymphoma')
pdf('Umap_of_Replicates_Lineage_Hematopoietic_Highly_Variable_Guides_Corrected_read_counts.pdf')
plot(umap_hv_rcc_plot, col = 'transparent', main = 'Hematopoietic Cancers')
points(umap_hv_rcc_plot_sub_nhem, col = rgb(0.7644444, 0.7511111, 0.7511111), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_blood, col = rgb(0.8, 0, 0), pch = 19, cex = 1.7)
points(umap_hv_rcc_plot_sub_lymphocyte, col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 19, cex = 1.7)
legend('bottomright', c('Blood', 'Lymphocyte'), bty = 'n',
       col = rgb(rbind(c(0.8, 0, 0), c(0.2117647, 0.3921569, 0.5450980))), pch = 19, pt.cex = 1.7)
dev.off()

write.xlsx(umap_c_plot, 'umap_c_plot.xlsx', rowNames = TRUE)
write.xlsx(umap_hv_c_plot, 'umap_hv_c_plot.xlsx', rowNames = TRUE)
write.xlsx(umap_rcc_plot, 'umap_rcc_plot.xlsx', rowNames = TRUE)
write.xlsx(umap_hv_rcc_plot, 'umap_hv_rcc_plot.xlsx', rowNames = TRUE)

## Overview over cell lines and cancer types
## count occurences of each type
lineage2 = lineage # data frame to change names

upper = c('AML', 'CLL', 'ATRT', 'NSCLC', 'SCLC')
lineage2$Lineage = str_to_sentence(lineage2$Lineage)
lineage2$Lineage.subtype = str_to_sentence(lineage2$Lineage.subtype)
lineage2$Lineage = gsub('_', ' ', lineage2$Lineage)
lineage2$Lineage.subtype = gsub('_', ' ', lineage2$Lineage.subtype)
lineage2$Lineage = gsub('Central nervous system', 'CNS', lineage2$Lineage)
lineage2$Lineage = gsub('Peripheral nervous system', 'PNS', lineage2$Lineage)

for (ii in 1:nrow(lineage2)) {
  if (lineage2$Lineage.subtype[ii] %in% str_to_sentence(upper)) {
    lineage2$Lineage.subtype[ii] = toupper(lineage2$Lineage.subtype[ii])
  }
}

mixed_colors_lin = rgb(t(col2rgb(colors_lin$inside)/255*matrix(nrow = 3, ncol = 16,
                                                               rep(c(rep(1.05, 8), rep(0.85, 8)), each = 3))))
names(mixed_colors_lin) = rownames(colors_lin)
mixed_colors_lin = mixed_colors_lin[order(names(mixed_colors_lin))]
## Change code from original function
PieDonutCustom <- function (data, mapping, start = getOption("PieDonut.start", 0), addPieLabel = TRUE,
                            addDonutLabel = TRUE, showRatioDonut = TRUE, showRatioPie = TRUE, ratioByGroup = TRUE,
                            showRatioThreshold = getOption("PieDonut.showRatioThreshold", 0.02),
                            labelposition = getOption("PieDonut.labelposition", 2), labelpositionThreshold = 0.1,
                            r0 = getOption("PieDonut.r0", 0.3), r1 = getOption("PieDonut.r1", 1),
                            r2 = getOption("PieDonut.r2", 1.2), explode = NULL, selected = NULL, explodePos = 0.1, 
                            color = "white", pieAlpha = 0.8, donutAlpha = 1, maxx = NULL, 
                            showPieName = TRUE, showDonutName = FALSE, title = NULL, 
                            pieLabelSize = 4, donutLabelSize = 3, titlesize = 5, explodePie = TRUE, 
                            explodeDonut = FALSE, use.label = TRUE, use.labels = TRUE, 
                            family = getOption("PieDonut.family", ""), palette_name="Dark2")
{
  (cols = colnames(data))
  if (use.labels) 
    data = moonBook::addLabelDf(data, mapping)
  count <- NULL
  if ("count" %in% names(mapping)) 
    count <- moonBook::getMapping(mapping, "count")
  count
  pies <- donuts <- NULL
  (pies = moonBook::getMapping(mapping, "pies"))
  if (is.null(pies)) 
    (pies = moonBook::getMapping(mapping, "pie"))
  if (is.null(pies)) 
    (pies = moonBook::getMapping(mapping, "x"))
  (donuts = moonBook::getMapping(mapping, "donuts"))
  if (is.null(donuts)) 
    (donuts = moonBook::getMapping(mapping, "donut"))
  if (is.null(donuts)) 
    (donuts = moonBook::getMapping(mapping, "y"))
  if (!is.null(count)) {
    df <- data %>% group_by(.data[[pies]]) %>% dplyr::summarize(Freq = sum(.data[[count]]))
    df
  }
  else {
    df = data.frame(table(data[[pies]]))
  }
  colnames(df)[1] = pies
  df$end = cumsum(df$Freq)
  df$start = dplyr::lag(df$end)
  df$start[1] = 0
  total = sum(df$Freq)
  df$start1 = df$start * 2 * pi/total
  df$end1 = df$end * 2 * pi/total
  df$start1 = df$start1 + start
  df$end1 = df$end1 + start
  df$focus = 0
  if (explodePie) 
    df$focus[explode] = explodePos
  df$mid = (df$start1 + df$end1)/2
  df$x = ifelse(df$focus == 0, 0, df$focus * sin(df$mid))
  df$y = ifelse(df$focus == 0, 0, df$focus * cos(df$mid))
  df$label = df[[pies]]
  df$ratio = df$Freq/sum(df$Freq)
  if (showRatioPie) {
    df$label = ifelse(df$ratio >= showRatioThreshold, paste0(df$label, 
                                                             "\n(", scales::percent(df$ratio), ")"), 
                      as.character(df$label))
  }
  df$labelx = (r0 + r1)/2 * sin(df$mid) + df$x
  df$labely = (r0 + r1)/2 * cos(df$mid) + df$y
  if (!is.factor(df[[pies]])) 
    df[[pies]] <- factor(df[[pies]])
  df
  mainCol = mixed_colors_lin # the line I had tochange
  df$radius = r1
  df$radius[df$focus != 0] = df$radius[df$focus != 0] + df$focus[df$focus != 
                                                                   0]
  df$hjust = ifelse((df$mid%%(2 * pi)) > pi, 1, 0)
  df$vjust = ifelse(((df$mid%%(2 * pi)) < (pi/2)) | (df$mid%%(2 * 
                                                                pi) > (pi * 3/2)), 0, 1)
  df$segx = df$radius * sin(df$mid)
  df$segy = df$radius * cos(df$mid)
  df$segxend = (df$radius + 0.05) * sin(df$mid)
  df$segyend = (df$radius + 0.05) * cos(df$mid)
  df
  if (!is.null(donuts)) {
    subColor = makeSubColor(mainCol, no = length(unique(data[[donuts]])))
    subColor
    data
    if (!is.null(count)) {
      df3 <- as.data.frame(data[c(donuts, pies, count)])
      colnames(df3) = c("donut", "pie", "Freq")
      df3
      df3 <- eval(parse(text = "complete(df3,donut,pie)"))
      df3$Freq[is.na(df3$Freq)] = 0
      if (!is.factor(df3[[1]])) 
        df3[[1]] = factor(df3[[1]])
      if (!is.factor(df3[[2]])) 
        df3[[2]] = factor(df3[[2]])
      df3 <- df3 %>% arrange(.data$pie, .data$donut)
      a <- df3 %>% spread(.data$pie, value = .data$Freq)
      a = as.data.frame(a)
      a
      rownames(a) = a[[1]]
      a = a[-1]
      a
      colnames(df3)[1:2] = c(donuts, pies)
    }
    else {
      df3 = data.frame(table(data[[donuts]], data[[pies]]), 
                       stringsAsFactors = FALSE)
      colnames(df3)[1:2] = c(donuts, pies)
      a = table(data[[donuts]], data[[pies]])
      a
    }
    a
    df3
    df3$group = rep(colSums(a), each = nrow(a))
    df3$pie = rep(1:ncol(a), each = nrow(a))
    total = sum(df3$Freq)
    total
    df3$ratio1 = df3$Freq/total
    df3
    if (ratioByGroup) {
      df3$ratio = scales::percent(df3$Freq/df3$group)
    }
    else {
      df3$ratio <- scales::percent(df3$ratio1)
    }
    df3$end = cumsum(df3$Freq)
    df3
    df3$start = dplyr::lag(df3$end)
    df3$start[1] = 0
    df3$start1 = df3$start * 2 * pi/total
    df3$end1 = df3$end * 2 * pi/total
    df3$start1 = df3$start1 + start
    df3$end1 = df3$end1 + start
    df3$mid = (df3$start1 + df3$end1)/2
    df3$focus = 0
    if (!is.null(selected)) {
      df3$focus[selected] = explodePos
    }
    else if (!is.null(explode)) {
      selected = c()
      for (i in 1:length(explode)) {
        start = 1 + nrow(a) * (explode[i] - 1)
        selected = c(selected, start:(start + nrow(a) - 
                                        1))
      }
      selected
      df3$focus[selected] = explodePos
    }
    df3
    df3$x = 0
    df3$y = 0
    df
    if (!is.null(explode)) {
      explode
      for (i in 1:length(explode)) {
        xpos = df$focus[explode[i]] * sin(df$mid[explode[i]])
        ypos = df$focus[explode[i]] * cos(df$mid[explode[i]])
        df3$x[df3$pie == explode[i]] = xpos
        df3$y[df3$pie == explode[i]] = ypos
      }
    }
    df3$no = 1:nrow(df3)
    df3$label = df3[[donuts]]
    if (showRatioDonut) {
      if (max(nchar(levels(df3$label))) <= 2) 
        df3$label = paste0(df3$label, "(", df3$ratio, 
                           ")")
      else df3$label = paste0(df3$label, "\n(", df3$ratio, 
                              ")")
    }
    df3$label[df3$ratio1 == 0] = ""
    df3$label[df3$ratio1 < showRatioThreshold] = ""
    df3$hjust = ifelse((df3$mid%%(2 * pi)) > pi, 1, 0)
    df3$vjust = ifelse(((df3$mid%%(2 * pi)) < (pi/2)) | (df3$mid%%(2 * 
                                                                     pi) > (pi * 3/2)), 0, 1)
    df3$no = factor(df3$no)
    df3
    labelposition
    if (labelposition > 0) {
      df3$radius = r2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$segx = df3$radius * sin(df3$mid) + df3$x
      df3$segy = df3$radius * cos(df3$mid) + df3$y
      df3$segxend = (df3$radius + 0.05) * sin(df3$mid) + 
        df3$x
      df3$segyend = (df3$radius + 0.05) * cos(df3$mid) + 
        df3$y
      if (labelposition == 2) 
        df3$radius = (r1 + r2)/2
      df3$labelx = (df3$radius) * sin(df3$mid) + df3$x
      df3$labely = (df3$radius) * cos(df3$mid) + df3$y
    }
    else {
      df3$radius = (r1 + r2)/2
      if (explodeDonut) 
        df3$radius[df3$focus != 0] = df3$radius[df3$focus != 
                                                  0] + df3$focus[df3$focus != 0]
      df3$labelx = df3$radius * sin(df3$mid) + df3$x
      df3$labely = df3$radius * cos(df3$mid) + df3$y
    }
    df3$segx[df3$ratio1 == 0] = 0
    df3$segxend[df3$ratio1 == 0] = 0
    df3$segy[df3$ratio1 == 0] = 0
    df3$segyend[df3$ratio1 == 0] = 0
    if (labelposition == 0) {
      df3$segx[df3$ratio1 < showRatioThreshold] = 0
      df3$segxend[df3$ratio1 < showRatioThreshold] = 0
      df3$segy[df3$ratio1 < showRatioThreshold] = 0
      df3$segyend[df3$ratio1 < showRatioThreshold] = 0
    }
    df3
    del = which(df3$Freq == 0)
    del
    if (length(del) > 0) 
      subColor <- subColor[-del]
    subColor
  }
  p <- ggplot() + ggforce::theme_no_axes() + coord_fixed()
  if (is.null(maxx)) {
    r3 = r2 + 0.3
  }
  else {
    r3 = maxx
  }
  p1 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", y0 = "y", 
                                             r0 = as.character(r0), r = as.character(r1), start = "start1", 
                                             end = "end1", fill = pies), alpha = pieAlpha, color = color, 
                                  data = df) + transparent() + scale_fill_manual(values = mainCol) + 
    xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
  if ((labelposition == 1) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                       y = "segy", xend = "segxend", yend = "segyend"), 
                            data = df) + geom_text(aes_string(x = "segxend", 
                                                              y = "segyend", label = "label", hjust = "hjust", 
                                                              vjust = "vjust"), size = pieLabelSize, data = df, 
                                                   family = family)
  }
  else if ((labelposition == 2) & (is.null(donuts))) {
    p1 <- p1 + geom_segment(aes_string(x = "segx", 
                                       y = "segy", xend = "segxend", yend = "segyend"), 
                            data = df[df$ratio < labelpositionThreshold, ]) + 
      geom_text(aes_string(x = "segxend", y = "segyend", 
                           label = "label", hjust = "hjust", 
                           vjust = "vjust"), size = pieLabelSize, 
                data = df[df$ratio < labelpositionThreshold, 
                ], family = family) + geom_text(aes_string(x = "labelx", 
                                                           y = "labely", label = "label"), size = pieLabelSize, 
                                                data = df[df$ratio >= labelpositionThreshold, ], 
                                                family = family)
  }
  else {
    p1 <- p1 + geom_text(aes_string(x = "labelx", y = "labely", 
                                    label = "label"), size = pieLabelSize, data = df, 
                         family = family)
  }
  if (showPieName) 
    p1 <- p1 + annotate("text", x = 0, y = 0, label = pies, 
                        size = titlesize, family = family)
  p1 <- p1 + theme(text = element_text(family = family))
  if (!is.null(donuts)) {
    if (explodeDonut) {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                 start = "start1", end = "end1", fill = "no", 
                                                 explode = "focus"), alpha = donutAlpha, 
                                      color = color, data = df3)
    }
    else {
      p3 <- p + ggforce::geom_arc_bar(aes_string(x0 = "x", 
                                                 y0 = "y", r0 = as.character(r1), r = as.character(r2), 
                                                 start = "start1", end = "end1", fill = "no"), 
                                      alpha = donutAlpha, color = color, data = df3)
    }
    p3 <- p3 + transparent() + scale_fill_manual(values = subColor) + 
      xlim(r3 * c(-1, 1)) + ylim(r3 * c(-1, 1)) + guides(fill = FALSE)
    p3
    if (labelposition == 1) {
      p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                         y = "segy", xend = "segxend", yend = "segyend"), 
                              data = df3) + geom_text(aes_string(x = "segxend", 
                                                                 y = "segyend", label = "label", hjust = "hjust", 
                                                                 vjust = "vjust"), size = donutLabelSize, 
                                                      data = df3, family = family)
    }
    else if (labelposition == 0) {
      p3 <- p3 + geom_text(aes_string(x = "labelx", 
                                      y = "labely", label = "label"), size = donutLabelSize, 
                           data = df3, family = family)
    }
    else {
      p3 <- p3 + geom_segment(aes_string(x = "segx", 
                                         y = "segy", xend = "segxend", yend = "segyend"), 
                              data = df3[df3$ratio1 < labelpositionThreshold, 
                              ]) + geom_text(aes_string(x = "segxend", 
                                                        y = "segyend", label = "label", hjust = "hjust", 
                                                        vjust = "vjust"), size = donutLabelSize, 
                                             data = df3[df3$ratio1 < labelpositionThreshold, 
                                             ], family = family) + geom_text(aes_string(x = "labelx", 
                                                                                        y = "labely", label = "label"), size = donutLabelSize, 
                                                                             data = df3[df3$ratio1 >= labelpositionThreshold, 
                                                                             ], family = family)
    }
    if (!is.null(title)) 
      p3 <- p3 + annotate("text", x = 0, y = r3, 
                          label = title, size = titlesize, family = family)
    else if (showDonutName) 
      p3 <- p3 + annotate("text", x = (-1) * r3, 
                          y = r3, label = donuts, hjust = 0, size = titlesize, 
                          family = family)
    p3 <- p3 + theme(text = element_text(family = family))
    grid::grid.newpage()
    print(p1, vp = grid::viewport(height = 1, width = 1))
    print(p3, vp = grid::viewport(height = 1, width = 1))
  }
  else {
    p1
  }
}

pdf('PieDonut_Lineages_with_All_Names.pdf', 15,15)
PieDonutCustom(lineage2, aes(Lineage, Lineage.subtype), explode = 1:16, color = 'black', showRatioDonut = FALSE,
               showRatioPie = FALSE)
dev.off()
pdf('PieDonut_Lineages_with_Subtype_Names.pdf', 15,15)
PieDonutCustom(lineage2, aes(Lineage, Lineage.subtype), explode = 1:16, color = 'black', addPieLabel = FALSE,
               addDonutLabel = FALSE, use.label = FALSE, use.labels = FALSE, showRatioDonut = FALSE,
               showRatioPie = FALSE, pieLabelSize = 0)
dev.off()
pdf('PieDonut_Lineages_without_Names.pdf', 15,15)
PieDonutCustom(lineage2, aes(Lineage, Lineage.subtype), explode = 1:16, color = 'black', addPieLabel = FALSE,
               addDonutLabel = FALSE, use.label = FALSE, use.labels = FALSE, showRatioDonut = FALSE,
               showRatioPie = FALSE, pieLabelSize = 0, donutLabelSize = 0, labelpositionThreshold = 0)
dev.off()

n_disease = table(lineage$Primary.disease)
pie(n_disease, col = c(brewer.pal(9, 'Set1'), brewer.pal(8, 'Pastel2')))
n_lineage = table(lineage$Lineage)
pie(n_lineage, col = c(brewer.pal(8, 'Set1'), brewer.pal(8, 'Pastel2')))
n_subtype = table(lineage$Lineage.subtype)
n_subsubtype = table(lineage$Lineage.sub.subtype)
frequency = list(Disease = n_disease, Lineage = n_lineage, Subtype = n_subtype, Subsubtype = n_subsubtype)

write.xlsx(frequency, 'Frequencies of cell lines.xlsx', colNames = FALSE)

library(scales)
library(plotrix)
## Read corrected read counts and LFCs only with samples that will be included in the final analysis
counts_gene_replicates_wns = read.xlsx("directory/Corrected read counts on gene level with nosites.xlsx", rowNames = TRUE)
counts_gene_replicates = read.xlsx("directory/Corrected read counts on gene level without nosites.xlsx", rowNames = TRUE)

counts_guide_replicates = read.xlsx("directory/Corrected read counts without bad cell lines.xlsx")
rownames(counts_guide_replicates) = counts_guide_replicates$sgRNA
counts_guide_gene = counts_guide_replicates[,(1:2)]
counts_guide_replicates = counts_guide_replicates[,-(1:2)]

lfc_gene_replicates = read.xlsx("directory/Corrected LFC on gene level without bad cell lines.xlsx", rowNames = TRUE)

lfc_guide_replicates = read.xlsx("directory/Corrected LFC without bad cell lines.xlsx")
rownames(lfc_guide_replicates) = lfc_guide_replicates$sgRNA
lfc_guide_gene = lfc_guide_replicates[,(1:2)]
lfc_guide_replicates = lfc_guide_replicates[,-(1:2)]

pos_cntrl = read.csv('directory/Control_essentials.csv', header = FALSE)

## Aggregate over cell lines
## Read counts on gene level with no sites
tcounts_gene_replicates_wns = as.data.frame(t(counts_gene_replicates_wns)) # transpose to bring replicates in row
tcounts_gene_replicates_wns = cbind(Cell_line = sub('_REP_.', '', rownames(tcounts_gene_replicates_wns)),
                                    tcounts_gene_replicates_wns) # generate col with line names to aggregate over
tcounts_gene_line_wns = aggregate(. ~ Cell_line , mean, data = tcounts_gene_replicates_wns) # aggregate
rownames(tcounts_gene_line_wns) = tcounts_gene_line_wns$Cell_line # make col with cell line names to rownames
tcounts_gene_line_wns = tcounts_gene_line_wns[,-1] # and delete this col to enable transposition
counts_gene_line_wns = as.data.frame(t(tcounts_gene_line_wns)) # transpose to get original format

## Read counts on gene level without no sites
tcounts_gene_replicates = as.data.frame(t(counts_gene_replicates))
tcounts_gene_replicates = cbind(Cell_line = sub('_REP_.', '', rownames(tcounts_gene_replicates)),
                                tcounts_gene_replicates)
tcounts_gene_line = aggregate(. ~ Cell_line , mean, data = tcounts_gene_replicates)
rownames(tcounts_gene_line) = tcounts_gene_line$Cell_line
tcounts_gene_line = tcounts_gene_line[,-1]
counts_gene_line = as.data.frame(t(tcounts_gene_line))

## Read counts on guide level
tcounts_guide_replicates = as.data.frame(t(counts_guide_replicates))
tcounts_guide_replicates = cbind(Cell_line = sub('_REP_.', '', rownames(tcounts_guide_replicates)),
                                 tcounts_guide_replicates)
tcounts_guide_line = aggregate(. ~ Cell_line , mean, data = tcounts_guide_replicates)
rownames(tcounts_guide_line) = tcounts_guide_line$Cell_line
tcounts_guide_line = tcounts_guide_line[,-1]
counts_guide_line = cbind.data.frame(counts_guide_gene,t(tcounts_guide_line))

## LFCs on gene level
tlfc_gene_replicates = as.data.frame(t(lfc_gene_replicates))
tlfc_gene_replicates = cbind(Cell_line = sub('_REP_.', '', rownames(tlfc_gene_replicates)),
                             tlfc_gene_replicates)
tlfc_gene_line = aggregate(. ~ Cell_line , mean, data = tlfc_gene_replicates)
rownames(tlfc_gene_line) = tlfc_gene_line$Cell_line
tlfc_gene_line = tlfc_gene_line[,-1]
lfc_gene_line = as.data.frame(t(tlfc_gene_line))

## LFCs on guide level
tlfc_guide_replicates = as.data.frame(t(lfc_guide_replicates))
tlfc_guide_replicates = cbind(Cell_line = sub('_REP_.', '', rownames(tlfc_guide_replicates)),
                              tlfc_guide_replicates)
tlfc_guide_line = aggregate(. ~ Cell_line , mean, data = tlfc_guide_replicates)
rownames(tlfc_guide_line) = tlfc_guide_line$Cell_line
tlfc_guide_line = tlfc_guide_line[,-1]
lfc_guide_line = cbind.data.frame(lfc_guide_gene, t(tlfc_guide_line))

## Sanity Checks
SanityCheck1 = all(colnames(counts_gene_line_wns) == colnames(counts_gene_line)) &&
  all(colnames(counts_guide_line)[-c(1:2, which(colnames(counts_guide_line) == "spike_in_plasmid_DNA_from_CP1608"))] ==
        colnames(lfc_gene_line)) && all(colnames(lfc_guide_line)[-(1:2)] == colnames(lfc_gene_line)) &&
  all(colnames(counts_guide_line)[-(1:2)] == colnames(counts_gene_line))
SanityCheck2 = all(rownames(lfc_guide_line) == lfc_guide_line$sgRNA) &&
  all(rownames(counts_guide_line) == counts_guide_line$sgRNA)
SanityCheck3 = all(rownames(lfc_gene_line) == rownames(counts_gene_line))

# SanityCheck1 # Are all data frames reduced to the same cell lines?
# SanityCheck2 # Do the guide names of the data frame on guide level after aggregation match them before aggregation?
# SanityCheck3 # Do the gene names of the data frames on gene level without nosites match?
all(SanityCheck1,SanityCheck2,SanityCheck3)

write.xlsx(lfc_guide_line, 'Log fold changes on guide and cell line level.xlsx', rowNames= TRUE)
write.xlsx(lfc_gene_line, 'Log fold changes on gene and cell line level.xlsx', rowNames= TRUE)
write.xlsx(counts_guide_line, 'Read counts on guide and cell line level.xlsx', rowNames= TRUE)
write.xlsx(counts_gene_line, 'Read counts on gene and cell line level.xlsx', rowNames= TRUE)

## Plot cumulative frequency of every cell line LFC
## LFC on gene level
pdf('Cumulative_Frequency_of_LFC_Cell_Lines_vs_Library_Gene_Level.pdf')
plot(ecdf(lfc_gene_line$A204), bty = 'n', xlim = c((min(lfc_gene_line)-0.1),(max(lfc_gene_line)+0.1)), main = '',
     pch = 20, lty = 0, ylab = 'Cumulative frequency',
     xlab = "Mean LFC of genes", col.01line = NULL)
for (ii in 1:ncol(lfc_gene_line)) {
  lines(ecdf(lfc_gene_line[,ii]), pch = 20, lty = 0, col = rgb(0.8, 0, 0), do.points = TRUE, cex = 0.4, 
        col.01line = NULL)
}
abline(v = 0, col = 'black', lwd = 2)
legend('bottomright', c('Plasmid library', 'Cell lines'), col = rgb(rbind(c(0, 0, 0), c(0.8, 0, 0))), pch = 20,
       bty = 'n', pt.cex = 1.3)
dev.off()

pdf('Cumulative_Frequency_of_LFC_Cell_Lines_vs_Library_Gene_Level_Zoom.pdf')
plot(ecdf(lfc_gene_line$A204), bty = 'n', xlim = c(-1,1), main = '', pch = 20, lty = 0,
     ylab = 'Cumulative frequency', xlab = "Mean LFC of genes", col.01line = NULL)
for (ii in 1:ncol(lfc_gene_line)) {
  lines(ecdf(lfc_gene_line[,ii]), pch = 20, lty = 0, col = rgb(0.8, 0, 0), do.points = TRUE, cex = 0.4, 
        col.01line = NULL)
}
abline(v = 0, col = 'black', lwd = 2)
legend('bottomright', c('Plasmid library', 'Cell lines'), col = rgb(rbind(c(0, 0, 0), c(0.8, 0, 0))), pch = 20,
       bty = 'n', pt.cex = 1.3)
dev.off()

## Plots with read counts on guide level

counts_log_guide = log2(counts_guide_line[,-(1:2)]+1)
rownames(counts_log_guide) = counts_guide_line$sgRNA
counts_plot_guide = counts_log_guide[,-(which(colnames(counts_log_guide) == 'spike_in_plasmid_DNA_from_CP1608'))]

pdf('Cumulative_Frequency_of_Read_Counts_Guide_Level_Cell_Lines_vs_Library.pdf')
plot(ecdf(counts_log_guide$spike_in_plasmid_DNA_from_CP1608), bty = 'n', main = '', pch = 20, lty = 0,
     xlim = c(min(counts_plot_guide), 20), ylab = 'Cumulative frequency',
     xlab = expression(paste(log[2], " read count of sgRNAs", sep = "")), col.01line = NULL)
for (ii in 3:ncol(counts_plot_guide)) {
  lines(ecdf(counts_plot_guide[,ii]), pch = 20, lty = 0, col = rgb(0.8,0,0), do.points = TRUE, cex = 0.4, 
        col.01line = NULL)
}
lines(ecdf(counts_log_guide$spike_in_plasmid_DNA_from_CP1608), pch = 20, lty = 0, col = c(1, 1, 1),
      do.points = TRUE, cex = 0.4, col.01line = NULL)
legend('bottomright', c('Plasmid library', 'Cell lines'), col = rgb(rbind(c(0, 0, 0), c(0.8, 0, 0))),
       pch = 20, bty = 'n', pt.cex = 1.3)
dev.off()

pdf('Cumulative_Frequency_of_Read_Counts_Guide_Level_Cell_Lines_vs_Library_With_Zoom.pdf')
plot(ecdf(counts_log_guide$spike_in_plasmid_DNA_from_CP1608), bty = 'n', main = '', pch = 20, lty = 0,
     xlim = c(min(counts_plot_guide), 20), ylab = 'Cumulative frequency',
     xlab = expression(paste(log[2], " read count of sgRNAs", sep = "")), col.01line = NULL)
for (ii in 3:ncol(counts_plot_guide)) {
  lines(ecdf(counts_plot_guide[,ii]), pch = 20, lty = 0, col = rgb(0.8,0,0), do.points = TRUE, cex = 0.4, 
        col.01line = NULL)
}
lines(ecdf(counts_log_guide$spike_in_plasmid_DNA_from_CP1608), pch = 20, lty = 0, col = c(1, 1, 1),
      do.points = TRUE, cex = 0.4, col.01line = NULL)
legend('bottomright', c('Plasmid library', 'Cell lines'), col = rgb(rbind(c(0, 0, 0), c(0.8, 0, 0))),
       pch = 20, bty = 'n', pt.cex = 1.3)
segments(c(8.2, 10.05, 10.05, 8.2, 8.2, 10.05), c(-0.01, -0.01, 0.6, 0.6, -0.01, 0.6),
         c(10.05, 10.05, 8.2, 8.2, 0, 8), c(-0.01, 0.6, 0.6, -0.01, 0.25, 1),
         col = rgb(0, 0, 0, alpha = 0.6))

par(fig = c(grconvertX(c(0, 8), from = "user", to = "ndc"), # define placement of the box inside the first plot
            grconvertY(c(0.25, 1), from = "user", to = "ndc")),
    mar = c(2.1,2.1,0.5,0.5), # define space between plot and edges of box
    new = TRUE) # plot on top of first plot

plot(ecdf(counts_log_guide$spike_in_plasmid_DNA_from_CP1608), main = '', bty = 'n', xlim = c(8.2, 10.05),
     ylim = c(0, 0.6), pch = 20, lty = 0, ylab = '', xlab = '', col.01line = NULL)
for (ii in 1:ncol(counts_plot_guide)) {
  lines(ecdf(counts_plot_guide[,ii]), pch = 20, lty = 0, col = rgb(0.8,0,0), do.points = TRUE, cex = 0.4, 
        col.01line = NULL)
}
lines(ecdf(counts_log_guide$spike_in_plasmid_DNA_from_CP1608), pch = 20, lty = 0, col = c(1, 1, 1), do.points = TRUE,
      cex = 0.4, col.01line = NULL)
box('figure')
dev.off()

## Plot quality of sgRNAs
genes = subset(lfc_guide_line, !grepl('ONE_INTERGENIC', lfc_guide_line$gene))
miRNA = subset(genes, !(gene %in% pos_cntrl$V1))
mult_targets = subset(names(table(miRNA$gene)), table(miRNA$gene) > 1)
mean = aggregate(miRNA[,-(1:2)], list(miRNA$gene), mean)
colnames(mean)[1] = 'gene'
firstq = aggregate(miRNA[,-(1:2)], list(miRNA$gene), FUN = 'quantile', probs = 0.25)
colnames(firstq)[1] = 'gene'
thirdq = aggregate(miRNA[,-(1:2)], list(miRNA$gene), FUN = 'quantile', probs = 0.75)
colnames(thirdq)[1] = 'gene'
mean_m = subset(mean, gene %in% mult_targets)
firstq_m = subset(firstq, gene %in% mult_targets)
thirdq_m = subset(thirdq, gene %in% mult_targets)
minimum = aggregate(miRNA[,-(1:2)], list(miRNA$gene), min)
colnames(minimum)[1] = 'gene'
minimum_m = subset(minimum, gene %in% mult_targets)
maximum = aggregate(miRNA[,-(1:2)], list(miRNA$gene), max)
colnames(maximum)[1] = 'gene'
maximum_m = subset(maximum, gene %in% mult_targets)
lab = c(1, 5, 10, 20, 50, 100, 200, 500, 1000, 2000)
for (ii in 3:ncol(mean)) {
  plot = data.frame(Gene = mean_m$gene, Mean = mean_m[,ii], quant25 = firstq_m[,ii],
                    quant75 = thirdq_m[,ii], minimum = minimum_m[,ii], maximum =  maximum_m[,ii])
  plot_order = plot[order(plot$minimum),]
  plot_order = cbind.data.frame(plot_order, n = 1:nrow(plot_order))
  pdf(paste('Quality_Library_', colnames(mean)[ii], '.pdf', sep =''))
  plot(cbind(plot_order$n, plot_order$quant25), col = 'transparent', bty = 'n', xlim = c(1,2000), xlab = 'Genes',
       ylab = 'LFC', main = colnames(mean)[ii], ylim = c(min(plot_order$minimum), max(plot_order$maximum)),
       log = 'x')
  # segments(log(plot_order$n, 20), plot_order$quant25, (log(plot_order$n, 20)), plot_order$quant75,
  # col = rgb(0.7291, 0.6588, 0.4549, alpha = 0.6), lwd = 2)
  # points(cbind(plot_order$n, plot_order$minimum), pch = 19, cex = 0.5, col = rgb(0, 0, 0,
  # alpha = 0.6))
  points(cbind(plot_order$n, plot_order$Mean), pch = 20, col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = 0.4))
  abline(h = 0, col = rgb(0.8, 0, 0, alpha = 0.7), lwd = 2)
  lines(plot_order$n, predict(loess(plot_order$minimum ~ plot_order$n, span = 0.15)),
        lwd = 3, col = 'black')
  # points(cbind(log(plot_order$n, 20), plot_order$maximum), pch = 20, cex = 0.5,
  # col = rgb((col2rgb('grey')/225)[1,1], (col2rgb('grey')/225)[2,1], (col2rgb('grey')/225)[3,1],
  # alpha = 0.6))
  lines(plot_order$n, predict(loess(plot_order$maximum ~ plot_order$n, span = 0.15)),
        col = 'grey', lwd = 3)
  
  lines(plot_order$n, predict(loess(plot_order$Mean ~ plot_order$n, span = 0.15)),
        col = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 3)
  abline(h = 0, col = rgb(0.8, 0, 0, alpha = 0.5), lwd = 2)
  legend('bottomright',
         c('Minimum LFC of miRNA genes', 'Mean LFC of miRNA genes', 'Maximum LFC of miRNA genes'),
         lty = 1, bty = 'n', lwd = 5,
         col = rgb(rbind(c(0, 0, 0), c(0.2117647, 0.3921569, 0.5450980), c(col2rgb('grey')/225))))
  dev.off()
  
  # pdf(paste('Quality_Library_Min_Max', colnames(mean)[ii], '.pdf', sep = ''))
  # plot(cbind(log(1:length(sortedmin), 20), sortedmin), col = 'transparent', bty = 'n',
  #      xlim = log(c(1,2000), 20), xlab = 'Genes', ylab = 'LFC', main = colnames(mean)[ii], xaxt = 'n', 
  #      ylim = c(min(sortedmin, sortedmax, sortedm$Mean), max(sortedmin, sortedmax, sortedm$Mean))) 
  # axis(1, at = log(lab, 20), labels = lab)
  # segments(log(1:length(sortedmin), 20), sortedmin, (log(1:length(sortedmax), 20)), sortedmax,
  #          col = rgb(0.7291, 0.6588, 0.4549, alpha = 0.6), lwd = 2)
  # points(cbind(log(1:nrow(sortedm), 20), sortedm$Mean), pch = 20, cex = 0.5,
  #        col = rgb(0.2117647, 0.3921569, 0.5450980))
  # abline(h = 0, col = rgb(0.8, 0, 0, alpha = 0.7), lwd = 2)
  # legend('bottomright', c('Mean LFC of miRNA genes', 'Span from max to min'), pch =  c(20, 124), bty = 'n',
  #        col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.7291, 0.6588, 0.4549))))
  # dev.off()
  
  pdf(paste('Inner_50_Quantile_Mean_', colnames(mean)[ii], '.pdf', sep = ''))
  plot(plot_order$Mean, (plot_order$quant25-plot_order$quant75), pch = 20, xlab = 'Mean LFC',
       ylab = 'Span of inner 50 % quantile', bty = 'n', main = colnames(mean_m)[ii])
  dev.off()
}


####Perform BAGEL2 and MAGeCK RRA in Python Jupyter Notebook



###continue with BAGEL2 and MAGeCK data here
rm(list = ls())
library(devtools)
library(ggplot2)
library(CoRe)
library(ggvenn)
library(ggridges)
library(data.table)
library(pheatmap)
library(Cairo)
getwd()
setwd("directory")
pos_cntrl = read.csv('directory/Control_essentials.csv', header = FALSE)
anno_clean = read.csv('directory/CP1608_anno_clean.csv', header = TRUE)
anno_genes = subset(anno_clean, !grepl('ONE_INTERGENIC_SITE_', GENES, fixed = TRUE)) # remove negative controls
LFC = read.xlsx("directory/Log fold changes on gene and cell line level.xlsx", rowNames = TRUE)
LFC_guide = read.xlsx("directory/Log fold changes on guide and cell line level.xlsx", rowNames = TRUE)

## Output from BAGEL2
setwd('directory to BAGEL2 pr results')
files_bagel = list.files() # get all files in output directory
bagel = vector(mode = "list") # initialise list to hold a data frame for each cell line
jj = 1 # set counter
for (ii in 1:length(files_bagel)) { # for all files in the folder CAUTION: ii does not count txt files but all files! 
  if (grepl('.txt', files_bagel[ii], fixed = TRUE) && # only the txt files
      (sub('_pr\\.txt', '', files_bagel[ii]) != 'Jurkat') && (sub('_pr\\.txt', '', files_bagel[ii]) != 'MDAMB453')) { # not excluded cell lines 
    bagel[[jj]] = read.table(files_bagel[ii], header = TRUE) # add data frame for each cell line to list
    names(bagel)[jj] = sub('_pr\\.txt', '', files_bagel[ii]) # name added data frame according to filename
    jj = jj + 1
  }
}

## Define essential genes by FDR < 0.1 and split essential and nonessential genes
lowest_BF = data.frame() # initialise data frame and lists to hold elements
bagel_essentials = vector(mode = 'list')
bagel_five = vector(mode = 'list')
bagel_ness = vector(mode = 'list')
for (ii in 1:length(bagel)) { # for every matrix in the list (for every cell line)
  bagel_essentials[[ii]] = subset(bagel[[ii]], FDR < 0.1) # take all genes with FDR < 0.1 as essentials
  bagel_five[[ii]] = subset(bagel[[ii]], FDR < 0.05)
  names(bagel_essentials)[ii] = names(bagel)[ii] # name the matrix after cell line
  names(bagel_five)[ii] = names(bagel)[ii]
  lowest_BF = rbind(lowest_BF, cbind(names(bagel)[ii], bagel_essentials[[ii]][nrow(bagel_essentials[[ii]]),])) # take last element of essentials and save the information
  bagel_ness[[ii]] = subset(bagel[[ii]], FDR >= 0.1) # take all other genes in the list of nonessentials
  names(bagel_ness)[ii] = names(bagel)[ii] # name the matrix after cell line
}
colnames(lowest_BF)[1] = 'Cell_Line'
lowest_BF = cbind.data.frame(lowest_BF, n = as.numeric(rownames(lowest_BF))) # rownames are the rownames of the row in the old matrix. Hence, it represents the nth gene in the ranking.

setwd('directory to MAGeCK RRA results')
files_mageck = list.files() # same for mageck
mageck = vector(mode = 'list')
jj = 1
for (ii in 1:length(files_mageck)) {
  if (grepl('gene_summary.txt', files_mageck[ii], fixed = TRUE)){ # excluded cell lines not run, don't have to be excluded
    mageck[[jj]] = read.table(files_mageck[ii], header = TRUE)
    if (sub('\\.gene_summary\\.txt', '', files_mageck[ii]) == 'BENMEN'){
      names(mageck)[jj] = sub('BENMEN', 'BEN-MEN', sub('\\.gene_summary\\.txt', '', files_mageck[ii]))
    }
    else{
      names(mageck)[jj] = sub('\\.gene_summary\\.txt', '', files_mageck[ii])
    }
    jj = jj + 1
  }
}

## Define essential genes by negative FDR < 0.1 and split essential and nonessential genes
lowest_mageck = data.frame()
mageck_essentials = vector(mode = 'list')
mageck_five = vector(mode = 'list')
mageck_ness = vector(mode = 'list')
for (ii in 1:length(mageck)) {
  mageck_essentials[[ii]] = subset(mageck[[ii]], neg.fdr < 0.1) # take negative FDR for essentials
  mageck_five[[ii]] = subset(mageck[[ii]], neg.fdr < 0.05)
  names(mageck_essentials)[ii] = names(mageck)[ii]
  names(mageck_five)[ii] = names(mageck)[ii]
  lowest_mageck = rbind(lowest_mageck,
                        cbind(names(mageck)[ii], mageck_essentials[[ii]][nrow(mageck_essentials[[ii]]),]))
  mageck_ness[[ii]] = subset(mageck[[ii]], neg.fdr >= 0.1 & pos.fdr >= 0.1) # take both FDRs for nonessentials since pos.fdr < 0.1 could correspond to tumor suppressor genes.
  names(mageck_ness)[ii] = names(mageck)[ii]
}
colnames(lowest_mageck)[1] = 'Cell_Line'
lowest_mageck = cbind.data.frame(lowest_mageck, n = as.numeric(rownames(lowest_mageck))) # Again, rowname comes from old matrix and corresponds to the nth gene 

## Get the intersect of the essential genes from BAGEL and MAGeCK
essentials = vector(mode = 'list')
five = vector(mode = 'list')
information_essentials = vector(mode = 'list')
information_five = vector(mode = 'list')
for (ii in 1:length(bagel_essentials)) {
  essentials[[ii]] = intersect(bagel_essentials[[ii]]$Gene,mageck_essentials[[ii]]$id) # regard essentials only if MAGeCK and BAGEL both identified them
  five[[ii]] = intersect(bagel_five[[ii]]$Gene, mageck_five[[ii]]$id) # same for 5% FDR
  names(essentials)[ii] = names(bagel_essentials)[ii] # name the matrix after the cell line
  names(five)[ii] = names(bagel_five)[ii]
  for (jj in 1:length(essentials[[ii]])) {
    if (jj == 1) { # has to be done in 2 steps because element has to be initialised with first element (? any better method?)
      information_essentials[[ii]] = # initialise
        cbind(bagel_essentials[[ii]][which(bagel_essentials[[ii]][,1] == essentials[[ii]][jj]),], # take rows of identified essentials with information from BAGEL
              mageck_essentials[[ii]][which(mageck_essentials[[ii]][,1] == essentials[[ii]][jj]),]) # and bind the cols with same from MAGeCK
    } else {
      information_essentials[[ii]][jj,] = # fill with rest of information, same logic as above
        cbind(bagel_essentials[[ii]][which(bagel_essentials[[ii]][,1] == essentials[[ii]][jj]),],
              mageck_essentials[[ii]][which(mageck_essentials[[ii]][,1] == essentials[[ii]][jj]),])
    }
  }
  names(information_essentials)[ii] = names(bagel_essentials)[ii] #name matrix after cell line
  for (kk in 1:length(five[[ii]])) {
    if (kk == 1) {
      information_five[[ii]] =
        cbind(bagel_five[[ii]][which(bagel_five[[ii]][,1] == five[[ii]][kk]),],
              mageck_five[[ii]][which(mageck_five[[ii]][,1] == five[[ii]][kk]),])
    } else {
      information_five[[ii]][kk,] = 
        cbind(bagel_five[[ii]][which(bagel_five[[ii]][,1] == five[[ii]][kk]),],
              mageck_five[[ii]][which(mageck_five[[ii]][,1] == five[[ii]][kk]),])
    }
  }
  names(information_five)[ii] = names(bagel_five)[ii]
}
common_essentials = essentials[[1]] # initialise data frame for intersects with first set
for (ii in 2:length(essentials)) { # intersects of essential genes from all cell lines, start with second because first was used to initialise
  common_essentials = intersect(common_essentials, essentials[[ii]]) # replace data frame with the intersect of itself with the next set
}

## Get genes not marked as essential in any screen by any method
nonessentials = intersect(bagel_ness$A204$Gene, mageck_ness$A204$id) # same logic as last loop
for (ii in 2:length(bagel_ness)) {
  nonessentials = intersect(nonessentials, intersect(bagel_ness[[ii]]$Gene, mageck_ness[[ii]]$id))
}
## Determine number of positive controls (known essentials) in essentials of each cell line
n_cntrl = vector() # initialise vector to hold numbers
for (ii in 1:length(essentials)) { # check each matrix (each cell line)
  counter = 0 # counts essential genes, is reset to 0 for every cell line here 
  for (jj in 1:length(essentials[[ii]])) { # check all elements of the current matrix with essential genes of this cell line
    if (essentials[[ii]][jj] %in% pos_cntrl$V1) { # if it is included in the data frame with positive controls
      counter = counter + 1 # count it
    }
  }
  n_cntrl[ii] = counter # after checking all elements of this one cell line, store the sum of detected positive controls
  names(n_cntrl)[ii] = names(essentials)[ii] # name matrix after cell line
}

## Binary matrix
## ADaM need a binary matrix as input where 1 marks an essential gene in the corresponding cell line
## Columns are cell lines and rows are genes
binmat = matrix(NA, nrow = length(unique(anno_genes$GENES)), ncol = length(essentials))
colnames(binmat) = names(essentials) # name cols as cell lines
rownames(binmat) = unique(anno_genes$GENES) # name rows as all genes
for (ii in 1:length(essentials)) { # ii for cell lines
  for (jj in 1:nrow(binmat)) { # jj for genes
    if (rownames(binmat)[jj] %in% essentials[[ii]]){
      binmat[jj,ii] = 1 # 1 if gene is essential in this cell line
    }
    else {
      binmat[jj,ii] = 0 # otherwise 0
    }
  }
}
binmat = as.data.frame(binmat) # format of data frame is needed. Initialised as matrix to initialise as empty matrix with predetermined dimensions

## Sanity Checks
# Is any of the nonessentials in any essentials list of the BAGEL output?

SanityCheck1_1 = all(names(bagel) == names(mageck))
SanityCheck1_2 = all(names(bagel_essentials) == names(mageck_essentials))
SanityCheck1_3 = all(names(bagel) == names(bagel_essentials))
SanityCheck1_4 = TRUE
for (ii in 1:length(bagel_essentials)){
  SanityCheck1_4 = SanityCheck1_4 && !any(nonessentials %in% bagel_essentials[[ii]]$Gene) ||
    !(any(nonessentials %in% mageck_essentials[[ii]]$id))
}

# SanityCheck1_1 # Are the cell lines the same for BAGEL and MAGeCK
# SanityCheck1_2 # Are the cell lines the same for the lists containing essential genes from BAGEL and MAGeCK
# SanityCheck1_3 # Are the cell lines the same for the complete lists and lists containing essential genes?
# SanityCheck1_4 # Are all nonessential genes not contained in the lists of essential genes of MAGeCK or BAGEL for any cell line?
all(SanityCheck1_1,SanityCheck1_2,SanityCheck1_3,SanityCheck1_4)

## CoRe identify essentials

## ADaM
pos_cntrl_ADaM = pos_cntrl$V1 # transform data frame into vector of characters (needed as input for ADaM)

pdf('Profiles_of_Fitness_Genes_CoRe.pdf') # save output figure in pdf
profile = CoRe.panessprofile(binmat) # profile with fitness genes in n cell lines
dev.off()
pdf('Null_Model_CoRe.pdf')
nullmod = CoRe.generateNullModel(binmat) # nullmodel by random permutations of binary matrix
dev.off()
odds = CoRe.empiricalOdds(observedCumSum = profile$CUMsums, simulatedCumSum = nullmod$nullCumSUM) # odds ratio
truePR = CoRe.truePositiveRate(binmat, pos_cntrl_ADaM) # TPR depending on the number of cell lines genes have to be essential 
# $P is the number of genes that are essential in >= n of the cell lines
# $TP is the number of positive controls that are in the sets that are essential in >= n cell lines
# $TPR is the TPR (detected positive controls) depending on n, is $TP/10 (number of true positives divided by number of positive controls)

## Change function so that plot does not contain label of 100 % covered
CoRe.tradeoffEO_TPR2<-function(EO,TPR,test_set_name,display=TRUE){
  x<-EO
  x[x==Inf]<-max(x[x<Inf])
  x<-(x-min(x))/(max(x)-min(x))
  
  y<-TPR
  y<-(y-min(y))/(max(y)-min(y))
  
  orEO<-EO
  orEO[orEO==Inf]<-max(orEO[orEO<Inf])
  orTPR<-TPR
  
  EO<-x
  TPR<-y
  point<-min(which(!y>x))
  
  if(display){
    CCOL<-'red'
    par(mar=c(4,4,4,4))
    par(mfrow=c(1,1))
    MAIN<-c('log10 (obs/Expct) n.genes [red, left]',
            paste('% covered ',test_set_name,' [blue, right]',sep=''))
    plot(EO,type='l',xlab='genes depleted in >= # cell lines',ylab='',axes=FALSE,lwd=4,main=MAIN,col=CCOL,cex.main=0.8,
         xlim=c(0,length(EO)))
    axis(2,at = seq(0,1,0.2),format(seq(min(orEO),max(orEO),(max(orEO)-min(orEO))/5),digits=2),col='red')
    axis(1)
    par(new=TRUE)
    plot(TPR,type='l',xlab='',ylab='',axes=FALSE,lwd=4,col='blue',ylim=c(0,1),xlim=c(0,length(EO)))
    axis(4,at = seq(0,1,0.2),format(seq(min(orTPR),max(orTPR),(max(orTPR)-min(orTPR))/5),digits=2),col='blue')
    
    
    
    abline(v=point)
    abline(h=y[point],lty=2)
    
    points(point,y[point],pch=16,cex=2)
    
    # legend('top',paste(format(100*orTPR[point],digits=2),'% covered',sep=''),bg = NULL,bty = 'n')
  }
  
  return(point)
}

pdf('Tradeoff_Empirical_Odds_TPR.pdf', 10, 9)
crossoverpoint = CoRe.tradeoffEO_TPR2(odds, truePR$TPR, test_set_name = 'CP1608 essentials') # tradeoff between odds ratio and TPR to determine threshold of n of cell lines in which a gene has to be depleted to be considered a CFG
dev.off()
CFG = CoRe.coreFitnessGenes(binmat, crossoverpoint) # get all genes that meet the ADaM threshold
CFG_cntrl = subset(CFG, CFG %in% pos_cntrl_ADaM) # split CFG into the previously known essential genes
CFG_other = subset(CFG, !(CFG %in% pos_cntrl_ADaM)) # and new essential genes
adam = CoRe.ADaM(binmat, TruePositives =  pos_cntrl_ADaM, display = TRUE) # can be done with one command but has less output, only stores CFGs

## 90 % percentile method
## Scale LFC to a median of negative controls to 0
LFC_nc = subset(LFC, grepl('ONE_INTERGENIC', rownames(LFC), fixed = TRUE)) # log fold change of negative controls
nc_median = apply(LFC_nc, 2, median)
nc_median_df = as.data.frame(matrix(rep(nc_median, each = nrow(LFC)), nrow = nrow(LFC))) # create data frame with median of each cell line to be able to substract the same value from each gene of the same cell line
colnames(nc_median_df) = colnames(LFC)
LFC_scale1 = LFC - nc_median_df # actual step 1 of scaling
LFC_scale1_test = apply(subset(LFC_scale1, grepl('ONE_INTERGENIC', rownames(LFC_scale1), fixed = TRUE)), 2, median)

## Scale LFC to a median of positive controls to -1
LFC_pc = subset(LFC_scale1, rownames(LFC_scale1) %in% pos_cntrl_ADaM)
pc_median = apply(LFC_pc, 2, median)
pc_median_df = as.data.frame(matrix(rep(pc_median, each = nrow(LFC_scale1)), nrow = nrow(LFC_scale1)))
colnames(pc_median_df) = colnames(LFC_scale1)
LFC_scale2 = LFC_scale1/(-pc_median_df) # actual step 2 of scaling
LFC_scale2_test = apply(subset(LFC_scale2, rownames(LFC_scale2) %in% pos_cntrl_ADaM), 2, median)

## Scale LFC on guide level to a median of negative controls to 0, same logic as above
LFC_g_nc = subset(LFC_guide, grepl('ONE_INTERGENIC', rownames(LFC_guide), fixed = TRUE))
nc_g_median = apply(LFC_g_nc[,-(1:2)], 2, median) # first 2 cols have information on genes and guides and are not numbers
nc_g_median_df = as.data.frame(matrix(rep(nc_g_median, each = nrow(LFC_guide)), nrow = nrow(LFC_guide)))
nc_g_median_df = cbind.data.frame(LFC_guide[,(1:2)], nc_g_median_df)
colnames(nc_g_median_df) = colnames(LFC_guide)
LFC_g_scale1 = LFC_guide[,-(1:2)] - nc_g_median_df[,-(1:2)]
LFC_g_scale1_test = apply(subset(LFC_g_scale1, grepl('ONE_INTERGENIC', rownames(LFC_g_scale1),
                                                     fixed = TRUE)), 2, median)
LFC_g_scale1 = cbind.data.frame(LFC_guide[,(1:2)], LFC_g_scale1)
## Scale LFC on guide level to a median of positive controls to -1
LFC_g_pc = subset(LFC_g_scale1, LFC_g_scale1$gene %in% pos_cntrl_ADaM)
pc_g_median = apply(LFC_g_pc[,-(1:2)], 2, median)
pc_g_median_df = as.data.frame(matrix(rep(pc_g_median, each = nrow(LFC_g_scale1)),
                                      nrow = nrow(LFC_g_scale1)))
pc_g_median_df = cbind.data.frame(LFC_guide[,(1:2)], pc_g_median_df)
colnames(pc_g_median_df) = colnames(LFC_g_scale1)
LFC_g_scale2 = LFC_g_scale1[,-(1:2)]/abs(pc_g_median_df[,-(1:2)])
LFC_g_scale2 = cbind.data.frame(LFC_guide[,(1:2)], LFC_g_scale2)
LFC_g_scale2_test = apply(subset(LFC_g_scale2[,-(1:2)],
                                 LFC_g_scale2$gene %in% pos_cntrl_ADaM), 2, median)
LFC_scale_guide = aggregate(LFC_g_scale2[,-(1:2)], list(LFC_g_scale2$gene), mean)
rownames(LFC_scale_guide) = LFC_scale_guide$Group.1
LFC_scale_guide = LFC_scale_guide[,-1]

## Calculations of 90th percentile method first with LFC, then with scaled LFC 
CEG_90 = CoRe.FiPer(LFC, display = TRUE, percentile = 0.9, method = 'fixed')
cairo_pdf('CEG_90th_percentile_Method_Average_Rank.pdf') # method that will be used to identify CEGs, save plot
CEG_avg = CoRe.FiPer(LFC, display = TRUE, percentile = 0.9, method = 'average')
dev.off()
CEG_slope = CoRe.FiPer(LFC, display = TRUE, percentile = 0.9, method = 'slope')
CEG_AUC = CoRe.FiPer(LFC, display = TRUE, percentile = 0.9, method = 'AUC')
consensus_set = intersect(CEG_90$cfgenes, intersect(CEG_AUC$cfgenes, CEG_slope$cfgenes))
# Since the average method identifies exactly the genes in the intersect of all methods,
# another intersect is calculated without the most stringent average method  

## Same with scaled LFC
CEG_90_scaled = CoRe.FiPer(LFC_scale2, display = TRUE, percentile = 0.9, method = 'fixed')
CEG_avg_scaled = CoRe.FiPer(LFC_scale2, display = TRUE, percentile = 0.9, method = 'average')
CEG_slope_scaled = CoRe.FiPer(LFC_scale2, display = TRUE, percentile = 0.9, method = 'slope')
CEG_AUC_scaled = CoRe.FiPer(LFC_scale2, display = TRUE, percentile = 0.9, method = 'AUC')
consensus_set_scaled = intersect(CEG_90_scaled$cfgenes, intersect(CEG_AUC_scaled$cfgenes,
                                                                  CEG_slope_scaled$cfgenes))
## Same with scaled LFC on guide level
CEG_90_scaled_g = CoRe.FiPer(LFC_scale_guide, display = TRUE, percentile = 0.9, method = 'fixed')
CEG_avg_scaled_g = CoRe.FiPer(LFC_scale_guide, display = TRUE, percentile = 0.9, method = 'average')
CEG_slope_scaled_g = CoRe.FiPer(LFC_scale_guide, display = TRUE, percentile = 0.9, method = 'slope')
CEG_AUC_scaled_g = CoRe.FiPer(LFC_scale_guide, display = TRUE, percentile = 0.9, method = 'AUC')
consensus_set_scaled_g = intersect(CEG_90_scaled_g$cfgenes, intersect(CEG_AUC_scaled_g$cfgenes,
                                                                      CEG_slope_scaled_g$cfgenes))

## Determine overlap of CFG and CEG
CFG_CEG = intersect(CFG_other, CEG_avg$cfgenes) # essential genes in this study only if they belong to both CFGs and CEGs
excl_CEG = CEG_avg$cfgenes[!(CEG_avg$cfgenes %in% CFG_CEG) & !(CEG_avg$cfgenes %in% pos_cntrl_ADaM)]

## Benchmarking
## Not enough controls / knwon genes?
# FP = CoRe.AssembleFPs()
# CoRe.CF_Benchmark(CEG_90, background = rownames(LFC), priorKnownSignatures = list(essentials = pos_cntrl_ADaM),
# falsePositives = subset(rownames(LFC), grepl('ONE_INTERGENIC_SITE', rownames(LFC))))

## Define nonessential genes by euclidean distance to 0
euclidean_zero = sqrt(apply((LFC^2), 1, sum))
euclidean_zero_neg = subset(euclidean_zero, grepl('ONE_INTERGENIC_SITE', names(euclidean_zero)))
euclidean_zero_genes = subset(euclidean_zero, !grepl('ONE_INTERGENIC_SITE', names(euclidean_zero)))
euclidean_zero_genes_sorted = euclidean_zero_genes[order(euclidean_zero_genes)]
euclidean_zero_ness = euclidean_zero_genes_sorted[1:5] # five elements with the lowest distance
names(euclidean_zero)[which(euclidean_zero == min (euclidean_zero))] # get name of lowest

## How many nonessentials can be extracted until one is not in the defined set of nonessentials from MAGeCK and BAGEL? 
n_zero = 0
while(all((names(euclidean_zero_genes_sorted)[(1:(n_zero + 1))] %in% nonessentials)) && # until the next entry is possibly an essential or tumor suppressor gene
      (n_zero <= length(euclidean_zero_genes_sorted))) { # break criterion to avoid endless loop
  euclidean_zero_ness_n = euclidean_zero_genes_sorted[(1:(n_zero + 1))]
  n_zero = n_zero + 1
}


## Same with scaled LFC
euclidean_zero_scaled = sqrt(apply((LFC_scale2^2), 1, sum))
euclidean_zero_neg_scaled = subset(euclidean_zero_scaled, grepl('ONE_INTERGENIC_SITE', names(euclidean_zero_scaled)))
euclidean_zero_genes_scaled = subset(euclidean_zero_scaled, !grepl('ONE_INTERGENIC_SITE',
                                                                   names(euclidean_zero_scaled)))
euclidean_zero_genes_scaled_sorted = euclidean_zero_genes_scaled[order(euclidean_zero_genes_scaled)]
euclidean_scaled_ness = euclidean_zero_genes_scaled_sorted[1:5]
names(euclidean_zero_scaled)[which(euclidean_zero_scaled == min(euclidean_zero_scaled))]

## How many nonessentials can be extracted until one is not in the defined set of nonessentials from MAGeCK and BAGEL? 
# euclidean_scaled_ness_n = euclidean_zero_genes_scaled_sorted[1]
n_scaled = 0
while((names(euclidean_zero_genes_scaled_sorted)[(n_scaled + 1)] %in% nonessentials) &&
      (n_scaled < length(euclidean_zero_genes_scaled_sorted))) { 
  euclidean_scaled_ness_n = euclidean_zero_genes_scaled_sorted[(1:(n_scaled + 1))]
  n_scaled = n_scaled + 1
}

## euclidean distance to median of negative control
euclidean_nc = sqrt(apply((LFC-nc_median_df)^2, 1, sum))
euclidean_nc_neg = subset(euclidean_nc, grepl('ONE_INTERGENIC_SITE', names(euclidean_nc)))
euclidean_nc_genes = subset(euclidean_nc, !grepl('ONE_INTERGENIC_SITE', names(euclidean_nc)))
euclidean_nc_genes_sorted = euclidean_nc_genes[order(euclidean_nc_genes)]
euclidean_nc_ness = euclidean_nc_genes_sorted[1:5]
names(euclidean_nc)[which(euclidean_nc == min (euclidean_nc))]

## How many nonessentials can be extracted until one is not in the defined set of nonessentials from MAGeCK and BAGEL? 
n_nc = 0
while((names(euclidean_nc_genes_sorted)[n_nc + 1] %in% nonessentials) &&
      (n_nc < length(euclidean_nc_genes_sorted))) {
  euclidean_nc_ness_n = euclidean_nc_genes_sorted[(1:(n_nc + 1))]
  n_nc = n_nc + 1
}

## Venn diagram of essentials

pdf('Venn_Diagramm_of_Identified_Essentials.pdf')
ggvenn(list(CFG = subset(CFG, !(CFG %in% pos_cntrl_ADaM)), CEG = subset(CEG_avg$cfgenes,
                                                                        !(CEG_avg$cfgenes %in% pos_cntrl_ADaM))),
       fill_color = c("red", "blue"), text_size = 2, auto_scale = TRUE)
dev.off()

## Venn diagram of CEG methods
#first get dataframe of all CEG from different criteria
library("qpcR")
CEG_overview = qpcR:::cbind.na(CEG_avg$cfgenes,CEG_90$cfgenes,CEG_AUC$cfgenes,CEG_slope$cfgenes)
write.csv(CEG_overview, file ="CEG_differing_criteria.csv")
#better subset vectors by "MIR"
library("stringr")  
CEG_avg_mir = CEG_avg$cfgenes[str_detect(CEG_avg$cfgenes, "MIR")]
CEG_fixed_mir = CEG_90$cfgenes[str_detect(CEG_90$cfgenes, "MIR")]
CEG_AUC_mir = CEG_AUC$cfgenes[str_detect(CEG_AUC$cfgenes, "MIR")]
CEG_SLOPE_mir = CEG_slope$cfgenes[str_detect(CEG_slope$cfgenes, "MIR")]

library(ggvenn)
##venn only MIR CEGs 
ggvenn(list(average = CEG_avg_mir, fixed = CEG_fixed_mir, AUC = CEG_AUC_mir,
            slope = CEG_SLOPE_mir))


#venn with all CEGs including essential controls
pdf('Comparison_of_Different_Methods_to_Identify_CEGs_Unscaled_LFC_clean.pdf')
ggvenn(list(average = CEG_avg$cfgenes, fixed = CEG_90$cfgenes, AUC = CEG_AUC$cfgenes, slope = CEG_slope$cfgenes),
       show_percentage = FALSE, text_size = 0)
dev.off()
pdf('Comparison_of_Different_Methods_to_Identify_CEGs_Scaled_LFC.pdf')
ggvenn(list(average = CEG_avg_scaled$cfgenes, fixed = CEG_90_scaled$cfgenes, AUC = CEG_AUC_scaled$cfgenes,
            slope = CEG_slope_scaled$cfgenes), show_percentage = FALSE)
dev.off()
## Different LFCs
ggvenn(list(unscaled = CEG_90$cfgenes, scaled_gene = CEG_90_scaled$cfgenes,
            scaled_guide = CEG_90_scaled_g$cfgenes)) +
  labs(title = 'Fixed')
ggvenn(list(unscaled = CEG_avg$cfgenes, scaled_gene = CEG_avg_scaled$cfgenes,
            scaled_guide = CEG_avg_scaled_g$cfgenes)) +
  labs(title = 'Average')
ggvenn(list(unscaled = CEG_slope$cfgenes, scaled_gene = CEG_slope_scaled$cfgenes,
            scaled_guide = CEG_slope_scaled_g$cfgenes)) +
  labs(title = 'Slope')
ggvenn(list(unscaled = CEG_AUC$cfgenes, scaled_gene = CEG_AUC_scaled$cfgenes,
            scaled_guide = CEG_AUC_scaled_g$cfgenes)) +
  labs(title = 'AUC')
ggvenn(list(unscaled = consensus_set, scaled_gene = consensus_set_scaled,
            scaled_guide = consensus_set_scaled_g)) +
  labs(title = 'Consensus Set')

## Venn diagram of nonessentials
ggvenn(list('FDR' = nonessentials, 'to 0' = names(euclidean_zero_ness), 'to median'
            = names(euclidean_nc_ness), 'scaled' = names(euclidean_scaled_ness)))
ggvenn(list('FDR' = nonessentials, 'to 0' = names(euclidean_zero_ness_n), 'to median'
            = names(euclidean_nc_ness_n), 'scaled' = names(euclidean_scaled_ness_n)))

## Sanity Checks
SanityCheck2_1 = all(CFG == adam)
SanityCheck2_2 = all(LFC_scale1_test == 0) && all.equal(as.numeric(LFC_scale2_test), rep(-1, 45))
SanityCheck2_3 = all.equal(as.numeric(LFC_g_scale1_test), rep(0, 45)) &&
  all.equal(as.numeric(LFC_g_scale2_test), rep(-1, 45))
SanityCheck2_4 = all(CEG_avg$cfgenes %in% consensus_set)
SanityCheck2_5 = all(names(euclidean_zero_ness) %in% nonessentials) &&
  all(names(euclidean_scaled_ness) %in% nonessentials) &&
  all(names(euclidean_nc_ness) %in% nonessentials)

# SanityCheck2_1 # Are the same CFGs determined by the step-by-step and the one-command method?
# SanityCheck2_2 # Are the medians of the scaled LFC equal/nearly equal to 0 (neg controls) / -1 (pos controls)
# SanityCheck2_3 # Are the medians of the scaled LFC on guide level nearly equal to 0 (nc) / -1 (pc)
# SanityCheck2_4 # Are all CEGs identified by the average method included in the consensus set?
# SanityCheck2_5 # Are none identified nonessentials by all three distance methods potential essential or tumor suppressor genes? 
all(SanityCheck2_1,SanityCheck2_2,SanityCheck2_3,SanityCheck2_4,SanityCheck2_5)

## Output
write.xlsx(LFC_scale2, 'Dependency score (scaled LFC).xlsx', rowNames = TRUE)
write.table(CFG, 'All core fitness genes determined by ADaM.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(CFG_other, 'miRNA core fitness genes determined by ADaM.txt', quote = FALSE, row.names = FALSE,
            col.names = FALSE)
write.table(common_essentials, 'Essential genes in all cell lines.txt', quote = FALSE, row.names = FALSE,
            col.names = FALSE)
write.table(crossoverpoint, 'ADaM threshold.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(CEG_avg$cfgenes, 'Core fitness genes determined by FiPer method.txt', quote = FALSE, 
            row.names = FALSE, col.names = FALSE)
write.table(CEG_avg$LocalMinRank, 'Rank threshold FiPer method.txt', quote = FALSE, row.names = FALSE,
            col.names = FALSE)
write.table(CFG_CEG, 'Overlap CFG and CEG.txt', quote = FALSE, row.names = FALSE, col.names = FALSE)
write.xlsx(binmat, 'Binary Matrix ADaM.xlsx', rowNames = TRUE, colNames = TRUE)
write.xlsx(bagel, 'BAGEL2 Output.xlsx')
write.xlsx(bagel_essentials, 'Essentials BAGEL2.xlsx')
write.xlsx(bagel_five, 'Essentials 5% FDR BAGEL2.xlsx')
write.xlsx(lowest_BF, 'Last included essential BAGEL2.xlsx')
write.xlsx(mageck, 'MAGeCK Output.xlsx')
write.xlsx(mageck_essentials, 'Essentials MAGeCK.xlsx')
write.xlsx(mageck_five, 'Essentials 5% FDR MAGeCK.xlsx')
write.xlsx(lowest_mageck, 'Last included essential MAGeCK.xlsx')
write.xlsx(information_essentials, 'BAGEL and MAGeCK Output for Essentials per Cell Line.xlsx')
write.xlsx(information_five, 'BAGEL and MAGeCK Output for Essentials with 5% FDR per Cell Line.xlsx')
write.xlsx(essentials, 'Essential genes per cell line according to BAGEL and MAGeCK.xlsx')
write.xlsx(five, 'Essential genes per cell line according to BAGEL and MAGeCK with 5% FDR.xlsx')
write.xlsx(as.data.frame(n_cntrl), 'Number of positive controls in identified essentials by BAGEL2 and MAGeCK.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(euclidean_nc_ness), 'Nonessential genes (top5) according to euclidean distance to median of negative controls.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(euclidean_zero_ness), 'Nonessential genes (top5) according to euclidean distance to 0.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(euclidean_scaled_ness), 'Nonessential genes (top5) according to euclidean distance of scaled LFC to 0.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(euclidean_zero), 'Euclidean distance of LFC to 0.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(euclidean_nc), 'Euclidean distance of LFC to median of negative controls.xlsx', rowNames = TRUE)
write.xlsx(as.data.frame(euclidean_zero_scaled), 'Euclidean distance of scaled LFC to 0.xlsx', rowNames = TRUE)

## Change function to plot line of 90th percentile
CoRe.VisCFness2<-function(depMat,gene,percentile=0.9,posControl='RPL8',negControl='MAP2K1',method='fixed'){
  gg<-gene
  depMat<-as.matrix(depMat)
  
  rankCL<-t(apply(depMat,1,function(x){
    sx<-order(x)
    x<-match(1:length(x),sx)
  }))
  
  rownames(rankCL)<- rownames(depMat)
  colnames(rankCL)<- colnames(depMat)
  
  rankG<-apply(depMat,2,function(x){
    sx<-order(x)
    x<-match(1:length(x),sx)})
  
  rownames(rankG)<- rownames(depMat)
  colnames(rankG)<- colnames(depMat)
  
  nG<-nrow(rankG)
  nCL<-ncol(rankG)
  
  par(mfrow=c(1,2))
  
  plot(rankG[negControl,names(sort(rankCL[negControl,]))],
       col='transparent',
       pch=16,ylim=c(0,nrow(depMat)),
       xlab='Cell line dependency rank per gene',
       ylab=''#,cex.lab=1.5,
       #cex.axis=1.5
  )
  title(ylab = expression(paste('Gene dependency rank in x' ^'th', ' dependent cell line', sep = '')), line = 2.5)
  
  abline(v=45*percentile,lwd=3,col='darkgrey',lty=2)
  
  points(rankG[negControl,names(sort(rankCL[negControl,]))],pch = 16,#cex=1.5,
         col=rgb(0.2117647, 0.3921569, 0.5450980,alpha = 0.95))
  
  points(rankG[posControl,names(sort(rankCL[posControl,]))],pch=16,#cex=1.5,
         col=rgb(0.8,0,0,alpha = 0.95))
  
  points(rankG[gg,names(sort(rankCL[gg,]))],pch=16,#cex=1.5,
         col=rgb(0.7291,0.6588,0.4549,alpha = 0.95))
  
  threshold = as.integer(nCL*percentile)
  
  if(method=='fixed'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){rankG[x,match(threshold,rankCL[x,])]}))
    Label = 'Gene rank in 90th perc. least dep cell line'
  }
  
  if(method=='average'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){mean(rankG[x,names(which(rankCL[x,]>=threshold))])}))
    Label = 'Gene average rank in ≥ 90th perc. of least dep cell lines'
  }
  
  if(method=='slope'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      b<-as.data.frame(a)
      p<-lm(a ~ seq(1:nCL) , data=b)
      coef(p)[2]
    }))
    Label = 'Slope of gene ranks across ranked dep cell lines'
  }
  
  if(method=='AUC'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      sum(a)
    }))
    Label = 'AUC of gene ranks across ranked dep cell lines'
  }
  
  cc<-c(rgb(0.8,0,0,alpha = 0.95),
        rgb(0.2117647, 0.3921569, 0.5450980,alpha = 0.95),
        rgb(0.7291,0.6588,0.4549,alpha = 0.95))
  legend('topleft',legend=c(posControl,negControl,gg),col=cc,pch=16,bty='n')#,cex=1.5)
  
  names(LeastDependentdf)<-rownames(rankG)
  
  hist(LeastDependentdf,main=Label,xlab='Average rank in least dependent cell lines',ylab='Frequency of genes',
       breaks=40)
  # ,cex.lab=1.5,cex.axis=1.5)
  # polygon(density(LeastDependentdf), col = 'red')
  abline(v = LeastDependentdf[negControl],lwd=4,col=rgb(0.2117647, 0.3921569, 0.5450980,alpha = 0.95))
  abline(v = LeastDependentdf[posControl],lwd=4,col=rgb(0.8,0,0,alpha = 0.95))
  abline(v = LeastDependentdf[gg],lwd=4,col=rgb(0.7291,0.6588,0.4549,alpha = 0.95))
  
  doR <- density(LeastDependentdf, bw = 'nrd0')
  
  localmin <- which(diff(-1*sign(diff(doR$y)))==-2)[1]+1
  myranks<- doR$x
  rankthreshold <- as.integer(myranks[localmin])+1
  
  abline(v = rankthreshold,
         lwd=3,
         col='darkgrey',
         lty=2)
  
}
## Visualise every essential gene
for (ii in 1:length(CFG_CEG)) {
  cairo_pdf(paste('Gene_Rank_', CFG_CEG[ii], '.pdf', sep = ''), 14, 7)
  CoRe.VisCFness2(LFC, CFG_CEG[ii], posControl = pos_cntrl_ADaM[5], negControl = names(euclidean_nc_ness)[1],
                  method = 'average')
  dev.off()
}

## Change function to plot two genes of interest
CoRe.VisCFness3<-function(depMat,gene1,gene2,percentile=0.9,posControl='RPL8',negControl='MAP2K1',method='fixed'){
  gg1<-gene1
  gg2<-gene2
  depMat<-as.matrix(depMat)
  
  rankCL<-t(apply(depMat,1,function(x){
    sx<-order(x)
    x<-match(1:length(x),sx)
  }))
  
  rownames(rankCL)<- rownames(depMat)
  colnames(rankCL)<- colnames(depMat)
  
  rankG<-apply(depMat,2,function(x){
    sx<-order(x)
    x<-match(1:length(x),sx)})
  
  rownames(rankG)<- rownames(depMat)
  colnames(rankG)<- colnames(depMat)
  
  nG<-nrow(rankG)
  nCL<-ncol(rankG)
  
  par(mfrow=c(1,2))
  
  plot(rankG[negControl,names(sort(rankCL[negControl,]))],
       col='transparent',
       pch=16,ylim=c(0,nrow(depMat)),
       xlab='Cell line dependency rank per gene',
       ylab=''#,cex.lab=1.5,
       #cex.axis=1.5
  )
  title(ylab = expression(paste('Gene dependency rank in x' ^'th', ' dependent cell line', sep = '')), line = 2.5)
  
  abline(v=45*percentile,lwd=3,col='darkgrey',lty=2)
  
  points(rankG[negControl,names(sort(rankCL[negControl,]))],pch = 16,#cex=1.5,
         col=rgb(0.2117647, 0.3921569, 0.5450980,alpha = 0.95))
  
  points(rankG[posControl,names(sort(rankCL[posControl,]))],pch=16,#cex=1.5,
         col=rgb(0.8,0,0,alpha = 0.85))
  
  points(rankG[gg1,names(sort(rankCL[gg1,]))],pch=16,#cex=1.5,
         col=rgb(0.7291,0.6588,0.4549,alpha = 0.85))
  
  points(rankG[gg2,names(sort(rankCL[gg2,]))],pch=16,
         col=rgb(0.3, 0.4, 0.35,alpha=0.85))
  
  threshold = as.integer(nCL*percentile)
  
  if(method=='fixed'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){rankG[x,match(threshold,rankCL[x,])]}))
    Label = 'Gene rank in 90th perc. least dep cell line'
  }
  
  if(method=='average'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){mean(rankG[x,names(which(rankCL[x,]>=threshold))])}))
    Label = 'Gene average rank in ≥ 90th perc. of least dep cell lines'
  }
  
  if(method=='slope'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      b<-as.data.frame(a)
      p<-lm(a ~ seq(1:nCL) , data=b)
      coef(p)[2]
    }))
    Label = 'Slope of gene ranks across ranked dep cell lines'
  }
  
  if(method=='AUC'){
    LeastDependentdf<-do.call(rbind,lapply(1:nG,function(x){
      a <- rankG[x,colnames(rankCL)[order(rankCL[x,])]]
      sum(a)
    }))
    Label = 'AUC of gene ranks across ranked dep cell lines'
  }
  
  cc<-c(rgb(0.8,0,0,alpha = 0.85),
        rgb(0.2117647, 0.3921569, 0.5450980,alpha = 0.85),
        rgb(0.7291,0.6588,0.4549,alpha = 0.85),
        rgb(0.3, 0.4, 0.35,alpha=0.85))
  legend('topleft',legend=c(posControl,negControl,gg1,gg2),col=cc,pch=16,bty='n')#,cex=1.5)
  
  names(LeastDependentdf)<-rownames(rankG)
  
  doR <- density(LeastDependentdf, bw = 'nrd0')
  
  localmin <- which(diff(-1*sign(diff(doR$y)))==-2)[1]+1
  myranks<- doR$x
  rankthreshold <- as.integer(myranks[localmin])+1
  
  hist(LeastDependentdf,main=Label,xlab='Average rank in least dependent cell lines',ylab='Frequency of genes',
       breaks=40)
  # ,cex.lab=1.5,cex.axis=1.5)
  axis(1,rankthreshold,lwd=0,cex.axis=0.7,line=-0.5)
  rug(rankthreshold,ticksize=-0.015,side=1)
  
  abline(v = LeastDependentdf[negControl],lwd=4,col=rgb(0.2117647, 0.3921569, 0.5450980,alpha = 0.85))
  abline(v = LeastDependentdf[posControl],lwd=4,col=rgb(0.8,0,0,alpha = 0.85))
  abline(v = LeastDependentdf[gg1],lwd=4,col=rgb(0.7291,0.6588,0.4549,alpha = 0.85))
  abline(v = LeastDependentdf[gg2],lwd=4,col=rgb(0.3, 0.4, 0.35,alpha = 0.85))
  
  
  
  abline(v = rankthreshold,
         lwd=3,
         col='darkgrey',
         lty=2)
  
}

cairo_pdf('Gene_Rank_Examples.pdf', 14, 7)
CoRe.VisCFness3(LFC, gene1 = "MIR708", gene2 = "MIR663A", posControl = pos_cntrl_ADaM[5],
                negControl = names(euclidean_nc_ness)[1], method = 'average')
dev.off()


## Ridge plot of some essentials and nonessentials

ridge = as.data.table(cbind(Gene = rownames(LFC_scale2), LFC_scale2))
essentials10 = c("MIR663A","MIR6859-3",'MIR8078',"MIR6859-2", 'MIR4448', "MIR3134",'MIR483', "MIR622",'MIR4454', 'MIR151A')
ridge_mean = data.frame(Gene = 'Mean', Cell_Line = names(apply(LFC_scale2, 1, mean)),
                        Dependency_Score = apply(LFC_scale2, 1, mean), category = 'mean')
colnames(ridge_mean)[c(2,3)] = c('Cell Line', 'Dependency Score')
one_inter_ridge = subset(ridge, grepl('ONE_INTERGENIC', Gene))
mean_one_inter = cbind('ONE_INTERGENIC',names(apply(one_inter_ridge[,-1], 2, mean)),
                       apply(one_inter_ridge[,-1], 2, mean))
colnames(mean_one_inter) = c('Gene', 'Cell Line', 'Dependency Score')
gene_ridge = subset(ridge, !grepl('ONE_INTERGENIC', Gene))
gene_long = melt(gene_ridge, id.vars = 'Gene', variable.name = 'Cell Line')
colnames(gene_long)[3] = 'Dependency Score'
ridge_long = rbind(gene_long, mean_one_inter)
ridge_neg = subset(ridge_long, grepl('ONE_INTERGENIC', Gene))
ridge_neg = cbind(ridge_neg, category = 'negative control')
ridge_pos = subset(ridge_long, Gene %in% pos_cntrl_ADaM[c(5,3,1)])
# ridge_pos = ridge_pos[order(ridge_pos$Gene),]
ridge_pos = cbind(ridge_pos, category = 'positive control')
ridge_ness = subset(ridge_long, Gene %in% names(euclidean_nc_ness)[1:3])
# ridge_ness = ridge_ness[order(ridge_ness$Gene),]
ridge_ness = cbind(ridge_ness, category = 'nonessentials')
ridge_ess = subset(ridge_long, Gene %in% essentials10)
# ridge_ess = ridge_ess[order(ridge_ess$Gene),]
ridge_ess = cbind(ridge_ess, category = 'essentials')
ridge_plot = rbind(ridge_neg, ridge_pos, ridge_ness, ridge_ess, ridge_mean)
ridge_plot$`Dependency Score` = as.numeric(ridge_plot$`Dependency Score`)
gorder = c(essentials10, names(euclidean_nc_ness)[1:3],pos_cntrl_ADaM[c(5,3,1)], 'ONE_INTERGENIC', 'Mean')
# ggplot(ridge_plot, aes(x = `Dependency Score`, y = `Gene`, fill = stat(x))) +
# geom_density_ridges_gradient(jittered_points = TRUE,
# position = position_points_jitter(width = 0.03, height = 0),
# point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.7) +
# scale_y_discrete(limits = unique(ridge_plot$Gene)) +
# scale_fill_viridis_c(name = 'DS', option = "F") +
# coord_cartesian(clip = "off") +
# geom_vline(xintercept = 0, linetype = "dashed", lwd = 1.2) +
# geom_vline(xintercept = -1, linetype = 'dashed', lwd = 1.2)+
# theme_ridges(grid= FALSE)
pdf('Ridgeplot_Selected_Genes.pdf')
ggplot(ridge_plot, aes(x = `Dependency Score`, y = `Gene`, fill = stat(x))) +
  geom_density_ridges_gradient(jittered_points = TRUE,
                               position = position_points_jitter(width = 0.03, height = 0),
                               point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.7) +
  scale_y_discrete(limits = gorder) +
  scale_fill_gradientn(name = 'DS', colours = rgb(rbind(c(1, 1, 0), c(0.8, 0,0),
                                                        c(0.2117647, 0.3921569, 0.5450980), c(0, 0, 0.4))),
                       values =
                         c(0, ((-1.1-min(ridge_plot$`Dependency Score`))/
                                 (max(ridge_plot$`Dependency Score`)-min(ridge_plot$`Dependency Score`))),
                           ((-0.1-min(ridge_plot$`Dependency Score`))/
                              (max(ridge_plot$`Dependency Score`)-min(ridge_plot$`Dependency Score`))), 1)) +
  coord_cartesian(clip = "off") +
  geom_vline(xintercept = 0, linetype = "dashed", color = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 1.2) +
  geom_vline(xintercept = -1, linetype = 'dashed', color = rgb(0.8, 0, 0), lwd = 1.2)+
  theme_ridges(grid= FALSE, center_axis_labels = TRUE)
dev.off()


## Ridge plot with distribution of all and all essentials for supplement

ridge_pos_all = subset(ridge_long, Gene %in% pos_cntrl_ADaM)
ridge_pos_all = cbind(ridge_pos_all, category = 'positive control')
ridge_ess_all = subset(ridge_long, Gene %in% CFG_CEG)
ridge_ess_all = cbind(ridge_ess_all, category = 'essentials')
ridge_ness_all = subset(ridge_long, Gene %in% names(euclidean_nc_ness))
ridge_ness_all = cbind(ridge_ness_all, category = 'nonessentials')
ridge_plot_all = rbind(ridge_neg, ridge_pos_all, ridge_ness_all, ridge_ess_all, ridge_mean)
ridge_plot_all$`Dependency Score` = as.numeric(ridge_plot_all$`Dependency Score`)
gorder_all = c(CFG_CEG[order(CFG_CEG, decreasing = TRUE)], names(euclidean_nc_ness), pos_cntrl_ADaM,
               'ONE_INTERGENIC', 'Mean')

pdf('Ridgeplot_All_Essential_Genes.pdf', 15, 40)
ggplot(ridge_plot_all, aes(x = `Dependency Score`, y = `Gene`, fill = stat(x))) +
  geom_density_ridges_gradient(jittered_points = TRUE,
                               position = position_points_jitter(width = 0.03, height = 0),
                               point_shape = '|', point_size = 2, point_alpha = 1, alpha = 0.7) +
  scale_y_discrete(limits = gorder_all) +
  scale_fill_gradientn(name = 'DS', colours = rgb(rbind(c(1, 1, 0.4), c(0.8, 0,0),
                                                        c(0.2117647, 0.3921569, 0.5450980), c(0, 0, 0.4))),
                       values =
                         c(0, ((-1.1-min(ridge_plot_all$`Dependency Score`))/
                                 (max(ridge_plot_all$`Dependency Score`)-min(ridge_plot_all$`Dependency Score`))),
                           ((-0.2
                             -min(ridge_plot_all$`Dependency Score`))/
                              (max(ridge_plot_all$`Dependency Score`)-min(ridge_plot_all$`Dependency Score`))), 1)) +
  coord_cartesian(clip = "off") +
  geom_vline(xintercept = 0, linetype = "dashed", color = rgb(0.2117647, 0.3921569, 0.5450980), lwd = 1.2) +
  geom_vline(xintercept = -1, linetype = 'dashed', color = rgb(0.8, 0, 0), lwd = 1.2)+
  theme_ridges(grid= FALSE, center_axis_labels = TRUE)
dev.off()

## Plot controls LFC
LFC_posc = subset(LFC, rownames(LFC) %in% pos_cntrl$V1)
LFC_negc = subset(LFC, grepl('ONE_INTER', rownames(LFC)))
mpc = apply(LFC_posc, 2, mean) 
mnc = apply(LFC_negc, 2, mean)
LFC_ess = subset(LFC, rownames(LFC) %in% CFG_CEG)
LFC_ness = subset(LFC, rownames(LFC) %in% names(euclidean_nc_ness))
mess = apply(LFC_ess, 2, mean)
mness = apply(LFC_ness, 2, mean)
plot_ess_ness = data.frame(p_cntrl = mpc, n_cntrl = mnc, essentials = mess, nonessentials = mness)
plot_ess_ness_order = plot_ess_ness[order(plot_ess_ness$p_cntrl),]
plot_essentiality = cbind(1:nrow(plot_ess_ness_order), plot_ess_ness_order)
pdf('LFC_Essential_and_Nonessential_Genes_in_Cell_lines.pdf')
plot(plot_essentiality$essentials, pch = 19, col = 'grey', bty = 'n', xlim = c(1, 50),
     ylim = c(min(plot_essentiality$p_cntrl, plot_essentiality$n_cntrl, plot_essentiality$essentials,
                  plot_essentiality$nonessentials), 
              max(plot_essentiality$essentials, plot_essentiality$nonessentials, plot_essentiality$n_cntrl,
                  plot_essentiality$p_cntrl)), xlab = 'Cell Lines',
     ylab = 'Mean LFC')
abline(h = 0, col = 'black')
points(plot_essentiality$p_cntrl, pch = 20, col = rgb(0.8, 0, 0))
points(plot_essentiality$nonessentials, pch = 19, col = rgb(0.7291, 0.6588, 0.4549))
points(plot_essentiality$n_cntrl, pch = 20, col = rgb(0.2117647, 0.3921569, 0.5450980))
legend('bottomright', c('One intergenic', 'Positive control', 'Identified essentials',
                        'Identified nonessentials'), bty = 'n', pch = c(20, 20, 19, 19),
       col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.8, 0, 0), c(col2rgb('grey')/225),
                       c(0.7291, 0.6588, 0.4549))))
dev.off()

## How many identified essentials are essential in less than 50 % of cell lines?
essential_in_cl = apply(binmat, 1, sum) # Number of cell lines in which this gene is essential
idessential_in_cl = subset(essential_in_cl, names(essential_in_cl) %in% CFG_CEG)
CEG_ordered_all = rownames(CEG_avg$geneRanks)[order(CEG_avg$geneRanks)] # order CEG according to their rank
CEG_ordered = subset(CEG_ordered_all, CEG_ordered_all %in% CFG_CEG)
CEG_ordered[which(CEG_ordered %in% CFG_CEG)] == CFG_other[which(CFG_other %in% CFG_CEG)]
CFG_ordered = subset(CFG_other, CFG_other %in% CFG_CEG)
order_idessential = idessential_in_cl[order(idessential_in_cl)]
pdf('Number_of_dependent_cell_lines_per_identified_essential_miRNA.pdf')
par(mar = c(5.8, 4.1, 4.1, 2.1))
plot(order_idessential, type = 'h', ylim = c(0,45), xlim = c(1,33), ylab = 'Number of dependent cell lines', bty = 'n',
     yaxt = 'n', xaxt = 'n', xlab = '')
axis(1, at = 1:33, labels = names(order_idessential), las = 3, tck = 0, pos = 0)
axis(2, at = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45), labels = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45), las = 2) 
abline(h = 0.75*45, col = rgb(1, 0, 0, alpha = 0.5), lwd = 2)
abline(h = 0.5*45, col = rgb(0.8, 0, 0, alpha = 0.8), lwd = 2)
abline(h = 0.25*45, col = rgb(0.5, 0, 0), lwd = 2)
abline(h = 45, col = 'grey', lwd = 2)
points(order_idessential, pch = 20)
dev.off()

## Number of dependent cell lines for every essential in at least one cell line
one_essential = subset(essential_in_cl, essential_in_cl != 0) # only the genes essential in at least one cell line
one_essential = one_essential[order(one_essential, decreasing = TRUE)]
one_essential_miRNA = subset(one_essential, !(names(one_essential) %in% pos_cntrl_ADaM))
noness = subset(one_essential_miRNA, !(names(one_essential_miRNA) %in% CFG_CEG))
idess = subset(one_essential_miRNA, names(one_essential_miRNA) %in% CFG_CEG)
above = subset(one_essential_miRNA, one_essential_miRNA >= 45/2)
pdf('Number_of_Dependent_Cell_Lines_for_Genes_Essential_in_at_Least_1_Cell_line.pdf')
plot(one_essential_miRNA, type = 'S',  bty = 'n', xlab = 'Number of genes', yaxt = 'n', xaxt = 'n', col = 'transparent',
     ylab = 'Number of dependent cell lines')
axis(1, at = c(1, length(above),50, 100, 150, 200, 230), labels = c(1, length(above),50, 100, 150, 200, 230),
     tck = 0)
rug(c(length(above), 230), tck = -0.01)
rug(c(1, 50, 100, 150, 200), tck = -0.03)
axis(2, at = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45), labels = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45), las = 2)
abline(v = (length(above)+0.5), col = 'grey', lty = 3, lwd = 3)
abline(h = crossoverpoint, col = rgb(0.7291, 0.6588, 0.4549), lwd = 3)
points(which(names(one_essential_miRNA) %in% names(noness)), noness, pch = 20,
       col = rgb(0.2117647, 0.3921569, 0.5450980))
points(which(names(one_essential_miRNA) %in% names(idess)), idess, pch = 20, col = rgb(0.8, 0, 0))
segments(-9, 45/2, (length(above)+0.5), 45/2, col = 'grey', lwd = 3)
legend('topright', c('Essential in at least one cell line', 'Identified fitness gene', '50 % of cell lines',
                     'ADaM threshold'), lty = c(NA,NA, 1, 1), pch = c(20, 20, NA, NA), lwd = 5, bty = 'n',
       col = rgb(rbind(c(0.2117647, 0.3921569, 0.5450980), c(0.8, 0, 0), (col2rgb('grey')/225)[,1],
                       c(0.7291, 0.6588, 0.4549))))
dev.off()

## Rebuild plots profile of fitness genes without positive controls
one_essential_woc = subset(one_essential, (!names(one_essential) %in% pos_cntrl_ADaM))
missing = rep(0, length(which(!(1:45 %in% names(table(one_essential_woc)))))) # fill the n of cell lines for which no gene is essential with 0
names(missing) = which(!(1:45 %in% names(table(one_essential_woc))))
ess_in_line = c(table(one_essential_woc), missing)
ess_in_line = ess_in_line[order(as.numeric(names(ess_in_line)))]
pdf('Essential_Genes_in_Exactly_n_Cell_Lines.pdf')
barplot(ess_in_line, col = rgb(0.2117647, 0.3921569, 0.5450980), xaxt = 'n', space = 0,
        ylab = 'Number of genes essential in n cell lines', xlab = 'n',
        main = paste(length(one_essential_woc), ' essential genes in at least 1 cell line', sep = ''))
axis(1, at = (c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)-0.5), labels = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45), pos = 0)
dev.off()

minima = density(ess_in_line)$x[which(diff(sign(diff(density(ess_in_line)$y))) == 2)]
write.xlsx(as.data.frame(minima), 'Local minima of density essential genes in exactly n cell lines for context-dependent essentials.xlsx')
x_plot = c((density(ess_in_line)$x)[1]-((density(ess_in_line)$x)[2]-(density(ess_in_line)$x)[1]),
           density(ess_in_line)$x)
y_plot = c(0, density(ess_in_line)$y)
pdf('Density_Essential_Genes_in_Exactly_n_Cell_Lines_with_first_Minimum_for_Conext_Dependent_Essentials.pdf')
plot(density(ess_in_line), bty = 'n', col = 'transparent', main = '', xlab = 'n',
     ylab = 'Density of genes essential in n cell lines')
polygon(x_plot, y_plot, col = rgb(0.2117647, 0.3921569, 0.5450980), border = 'transparent')
abline(v = minima[1], col = rgb(0.8, 0, 0), lwd = 2)
dev.off()
ess_at_least_line = ess_in_line
for (ii in as.numeric(names(ess_at_least_line))) {
  ess_at_least_line[ii] = sum(ess_at_least_line[ii:length(ess_at_least_line)])
}
perc_at_least = ess_at_least_line/length(one_essential_woc)
yaxis = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
pdf('Percentage_of_Genes_Essential_in_at_Least_n_Cell_Lines_With_Percentages_Printed.pdf')
par(xpd = TRUE)
barplot(perc_at_least, col = rgb(0.2117647, 0.3921569, 0.5450980), xaxt = 'n', space = 0, axes = FALSE,
        ylab = 'Percentage of genes essential in at least n cell lines', xlab = 'n',
        main = paste(length(one_essential_woc), ' essential genes in at least 1 cell line', sep = ''))
axis(1, at = (c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)-0.5), labels = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45), pos = 0)
axis(2, at = yaxis/100, labels = paste(yaxis, ' %', sep = ''), las = 2, pos = -1)
text(x = ((1:45)-0.47), y = perc_at_least+0.015, round(perc_at_least*100, digits = 0), adj = 0.5, cex = 0.6)
dev.off()

pdf('Percentage_of_Genes_Essential_in_at_Least_n_Cell_Lines.pdf')
barplot(perc_at_least, col = rgb(0.2117647, 0.3921569, 0.5450980), xaxt = 'n', space = 0, axes = FALSE,
        ylab = 'Percentage of genes essential in at least n cell lines', xlab = 'n',
        main = paste(length(one_essential_woc), ' essential genes in at least 1 cell line', sep = ''))
axis(1, at = (c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45)-0.5), labels = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45), pos = 0)
axis(2, at = yaxis/100, labels = paste(yaxis, ' %', sep = ''), las = 2, pos = -1)
dev.off()

## Vulcano plot 
## get p-values from every gene in a cell line, take the lowest from positive and negative p value
LFC_woc = subset(LFC, !(rownames(LFC) %in% pos_cntrl_ADaM | grepl('NO_SITE', rownames(LFC)) |
                          grepl('ONE_INTER', rownames(LFC))))
for (ii in 1:length(mageck)) {
  mageck_woc = subset(mageck[[ii]], !(id %in% pos_cntrl_ADaM | grepl('NO_SITE', id) | grepl('ONE_INTER', id)))
  mageck_woc = mageck_woc[order(mageck_woc$id),]
  p_val = data.frame(id = mageck_woc$id, p_neg = mageck_woc$neg.p.value, p_pos = mageck_woc$pos.p.value, 
                     p_value = apply(cbind(mageck_woc$neg.p.value, mageck_woc$pos.p.value), 1, min),
                     LFC = LFC_woc[,ii][order(rownames(LFC_woc))])
  p_val_order = p_val[order(p_val$LFC),]
  p_val_order_ess = subset(p_val_order, p_neg < 0.01)
  p_val_order_ts = subset(p_val_order, p_pos < 0.01)
  pdf(paste('Vulcano_Plot_LFC_p_value_MAGeCK_', names(mageck)[ii], '.pdf', sep = ''))
  plot(p_val_order$LFC, p_val_order$p_value, pch = 20, ylim = c(1,0), col = 'grey', xlab = 'log fold change',
       ylab = 'p value according to MAGeCK', bty = 'n', main = names(mageck)[ii])
  points(p_val_order_ess$LFC, p_val_order_ess$p_value, pch = 20, col = rgb(0.8, 0, 0))
  points(p_val_order_ts$LFC, p_val_order_ts$p_value, pch = 20, col = rgb(0.2117647, 0.3921569, 0.5450980))
  dev.off()
}

## Build heatmap
LFC_nonessentials = subset(LFC_scale2, rownames(LFC_scale2) %in% names(euclidean_nc_ness))
LFC_CFG_CEG = subset(LFC_scale2, rownames(LFC_scale2) %in% CFG_CEG)
LFC_essentials = subset(LFC_scale2, rownames(LFC_scale2) %in% pos_cntrl_ADaM)

LFC_other = subset(LFC_scale2, !(rownames(LFC_scale2) %in% names(euclidean_nc_ness) |
                                   rownames(LFC_scale2) %in% CFG_CEG | rownames(LFC_scale2) %in% pos_cntrl_ADaM))
LFC_other_v = cbind(LFC_other, Variance = apply(LFC_other, 1, var))
LFC_other_v = LFC_other_v[order(LFC_other_v$Variance, decreasing = TRUE),]
LFC_other_hv = LFC_other_v[1:30,-ncol(LFC_other_v)]

LFC_v = cbind(LFC_scale2, variance = apply(LFC, 1, var))
LFC_v = LFC_v[order(LFC_v$variance, decreasing = TRUE),]
LFC_hv = LFC_v[1:30,-ncol(LFC_v)]

color1 = colorRampPalette(rgb(rbind(c(1, 1, 0.4), c(0.8,0,0), c(0.2117647, 0.3921569, 0.5450980))), bias = 0.8)(1500)
color2 = colorRampPalette(rgb(rbind(c(1, 0.8, 0), c(0.8,0,0), c(0.2117647, 0.3921569, 0.5450980),
                                    c(0.4, 1, 1))), bias = 1.5)(1500)

## Heatmap A: Overlap of CFEGs and CEGs with scaled LFC
pdf('Heatmap_of_Identified_Essentials.pdf', 10, 10)
pheatmapA = pheatmap(rbind(LFC_nonessentials, LFC_CFG_CEG, LFC_essentials), cluster_row = FALSE,
                     color = color1, gaps_row = c(nrow(LFC_nonessentials), nrow(LFC_nonessentials)+nrow(LFC_CFG_CEG)),
                     border_color = NA, breaks = NULL)
grob_classesA = purrr::map(pheatmapA$gtable$grobs, class)
idx_grobA = which(purrr::map_lgl(grob_classesA, function(cl) 'gTree' %in% cl))[1]
grob_namesA = names(pheatmapA$gtable$grobs[[idx_grobA]]$children)
idx_rectA = grob_namesA[grep('rect', grob_namesA)][1]
pheatmapA$gtable$grobs[[idx_grobA]]$children[[idx_rectA]]$gp$col =
  pheatmapA$gtable$grobs[[idx_grobA]]$children[[idx_rectA]]$gp$fill
pheatmapA$gtable$grobs[[idx_grobA]]$children[[idx_rectA]]$gp$lwd = 3
pheatmapA
dev.off()

## heatmap B: Highly variable genes
pdf('Heatmap_of_Highly_Variable_Genes.pdf', 10, 10)
pheatmapO = pheatmap(rbind(LFC_nonessentials, LFC_other_hv, LFC_essentials), cluster_row = FALSE,
                     color = color2, gaps_row = c(nrow(LFC_nonessentials), nrow(LFC_nonessentials)+nrow(LFC_other_hv)),
                     border_color = NA)
grob_classesO = purrr::map(pheatmapO$gtable$grobs, class)
idx_grobO = which(purrr::map_lgl(grob_classesO, function(cl) 'gTree' %in% cl))[1]
grob_namesO = names(pheatmapO$gtable$grobs[[idx_grobO]]$children)
idx_rectO = grob_namesO[grep('rect', grob_namesO)][1]
pheatmapO$gtable$grobs[[idx_grobO]]$children[[idx_rectO]]$gp$col =
  pheatmapO$gtable$grobs[[idx_grobO]]$children[[idx_rectO]]$gp$fill
pheatmapO$gtable$grobs[[idx_grobO]]$children[[idx_rectO]]$gp$lwd = 3
pheatmapO
dev.off()


#check confounding factors in dataset

rm(list = ls())
library(openxlsx)
library(stringr)
library(Cairo)
library(psych)
library(data.table)

colors_lin = read.xlsx('directory/Color coding lineage.xlsx')
cl_info = read.xlsx("directory/B41_miRNA_screen_cell_line_overview_stripped.xlsx")
LFC = read.xlsx("directory/Log fold changes on gene and cell line level.xlsx", rowNames = TRUE)
pos_cntrl = read.csv('directory/Control_essentials.csv', header = FALSE)
correlation = read.xlsx("directory/Intra-line correlations all genes corrected LFC.xlsx", rowNames = TRUE)
inline = read.xlsx("directory/Inline_Assays_stripped.xlsx")
run = read.xlsx("directory/Information on runs good cell lines.xlsx", rowNames = TRUE)

## plot with inner 50% quantile for controls
LFC_posc = subset(LFC, rownames(LFC) %in% pos_cntrl$V1)
LFC_negc = subset(LFC, grepl('ONE_INTER', rownames(LFC)))
qpc = apply(LFC_posc, 2, FUN = 'quantile', probs = c(0.25, 0.5, 0.75))
qnc = apply(LFC_negc, 2, FUN = 'quantile', probs = c(0.25, 0.5, 0.75))
rownames(qpc) = paste(rownames(qpc), 'positive')
rownames(qnc) = paste(rownames(qnc), 'negative')
quantiles = rbind(qpc, qnc)
quantiles = quantiles[,order(quantiles[2,])]
pdf('Mean_LFC_of_controls.pdf')
plot(quantiles[2,], ylim = c(min(quantiles), max(quantiles)), col = 'transparent', bty = 'n', xaxt = 'n',
     xlab = 'Cell Lines', ylab = 'Mean LFC')
axis(1, at = c(1, 5, 10, 15, 20, 25, 30, 35, 40, 45))
# segments(1:45, quantiles[1,], 1:45, quantiles[3,], col = rgb(0.8, 0, 0, alpha = 0.5), lwd = 10)
# polygon(c(1:45, 45:1), c(quantiles[1,],  rev(quantiles[3,])), col = rgb(0.8, 0, 0, alpha = 0.5),
# border = 'transparent')
for (ii in 1:45) {
  polygon(c(ii-0.5, ii+0.5, ii+0.5, ii - 0.5), c(quantiles[1,ii], quantiles[1,ii], quantiles[3,ii], quantiles[3,ii]),
          col = rgb(0.8, 0, 0, alpha = 0.5), border = 'transparent')
}
points(quantiles[2,], col = rgb(0.8, 0, 0), pch = 20)
# segments(1:45, quantiles[4,], 1:45, quantiles[6,], col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = 0.5), lwd = 10)
# polygon(c(1:45, 45:1), c(quantiles[4,],  rev(quantiles[6,])), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = 0.5),
# border = 'transparent')
for (ii in 1:45) {
  polygon(c(ii-0.5, ii+0.5, ii+0.5, ii - 0.5), c(quantiles[4,ii], quantiles[4,ii], quantiles[6,ii], quantiles[6,ii]),
          col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = 0.5), border = 'transparent')
}
points(quantiles[5,], col = rgb(0.2117647, 0.3921569, 0.5450980), pch = 20)
abline(h = median(quantiles[2,]), col = rgb(0.8, 0, 0, alpha = 0.8), lwd = 3)
abline(h = median(quantiles[5,]), col = rgb(0.2117647, 0.3921569, 0.5450980, alpha = 0.8), lwd = 3)
legend('bottomright', c('Median', 'Inner 50 % quantile', 'Grand median', '', 'Negative controls', 
                        'Positive controls'), col = rgb(rbind(t(replicate(4, (col2rgb('grey')/255)[,1])),
                                                              c(0.2117647, 0.3921569, 0.5450980), c(0.8, 0, 0))),
       lty = c(rep(NA, 2), 1, rep(NA, 3)), pch = c(20, 15, NA, NA, 18, 18), lwd = 3, bty = 'n',
       pt.cex = c(1, 1.8, 1, 1, 2.2, 2.2))
dev.off()

## Confounding factors

## Determine cell lines with a median LFC of positive controls above or below the grand median
below = names(which(quantiles[2,] < median(quantiles[2,])))
above = names(which(quantiles[2,] >= median(quantiles[2,])))

## intraline correlation 
mcor = apply(correlation, 1, mean, na.rm = TRUE)
## Add column with information if median LFC of positive controls is above the grand median 
mcor_ab = data.frame(mcor = mcor, above_med = as.numeric(names(mcor) %in% above),
                     median = quantiles[2,][order(names(quantiles[2,]))])

## Prepare colors
## Format the names of lineages so that they match
cl_info$Lineage = gsub('_', ' ', cl_info$Lineage)
cl_info$Lineage = str_to_sentence(cl_info$Lineage)
cl_info$Lineage = gsub("Central nervous system", 'CNS', cl_info$Lineage)
cl_info$Lineage = gsub("Peripheral nervous system", 'PNS', cl_info$Lineage)
colnames(colors_lin)[1] = 'Lineage'
colors_lin$Lineage = gsub('Duct', 'duct', colors_lin$Lineage)

## Assign colors to cell lines according to their lineage
colors_cl = data.frame(Cell_line = cl_info$Cell.line, Lineage = cl_info$Lineage, inside = NA, outside = NA,
                       above_med = mcor_ab$above_med)
for (ii in 1:nrow(colors_cl)) {
  colors_cl[ii,(3:4)] = colors_lin[which(toupper(colors_lin$Lineage) == toupper(colors_cl$Lineage[ii])),2:3]
}

## Make legend to insert in a larger figure if needed
pdf('Legend_for_Lineage.pdf')
par(xpd = TRUE)
plot(1, axes = FALSE, bty = 'n', col = 'transparent', xlab = '', ylab = '')
legend(0.445, 1.6, colors_lin$Lineage[1:8], pch = 21, col = colors_lin$outside[1:8], pt.bg = colors_lin$inside[1:8],
       pt.lwd = 3, pt.cex = 4.2, cex = 3, bty = 'n')
legend(0.945, 1.6, colors_lin$Lineage[9:16], pch = 21, col = colors_lin$outside[9:16], pt.bg = colors_lin$inside[9:16],
       pt.lwd = 3, pt.cex = 4.2, cex = 3, bty = 'n')
dev.off()

## statistical tests
t_cor = t.test(subset(mcor_ab$mcor, mcor_ab$above_med == 0), subset(mcor_ab$mcor, mcor_ab$above_med == 1))
r_cor = cor.test(mcor_ab$median, mcor_ab$mcor)

## Generate plot with boxplots of the two sets (above and below grand median) on the left and scatterplot to
## illustrate a possible linear correlation between median LFC of positive controls and the mean correlation
## coefficient of replicates in a cell line
cairo_pdf('Confounding_Factor_Correlation.pdf', 20, 10)
par(mfrow=c(1,2))
boxplot(mcor ~ above_med, data = mcor_ab, col = 'transparent', axes = FALSE, 
        main = substitute(paste(italic(p), ' = ', x, sep = ''), list(x = format(round(t_cor$p.value, 4), nsmall = 4))),
        ylab = substitute(paste('Average pair-wise replicate LFC ', italic(r), sep = '')),
        xlab = 'Median LFC of positive controls')
axis(1, at = c(1, 2), labels = c('< grand median', '\u2265 grand median'), lwd = 0)
axis(2)
for (ii in 1:length(unique(colors_cl$Lineage))) {
  plotcol = subset(colors_cl, Lineage == unique(colors_cl$Lineage)[ii])
  plotcor = subset(mcor_ab, rownames(mcor_ab) %in% plotcol$Cell_line)
  if(all(plotcor$above_med == 1)) {
    plotcor = rbind(plotcor, dummy = c(0,0,0))
    stripchart(mcor ~ above_med, data = plotcor, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE, lwd = 3,
               add = TRUE, col = c('transparent', plotcol$outside[1]), bg = plotcol$inside[1], cex = 1.4, jitter = 0.3)
  } else {
    stripchart(mcor ~ above_med, data = plotcor, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE, lwd = 3,
               add = TRUE, col = plotcol$outside[1], bg = plotcol$inside[1], cex = 1.4, jitter = 0.3)
  }
}
plot(mcor_ab$median, mcor_ab$mcor, pch = 21, col = colors_cl$outside, bg = colors_cl$inside, bty = 'n', lwd = 3,
     main = substitute(paste(italic(r), ' = ', x, ', ', italic(p), ' = ', y, sep = ''),
                       list(x = format(round(r_cor$estimate, 4), nsmall = 4),
                            y = format(round(r_cor$p.value, 4), nsmall = 4))), yaxt = 'n', ylab = '', cex = 1.4,
     xlab = 'Median LFC of positive controls')
dev.off()
cairo_pdf('Confounding_Factor_Correlation_BIG.pdf', 20, 10)
par(mfrow=c(1,2))
boxplot(mcor ~ above_med, data = mcor_ab, col = 'transparent', axes = FALSE, cex.main = 1.5, cex.sub = 1.5,
        main = substitute(paste(italic(p), ' = ', x, sep = ''), list(x = format(round(t_cor$p.value, 4), nsmall = 4))),
        ylab = '', cex.lab = 1.5, cex.axis = 1.5,
        xlab = 'Median LFC of positive controls')
axis(1, at = c(1, 2), labels = c('< grand median', '\u2265 grand median'), lwd = 0, cex.axis = 1.5)
axis(2, cex.axis = 1.5, line = -1.5)
title(ylab = substitute(paste('Average pair-wise replicate LFC ', italic(r), sep = '')), line = 1.3,
      cex.lab = 1.5)
for (ii in 1:length(unique(colors_cl$Lineage))) {
  plotcol = subset(colors_cl, Lineage == unique(colors_cl$Lineage)[ii])
  plotcor = subset(mcor_ab, rownames(mcor_ab) %in% plotcol$Cell_line)
  if(all(plotcor$above_med == 1)) {
    plotcor = rbind(plotcor, dummy = c(0,0,0))
    stripchart(mcor ~ above_med, data = plotcor, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE, lwd = 3,
               add = TRUE, col = c('transparent', plotcol$outside[1]), bg = plotcol$inside[1], cex = 2.2, jitter = 0.3)
  } else {
    stripchart(mcor ~ above_med, data = plotcor, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE, lwd = 3,
               add = TRUE, col = plotcol$outside[1], bg = plotcol$inside[1], cex = 2.2, jitter = 0.3)
  }
}
plot(mcor_ab$median, mcor_ab$mcor, pch = 21, col = colors_cl$outside, bg = colors_cl$inside, bty = 'n', lwd = 3,
     main = substitute(paste(italic(r), ' = ', x, ', ', italic(p), ' = ', y, sep = ''),
                       list(x = format(round(r_cor$estimate, 4), nsmall = 4),
                            y = format(round(r_cor$p.value, 4), nsmall = 4))), yaxt = 'n', ylab = '', cex = 2.2,
     xlab = 'Median LFC of positive controls', cex.main = 1.5, cex.sub = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

## Transduction efficiency
inline = inline[order(inline$Cell.line),] # order to fit with other tables
## Build data frame with transduction efficiency and information on whether the median LFC of positive controls is
## the grand median
efficiency = data.frame(Cell_line = inline$Cell.line, efficiency = inline$Transduction.as.determined.by.inline.assay,
                        median = quantiles[2,][order(names(quantiles[2,]))])
efficiency = cbind.data.frame(efficiency, above_med = as.numeric(efficiency$Cell_line %in% above))
## Statistical tests
t_eff = t.test(subset(efficiency$efficiency, efficiency$above_med == 0),
               subset(efficiency$efficiency, efficiency$above_med == 1))
r_eff = cor.test(efficiency$median, efficiency$efficiency)

## Generate plot with boxplots of the two sets (above and below grand median) on the left and scatterplot to
## illustrate a possible linear correlation between median LFC of positive controls and the transduction efficiency
## for the library in a cell line
cairo_pdf('Confounding_Factor_Transduction_Efficiency.pdf', 20, 10)
par(mfrow=c(1,2))
boxplot(efficiency ~ above_med, data = efficiency, col = 'transparent', axes = FALSE, 
        main = substitute(paste(italic(p), ' = ', x, sep = ''), list(x = format(round(t_eff$p.value, 4), nsmall = 4))),
        ylab = 'Average transduction efficiency',
        xlab = 'Median LFC of positive controls')
axis(1, at = c(1, 2), labels = c('< grand median', '\u2265 grand median'), lwd = 0)
axis(2)
for (ii in 1:length(unique(colors_cl$Lineage))) {
  plotcol = subset(colors_cl, Lineage == unique(colors_cl$Lineage)[ii])
  ploteff = subset(efficiency, Cell_line %in% plotcol$Cell_line)
  if(all(ploteff$above_med == 1)) {
    ploteff[nrow(ploteff)+1,] = cbind('dummy', 0, 0, 0)
    stripchart(efficiency ~ above_med, data = ploteff, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = c('transparent', plotcol$outside[1]), bg = plotcol$inside[1], cex = 1.4,
               jitter = 0.3)
  } else {
    stripchart(efficiency ~ above_med, data = ploteff, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = plotcol$outside[1], bg = plotcol$inside[1], cex = 1.4, jitter = 0.3)
  }
}
plot(efficiency$median, efficiency$efficiency, pch = 21, col = colors_cl$outside, bg = colors_cl$inside, bty = 'n',
     lwd = 3, main = substitute(paste(italic(r), ' = ', x, ', ', italic(p), ' = ', y, sep = ''),
                                list(x = format(round(r_eff$estimate, 4), nsmall = 4),
                                     y = format(round(r_eff$p.value, 4), nsmall = 4))), yaxt = 'n', ylab = '', cex = 1.4,
     xlab = 'Median LFC of positive controls')
dev.off()
cairo_pdf('Confounding_Factor_Transduction_Efficiency_BIG.pdf', 20, 10)
par(mfrow=c(1,2))
boxplot(efficiency ~ above_med, data = efficiency, col = 'transparent', axes = FALSE, 
        main = substitute(paste(italic(p), ' = ', x, sep = ''), list(x = format(round(t_eff$p.value, 4), nsmall = 4))),
        ylab = '', cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
        xlab = 'Median LFC of positive controls')
axis(1, at = c(1, 2), labels = c('< grand median', '\u2265 grand median'), lwd = 0, cex.axis = 1.5)
axis(2, cex.axis = 1.5, line = -1.5)
title(ylab = 'Average transduction efficiency', line = 1.3, cex.lab = 1.5)
for (ii in 1:length(unique(colors_cl$Lineage))) {
  plotcol = subset(colors_cl, Lineage == unique(colors_cl$Lineage)[ii])
  ploteff = subset(efficiency, Cell_line %in% plotcol$Cell_line)
  if(all(ploteff$above_med == 1)) {
    ploteff[nrow(ploteff)+1,] = cbind('dummy', 0, 0, 0)
    stripchart(efficiency ~ above_med, data = ploteff, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = c('transparent', plotcol$outside[1]), bg = plotcol$inside[1], cex = 2.2,
               jitter = 0.3)
  } else {
    stripchart(efficiency ~ above_med, data = ploteff, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = plotcol$outside[1], bg = plotcol$inside[1], cex = 2.2, jitter = 0.3)
  }
}
plot(efficiency$median, efficiency$efficiency, pch = 21, col = colors_cl$outside, bg = colors_cl$inside, bty = 'n',
     lwd = 3, main = substitute(paste(italic(r), ' = ', x, ', ', italic(p), ' = ', y, sep = ''),
                                list(x = format(round(r_eff$estimate, 4), nsmall = 4),
                                     y = format(round(r_eff$p.value, 4), nsmall = 4))), yaxt = 'n', ylab = '', cex = 2.2,
     xlab = 'Median LFC of positive controls', cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

##Library coverage at transduction
## Build data frame with coverage at transduction and information on whether the median LFC of positive controls is
## the grand median
coverage = data.frame(Cell_line = inline$Cell.line, coverage = inline$Coverage.at.transduction,
                      median = quantiles[2,][order(names(quantiles[2,]))])
coverage = cbind.data.frame(coverage, above_med = as.numeric(coverage$Cell_line %in% above))
## Statistical tests
t_cov = t.test(subset(coverage$coverage, coverage$above_med == 0),
               subset(coverage$coverage, coverage$above_med == 1))
r_cov = cor.test(coverage$median, coverage$coverage)

## Generate plot with boxplots of the two sets (above and below grand median) on the left and scatterplot to
## illustrate a possible linear correlation between median LFC of positive controls and the coverage at transduction
## of the library in a cell line
cairo_pdf('Confounding_Factor_Coverage.pdf', 20, 10)
par(mfrow=c(1,2))
boxplot(coverage ~ above_med, data = coverage, col = 'transparent', axes = FALSE, 
        main = substitute(paste(italic(p), ' = ', x, sep = ''), list(x = format(round(t_cov$p.value, 4), nsmall = 4))),
        ylab = 'Average library coverage (x) at transduction',
        xlab = 'Median LFC of positive controls')
axis(1, at = c(1, 2), labels = c('< grand median', '\u2265 grand median'), lwd = 0)
axis(2)
for (ii in 1:length(unique(colors_cl$Lineage))) {
  plotcol = subset(colors_cl, Lineage == unique(colors_cl$Lineage)[ii])
  plotcov = subset(coverage, Cell_line %in% plotcol$Cell_line)
  if(all(plotcov$above_med == 1)) {
    plotcov[nrow(plotcov)+1,] = cbind('dummy', 0, 0, 0)
    stripchart(coverage ~ above_med, data = plotcov, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = c('transparent', plotcol$outside[1]), bg = plotcol$inside[1], cex = 1.4,
               jitter = 0.3)
  } else {
    stripchart(coverage ~ above_med, data = plotcov, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = plotcol$outside[1], bg = plotcol$inside[1], cex = 1.4, jitter = 0.3)
  }
}
plot(coverage$median, coverage$coverage, pch = 21, col = colors_cl$outside, bg = colors_cl$inside, bty = 'n',
     lwd = 3, main = substitute(paste(italic(r), ' = ', x, ', ', italic(p), ' = ', y, sep = ''),
                                list(x = format(round(r_cov$estimate, 4), nsmall = 4),
                                     y = format(round(r_cov$p.value, 4), nsmall = 4))), yaxt = 'n', ylab = '', cex = 1.4,
     xlab = 'Median LFC of positive controls')
dev.off()
cairo_pdf('Confounding_Factor_Coverage_BIG.pdf', 20, 10)
par(mfrow=c(1,2))
boxplot(coverage ~ above_med, data = coverage, col = 'transparent', axes = FALSE, 
        main = substitute(paste(italic(p), ' = ', x, sep = ''), list(x = format(round(t_cov$p.value, 4), nsmall = 4))),
        ylab = '', cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5,
        xlab = 'Median LFC of positive controls')
axis(1, at = c(1, 2), labels = c('< grand median', '\u2265 grand median'), lwd = 0, cex.axis = 1.5)
axis(2, cex.axis = 1.5, line = -1.5)
title(ylab = 'Average library coverage (x) at transduction', line = 1.3, cex.lab = 1.5)
for (ii in 1:length(unique(colors_cl$Lineage))) {
  plotcol = subset(colors_cl, Lineage == unique(colors_cl$Lineage)[ii])
  plotcov = subset(coverage, Cell_line %in% plotcol$Cell_line)
  if(all(plotcov$above_med == 1)) {
    plotcov[nrow(plotcov)+1,] = cbind('dummy', 0, 0, 0)
    stripchart(coverage ~ above_med, data = plotcov, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = c('transparent', plotcol$outside[1]), bg = plotcol$inside[1], cex = 2.2,
               jitter = 0.3)
  } else {
    stripchart(coverage ~ above_med, data = plotcov, method = 'jitter', pch = 21, vertical = TRUE, axes = FALSE,
               lwd = 3, add = TRUE, col = plotcol$outside[1], bg = plotcol$inside[1], cex = 2.2, jitter = 0.3)
  }
}
plot(coverage$median, coverage$coverage, pch = 21, col = colors_cl$outside, bg = colors_cl$inside, bty = 'n',
     lwd = 3, main = substitute(paste(italic(r), ' = ', x, ', ', italic(p), ' = ', y, sep = ''),
                                list(x = format(round(r_cov$estimate, 4), nsmall = 4),
                                     y = format(round(r_cov$p.value, 4), nsmall = 4))), yaxt = 'n', ylab = '', cex = 2.2,
     xlab = 'Median LFC of positive controls', cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5)
dev.off()

## ANOVA with confounding factors
## Merge replicates from run (information on sequencing and DNA isolation run)
run = run[,order(colnames(run))]
run_cl = data.frame('NA' = rep(NA, 3))
rownames(run_cl) = rownames(run)
counter_out = 1
ii = 1
while (ii <= ncol(run)) {
  counter_in = 0
  while (gsub('_REP_.', '', colnames(run)[ii]) == gsub('_REP_.', '', colnames(run)[ii + counter_in + 1]) &&
         (ii + counter_in + 1) <= ncol(run)) {
    counter_in = counter_in + 1
  }
  if (table(as.matrix(run)[1,ii:(ii+counter_in)])[1] == length(run[1,ii:(ii+counter_in)])) {
    run_cl[1,counter_out] = as.numeric(run[1,ii])
  } else if (max(table(as.matrix(run)[1,ii:(ii+counter_in)])) == 2) {
    run_cl[1,counter_out] = as.numeric(names(table(as.matrix(run)[1,ii:(ii+counter_in)]))[
      which(table(as.matrix(run)[1,ii:(ii+counter_in)]) == max(table(as.matrix(run)[1,ii:(ii+counter_in)])))])
  } else {
    run_cl[1,counter_out] = as.numeric(names(table(as.matrix(run)[1,ii:(ii+counter_in)]))[
      which(table(as.matrix(run)[1,ii:(ii+counter_in)]) == min(table(as.matrix(run)[1,ii:(ii+counter_in)])))])
  }
  if (table(as.matrix(run)[2,ii:(ii+counter_in)])[1] == length(run[2,ii:(ii+counter_in)])) {
    run_cl[2,counter_out] = run[2,ii]
  } else if (max(table(as.matrix(run)[2,ii:(ii+counter_in)])) == 2) {
    run_cl[2,counter_out] = names(table(as.matrix(run)[2,ii:(ii+counter_in)]))[
      which(table(as.matrix(run)[2,ii:(ii+counter_in)]) == max(table(as.matrix(run)[2,ii:(ii+counter_in)])))]
  } else {
    run_cl[2,counter_out] = names(table(as.matrix(run)[2,ii:(ii+counter_in)]))[
      which(table(as.matrix(run)[2,ii:(ii+counter_in)]) == min(table(as.matrix(run)[2,ii:(ii+counter_in)])))]
  }
  if (table(as.matrix(run)[3,ii:(ii+counter_in)])[1] == length(run[3,ii:(ii+counter_in)])) {
    run_cl[3,counter_out] = as.numeric(run[3,ii])
  } else if (max(table(as.matrix(run)[3,ii:(ii+counter_in)])) == 2) {
    run_cl[3,counter_out] = as.numeric(names(table(as.matrix(run)[3,ii:(ii+counter_in)]))[
      which(table(as.matrix(run)[3,ii:(ii+counter_in)]) == max(table(as.matrix(run)[3,ii:(ii+counter_in)])))])
  } else {
    run_cl[3,counter_out] = as.numeric(names(table(as.matrix(run)[3,ii:(ii+counter_in)]))[
      which(table(as.matrix(run)[3,ii:(ii+counter_in)]) == min(table(as.matrix(run)[3,ii:(ii+counter_in)])))])
  }
  colnames(run_cl)[counter_out] = gsub('_REP_.', '', colnames(run)[ii])
  ii = ii + counter_in + 1
  counter_out = counter_out + 1
}
Seq_run = as.numeric(run_cl[1,])
DNA_run = as.numeric(run_cl[3,])
anova_run = aov(apply(LFC_posc, 2, median) ~ Seq_run * DNA_run)
anova_lin = aov(apply(LFC_posc, 2, median) ~ cl_info$Primary.disease * cl_info$Lineage *
                  cl_info$Lineage.subtype)
anova_rest = aov(apply(LFC_posc, 2, median) ~ mcor * efficiency$efficiency * coverage$coverage)
summary(anova_run) # Interaction barely significant, not anymore after correction for multiple testing
summary(anova_lin) # Nothing significant
summary(anova_rest) # Mean correlation barely significant, not anymore after correction for multiple testing
short = cbind(t(LFC_posc), Seq_run, DNA_run, mcor, efficiency = efficiency$efficiency, coverage = coverage$coverage,
              Primary_disease = cl_info$Primary.disease, Lineage = cl_info$Lineage,
              Lineage_subtype = cl_info$Lineage.subtype)
long = melt(as.data.table(short), id.vars = 11:ncol(short), variable.name = 'pos_cntrl')
# summary(aov(value ~ Seq_run * DNA_run * mcor * efficiency * coverage * Primary_disease * Lineage *
# Lineage_subtype , data = long)) # too many

write.xlsx(run_cl, 'Information on runs on cell line level.xlsx', rowNames = TRUE)

## Sanity Checks
SanityCheck1 = all(c(below, above) %in% colnames(quantiles)) && all(colnames(quantiles) %in% c(below, above))
SanityCheck2 = all(cl_info$Lineage %in% colors_lin[,1])
SanityCheck3 = all(rownames(mcor_ab) == colors_cl$Cell_line)#cl_info$Cell.line)
SanityCheck4 = all(inline$Cell.line == colors_cl$Cell_line)


# SanityCheck1 # Are all names from the split set of cell lines above and below the grand median preserved?
# SanityCheck2 # Are the names for lineages formatted in the same manner?
# SanityCheck3 # Are the data frames for the color coding of cell lines and for the plot of correlation coefficients ordered in the same manner?
all(SanityCheck1,SanityCheck2,SanityCheck3)



