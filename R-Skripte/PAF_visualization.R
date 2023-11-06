#####SKRIPT TO VISUALIZE PAF#####
install.packages("remotes")
remotes::install_github("trevorld/r-optparse")
library(pafr)
library(optparse)
library(ggplot2)
paf <- read_paf("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/SE02_scaffolds_witch_correction_LR_updated_ref.paf")
paf_corrected_Scaff_against_modif_ref <- read_paf("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/SE02_scaffolds_witch_correction_LR_updated_ref.paf")

dotplot(paf_corrected_Scaff_against_modif_ref)

dotplot(paf_corrected_Scaff_against_modif_ref, label_seqs=TRUE, order_by="qstart", alignment_colour="blue")
paf$tname
ggplot(data=paf, aes(x=V3, xend=V4, y=V8, yend=V9, color=factor(V1)),show.legend = FALSE ) + geom_segment() + labs(x="Reference coordinate", y="Query coordinate")
paf$qname
plot_coverage(paf, fill='qname', direct_label = F)
to_keep <- list(paf$qname,
                c("chr1", "chr2", "chr3", "chr4","chr5", "chr6", "chr7", "chr8")
)
dotplot(paf, label_seqs=TRUE,order_by="qstart")
to_keep <- list(
  c("chr1_RagTag", "chr2_RagTag", "chr3_RagTag", "chr4_RagTag","chr5_RagTag", "chr6_RagTag", "chr7_RagTag", "chr8_RagTag"),
  c("chr1", "chr2", "chr3", "chr4","chr5", "chr6", "chr7", "chr8")
)
dotplot(paf, label_seqs=TRUE, order_by="provided", ordering=to_keep, line_size = 0.5,
        xlab = "Pseudochromosomes SE02_002 corrected",
        ylab = "Reference v5.1")
plot_synteny(paf, q_chrom="chr1_RagTag", t_chrom="chr1", centre=TRUE)
plot_synteny(paf, q_chrom="chr4_RagTag", t_chrom="chr1", centre=TRUE)
plot_synteny(paf, q_chrom="chr4_RagTag", t_chrom="chr4", centre=TRUE)
plot_synteny(paf, q_chrom="chr5_RagTag", t_chrom="chr5", centre=TRUE, xlab ="SE02_corrected")
plot_synteny()

plot_coverage(paf, fill='qname') + scale_fill_brewer(palette="Set1")

?dotplot

#Asynt
#plotsr
################################OWN FUNCTION#########################################

library(pafr)
library(optparse)
library(ggplot2)


option_list <- list(
  make_option(c("-o","--output"), type="character",
              help="output filename prefix [input.paf]", 
              dest="output_filename"),
  make_option(c("-p","--plot-size"), type="numeric", default=15,
              help="plot size X inches [%default]",
              dest="plot_size"),
  make_option(c("-f","--flip"), action="store_true", default=FALSE,
              help="flip query if most alignments are in reverse complement [%default]",
              dest="flip"),
  make_option(c("-b","--break-point"), action="store_true", default=FALSE,
              help="show break points [%default]",
              dest="break_point"),
  make_option(c("-s", "--sort-by-refid"), action="store_true", default=TRUE,
              help="sort reference IDs in alphabetical order, default by length [%default]",
              dest="sortbyID"),
  make_option(c("-q", "--min-query-length"), type="numeric", default=40000,
              help="filter queries with total alignments less than cutoff X bp [%default]",
              dest="min_query_aln"),
  make_option(c("-m", "--min-alignment-length"), type="numeric", default=10000,
              help="filter alignments less than cutoff X bp [%default]",
              dest="min_align"),
  make_option(c("-r", "--min-ref-len"), type="numeric", default=1000000,
              help="filter references with length less than cutoff X bp [%default]",
              dest="min_ref_len"),
  make_option(c("-i", "--reference-ids"), type="character", default=NULL,
              help="comma-separated list of reference IDs to keep and order [%default]",
              dest="refIDs")
)

options(error=traceback)
parser <- OptionParser(usage = "%prog [options] input.paf\n\nFor more information, see https://github.com/moold/paf2dotplot", option_list = option_list)
opts = parse_args(parser, positional_arguments = c(0, 1))
opt = opts$options
input_file = opts$args
#if(length(input_file) <= 0){
#  cat(sprintf("Error: missing input file: input.paf!\n\n"))
#  print_help(parser)
#  quit()
#}else if (file.access(input_file, mode=4) == -1){
#  cat(sprintf("Error: input file: %s does not exist or cannot be read!\n\n", input_file))
#  print_help(parser)
#  quit()
#}
if(is.null(opt$output_filename)){
  opt$output_filename = input_file
}


# save
#ggsave(filename = paste0(opt$output_filename, ".pdf"), width = opt$plot_size, height = opt$plot_size * 0.8, units = "in", dpi = 300, limitsize = F)
#options(warn=0) # turn on warnings

Dotplot_function <- function(paf_file, name_of_x){
  library(pafr)
  library(optparse)
  library(ggplot2)
  
  # read in alignments
  alignments = paf_file
  # avoid inter overflow
  alignments[, c(2:4, 7:12)] = apply(alignments[, c(2:4, 7:12)], 2, as.double)
  
  # set column names
  # PAF IS ZERO-BASED - CHECK HOW CODE WORKS
  colnames(alignments)[1:12] = c("queryID","queryLen","queryStart","queryEnd","strand","refID","refLen","refStart","refEnd","numResidueMatches","lenAln","mapQ")
  
  # caculate similarity
  alignments$percentID = alignments$numResidueMatches / alignments$lenAln
  queryStartTemp = alignments$queryStart
  # Flip starts, ends for negative strand alignments
  alignments$queryStart[which(alignments$strand == "-")] = alignments$queryEnd[which(alignments$strand == "-")]
  alignments$queryEnd[which(alignments$strand == "-")] = queryStartTemp[which(alignments$strand == "-")]
  rm(queryStartTemp)
  
  cat(paste0("\nNumber of alignments: ", nrow(alignments), "\n"))
  cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)), "\n"))
  
  # sort by ref chromosome sizes, keep top X chromosomes OR keep specified IDs
  if(is.null(opt$refIDs)){
    if (opt$sortbyID){
      refIDsToKeepOrdered = unique(sort(alignments$refID))
    }else{
      chromMax = tapply(alignments$refLen, alignments$refID, max)
      refIDsToKeepOrdered = names(sort(chromMax, decreasing = T))
    }
  }else{
    refIDsToKeepOrdered = unlist(strsplit(opt$refIDs, ","))
    alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
  }
  
  # filter queries by alignment length, for now include overlapping intervals
  queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
  alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]
  # filter alignment by length
  alignments = alignments[which(alignments$lenAln > opt$min_align),]
  # filter alignment by ref length
  alignments = alignments[which(alignments$refLen > opt$min_ref_len),]
  # re-filter queries by alignment length, for now include overlapping intervals
  queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
  alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]
  
  cat(paste0("\nAfter filtering... Number of alignments: ", nrow(alignments),"\n"))
  cat(paste0("After filtering... Number of query sequences: ", length(unique(alignments$queryID)),"\n\n"))
  
  # filter refIDsToKeepOrdered if some refID are filtered
  refIDsToKeepOrdered = refIDsToKeepOrdered[which(refIDsToKeepOrdered %in% alignments$refID)]
  
  # sort df on ref
  alignments$refID = factor(alignments$refID, levels = refIDsToKeepOrdered) # set order of refID
  alignments = alignments[with(alignments,order(refID,refStart)),]
  chromMax = tapply(alignments$refLen, alignments$refID, max)
  # make new ref alignments for dot plot
  
  alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
  alignments$refEnd2 = alignments$refEnd + sapply(as.character(alignments$refID), function(x) ifelse(x == names(chromMax)[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
  
  ## queryID sorting step 1/2
  # sort levels of factor 'queryID' based on longest alignment
  alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID)))
  queryMaxAlnIndex = tapply(alignments$lenAln, alignments$queryID, which.max, simplify = F)
  alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID))[order(mapply(
    function(x, i)
      alignments$refStart2[which(i == alignments$queryID)][x],
    queryMaxAlnIndex,
    names(queryMaxAlnIndex)
  ))])
  
  ## queryID sorting step 2/2
  ## sort levels of factor 'queryID' based on longest aggregrate alignmentst to refID's
  # per query ID, get aggregrate alignment length to each refID 
  queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x) tapply(alignments$lenAln[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum) )
  if(length(levels(alignments$refID)) > 1){
    queryID_Ref = apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
  } else {
    queryID_Ref = sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])
  }
  # set order for queryID
  alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])
  queryMax = tapply(alignments$queryLen, alignments$queryID, max)
  
  if(opt$flip){
    #  flip query starts stops to forward if most align are in reverse complement
    queryRevComp = tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)) < 0
    queryRevComp = names(queryRevComp)[which(queryRevComp)]
    alignments$queryStart[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
    alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1
  }
  ## make new query alignments for dot plot
  alignments$queryStart2 = alignments$queryStart + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
  alignments$queryEnd2 = alignments$queryEnd + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
  
  # plot break points
  if (opt$break_point) {
    break_size = 1;
    alignments$break_col = rep(0, length(alignments$percentID));
  }else{
    break_size = 0.009;
    alignments$break_col = alignments$percentID;
  }
  
  options(warn = -1) # turn off warnings
  gp = ggplot(alignments) + 
    theme_bw() + 
    theme(
      text = element_text(size = 12),
      panel.grid.minor.y = element_blank(),
      panel.grid.minor.x = element_blank(),
      axis.text.y = element_text(angle = 15),
      axis.text.x = element_text(hjust = 1, angle = 45)
    ) +
    scale_color_distiller(palette = "Spectral") 
  
  # x
  if (length(unique(alignments$refID)) == 1){
    reflen = unique(alignments$refLen)
    xbreaks = seq(0, reflen, reflen/10)
    if (reflen/10 > 1e9){
      xlables = paste(round(xbreaks/1e9), "GB")
    }else if (reflen/10 > 1e6) {
      xlables = paste(round(xbreaks/1e6), "MB")
    }else if(reflen/10 > 1e3){
      xlables = paste(round(xbreaks/1e3), "KB")
    }else{
      xlables = paste(round(xbreaks), "bp")
    }
    
    gp = gp + scale_x_continuous(expand = c(0, 0), limits = c(0, reflen + 0.1), breaks = xbreaks, labels = xlables) +
      xlab(unique(alignments$refID))
  }else{
    gp = gp + 
      theme(panel.grid.major.x=element_blank()) +
      geom_vline(xintercept = cumsum(as.numeric(chromMax)), col="#ebebeb") + 
      scale_x_continuous(expand = c(0, 0), limits = c(0, sum(as.numeric(chromMax)) + 0.1), 
                         breaks = cumsum(as.numeric(chromMax)) - chromMax/2, labels = substr(levels(alignments$refID), start = 1, stop = 20)) + 
      xlab("original Reference v5.1")
  }
  
  # y
  if (length(unique(alignments$queryID)) == 1){
    queryLen = unique(alignments$queryLen)
    ybreaks = seq(0, queryLen, queryLen/10)
    if (queryLen/10 > 1e9){
      ylables = paste(round(ybreaks/1e9), "GB")
    }else if (queryLen/10 > 1e6) {
      ylables = paste(round(ybreaks/1e6), "MB")
    }else if(queryLen/10 > 1e3){
      ylables = paste(round(ybreaks/1e3), "KB")
    }else{
      ylables = paste(round(ybreaks), "bp")
    }
    
    gp = gp + scale_y_continuous(expand = c(0, 0), limits = c(0, queryLen + 0.1), breaks = ybreaks, labels = ylables) +
      ylab(unique(alignments$queryID))
  }else{
    gp = gp + 
      theme(panel.grid.major.y=element_blank()) +
      geom_hline(yintercept = cumsum(as.numeric(queryMax)), col="#ebebeb") + 
      scale_y_continuous(expand = c(0, 0), limits = c(0, sum(as.numeric(queryMax)) + 0.1), 
                         breaks = cumsum(as.numeric(queryMax)) - queryMax/2, labels = substr(levels(alignments$queryID), start = 1, stop = 20)) + 
      ylab(name_of_x)
  }
  
  # co-line
  gp = gp + 
    geom_point(mapping = aes(x = refStart2, y = queryStart2, color = break_col), size = break_size) +
    geom_point(mapping = aes(x = refEnd2, y = queryEnd2, color = break_col), size = break_size) +
    geom_segment(aes(x = refStart2, xend = refEnd2, y = queryStart2, yend = queryEnd2, color = percentID))
  gp
  
}

#SE02_002
#with correction
read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/SE02_scaffolds_witch_correction_LR_updated_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T)
paf_file <- read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/SE02_scaffolds_witch_correction_LR_updated_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T)
name_of_x <- "SE02_with_LR_correction"
Dotplot_function(paf_file, name_of_x)

#without correction
Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/SE02_scaffolds_wo_correction_updated_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "SE02_without_correction")

#orig ref
Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/SE02_chr_against_v5.1.paf", stringsAsFactors = F, row.names=NULL, fill = T), "SE02_with_correction")


#NO02_001
#with correction 
Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/nor02_001_hifi/NO02_scaffolds_witch_correction_LR_updated_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "NO02_001_with_correction")

############ FINAL ASSEMBLIES ###########################
library(ggpubr)
?ggarrange

AT <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/AT03_004_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "AT03_004 pseudochromosomes")


CH <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/CH05_002_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "CH05_002 pseudochromosomes")


CZ <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/CZ01_005_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "CZ01_005 pseudochromosomes")


DE <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/DE01_001_scaffolds_against_ref_ragtag_minimap.paf", stringsAsFactors = F, row.names=NULL, fill = T), "DE01_001 pseudochromosomes")

ggarrange(ncol = 2, nrow = 2, AT, CH, CZ, DE)


ES03 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/ES03_014_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "ES03_014 pseudochromosomes")


ES04 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/ES04_014_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "ES04_014 pseudochromosomes")


ES06 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/ES06_006_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "ES06_006 pseudochromosomes")


ES18 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/ES18_040_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "ES18_040 pseudochromosomes")
ggarrange(ncol = 2, nrow = 2, ES03, ES04, ES06, ES18)


ES20 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/ES20_028_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "ES20_028 pseudochromosomes")


ES24 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/ES24_001_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "ES24_001 pseudochromosomes")


ES27 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/ES27_024_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "ES27_024 pseudochromosomes")


IS01 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/IS01_007_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "IS01_007 pseudochromosomes")



NO01 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/NO01_005_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "NO01_005 pseudochromosomes")


NO02 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/NO02_001_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "NO02_001 pseudochromosomes")


NO04 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/NO04_015_scaffolds_against_ref_ragtag_minimap.paf", stringsAsFactors = F, row.names=NULL, fill = T), "NO04_015 pseudochromosomes")


PL <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/PL01_004_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "PL01_004 pseudochromosomes")



PT <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/PT01_005_scaffolds_against_ref_ragtag_minimap.paf", stringsAsFactors = F, row.names=NULL, fill = T), "PT01_005 pseudochromosomes")


SE02 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/SE02_002_scaffolds_against_ref_ragtag_minimap.paf", stringsAsFactors = F, row.names=NULL, fill = T), "SE02_002 pseudochromosomes")


SE06 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/SE06_020_scaffolds_against_ref_ragtag_minimap.paf", stringsAsFactors = F, row.names=NULL, fill = T), "SE06_020 pseudochromosomes")


SE07 <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap/SE07_013_scaffolds_against_ref_ragtag_minimap.paf", stringsAsFactors = F, row.names=NULL, fill = T), "SE07_013 pseudochromosomes")
ggarrange(ncol = 2, nrow = 2, PT, SE02, SE06,SE07)


SI <- Dotplot_function(read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/paf_files_minimap_purged/SI01_005_scaffolds_against_ref.paf", stringsAsFactors = F, row.names=NULL, fill = T), "SL01_005 pseudochromosomes")

ggarrange(ncol = 2, nrow = 5, AT, CH, CZ, DE, ES03, ES04, ES06, ES18, ES20, ES27)
ggarrange(ncol = 2, nrow = 2, ES20, ES27, IS01, SI)
ggarrange(ncol = 2, nrow = 2, NO01, NO02, NO04,PL)

ggarrange(ncol = 2, nrow = 5 ,  IS01, NO01, NO02, NO04,PL,  PT, SE02, SE06,SE07,SI)
ggarrange(ncol = 2, nrow = 2, ES20, ES27, IS01, SI)

################## COVERAGE ##############
devtools::install_github("hadley/reshape")
library(reshape) # loads the library to rename the column names
install.packages("reshape2")
library(reshape2)

SE02_chr1_coverage <- read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/pilon/SE02_002_SR_chr1.coverage",
              header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
names(SE02_chr1_coverage) <- c("Chr", "Nucleotide", "depth")


# Read the data from the coverage file
coverage_data <- read.table("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/SE02_002_LR_chr1.coverage",
                            header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
#Read in coverage table
coverage_data_2 <- read.table("SE02_002_LR_chr2.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coverage_data_3 <- read.table("SE02_002_LR_chr3.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coverage_data_4 <- read.table("SE02_002_LR_chr4.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coverage_data_5 <- read.table("SE02_002_LR_chr5.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coverage_data_6 <- read.table("SE02_002_LR_chr6.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coverage_data_7 <- read.table("SE02_002_LR_chr7.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
coverage_data_8 <- read.table("SE02_002_LR_chr8.coverage",header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
#change headers of table
names(coverage_data_2) <- c("Chr", "Nucleotide", "depth")
names(coverage_data_3) <- c("Chr", "Nucleotide", "depth")
names(coverage_data_4) <- c("Chr", "Nucleotide", "depth")
names(coverage_data_5) <- c("Chr", "Nucleotide", "depth")
names(coverage_data_6) <- c("Chr", "Nucleotide", "depth")
names(coverage_data_7) <- c("Chr", "Nucleotide", "depth")
names(coverage_data_8) <- c("Chr", "Nucleotide", "depth")


#Compute windows of 1000b
# Assuming you have a data frame named coverage_data
# with columns "Chr," "Nucleotide," and "depth"

# Load the necessary library
library(dplyr)

# Define the window size
window_size <- 1000

# Create a function to calculate the mean coverage for a window
calculate_mean_coverage <- function(data) {
  return(data %>%
           summarise(Mean_Coverage = mean(depth)))
}

# Group the data by Chromosome and calculate mean coverage within each window
result <- SE02_chr1_coverage %>%
  group_by(Chr, Window = ((Nucleotide - 1) %/% window_size) * window_size + 1) %>%
  do(calculate_mean_coverage(.)) %>%
  ungroup()

# Print or save the result
View(result)
View(result[40000:50000,])
sorted_cov <- result[order(result$Mean_Coverage), ]
View(result[90748:90760,])
# You can also save the result to a file if needed
# write.table(result, "mean_coverage_per_window.txt", sep = "\t", row.names = FALSE, col.names = TRUE)


# Define the window size
View(result[27000:290000,])
no_zoom <- ggplot(result, aes(x=Window, y=Mean_Coverage)) +
  geom_point(size=0.5)+
  ggtitle(label="Coverage of reads against Scaffolds", subtitle = "SE02-002: Chr1" )+
  geom_vline(xintercept = c(11280001), color = "lightgrey")+
  geom_vline(xintercept = c(13332001), color = "lightgrey") +
  annotate('rect', xmin=11280001, xmax=13332001, ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgrey')+
  annotate('rect', xmin=4753001, xmax=4772001, ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen')+
  annotate('rect', xmin=7954001, xmax=7970001, ymin=-Inf, ymax=Inf, alpha=.2, fill='lightgreen')+
  annotate('rect', xmin=15988001, xmax=16003001, ymin=-Inf, ymax=Inf, alpha=.2, fill='lightblue')+
  annotate('rect', xmin=29000001, xmax=29018001, ymin=-Inf, ymax=Inf, alpha=.2, fill='red')+
  scale_x_continuous(labels = scales::comma)+
  xlab(label = "Chromosome position [bp]")+
  scale_fill_manual(values = c('lightgrey' = 'lightgrey', 'lightgreen' = 'lightgreen', 'blue' = 'blue', 'red' = 'red'),
                    labels = c('Centromere', 'Plastid', 'rDNA', 'Duplication'))

# Create a data frame for the annotations
# Create a data frame for the annotations
annotations <- data.frame(
  xmin = c(11280001, 4753001, 7954001, 15988001, 29000001),
  xmax = c(13332001, 4772001, 7970001, 16003001, 29018001),
  ymin = -Inf,
  ymax = Inf,
  alpha = c(0.2, 0.2, 0.2, 0.2, 0.2),
  fill = c('lightgrey', 'lightgreen', 'lightgreen', 'lightblue', 'red')
)
coverage_vector <- rep("Coverage", times = 29832)
result$Annotation <- coverage_vector
length()
unique(result$Annotation)
result[11280,4] <- "Centromere"
result[4603,4] <- "Plastid"
result[11280,4] <- "Centromere"
result[15878,4] <- "rDNA"
result[28900,4] <- "Duplication"
# Create the ggplot
zoom <- ggplot(result, aes(x = Window, y = Mean_Coverage, color=Annotation))+
  scale_color_manual(values=c("Centromere"="lightgrey", "Plastid"="lightgreen","rDNA"="blue","Duplication"="red", "Coverage"="black"),
                     guide = guide_legend(override.aes = list(size = 3) ))+ 
  geom_point(size = 0.5)  +
  annotate('rect', xmin = c(11280001, 4603001, 7844001, 15878001, 28900000), xmax = c(13332001, 4950001, 8045001, 16133001, 29068001), ymin = -Inf, ymax = Inf, fill = c('lightgrey', 'lightgreen', 'lightgreen', 'blue', 'red'), alpha = 0.2) +
  scale_x_continuous(labels = scales::comma) +
  xlab(label = "Chromosome position [bp]")+
  theme_minimal() +
  ylab(label = "Coverage")+
  ylim(0,200)


no_zoom <- ggplot(result, aes(x = Window, y = Mean_Coverage,color=Annotation))+
  scale_color_manual(values=c("Centromere"="lightgrey", "Plastid"="lightgreen","rDNA"="blue","Duplication"="red", "Coverage"="black"),
                     guide = guide_legend(override.aes = list(size = 3) ))+ 
  geom_point(size = 0.5) +
  ggtitle(label = "SE02-002: Chromosome 1")  +
  annotate('rect', xmin = c(11280001, 4603001, 7844001, 15878001, 28900000), xmax = c(13332001, 4950001, 8045001, 16133001, 29068001), ymin = -Inf, ymax = Inf, fill = c('lightgrey', 'lightgreen', 'lightgreen', 'blue', 'red'), alpha = 0.2) +
  scale_x_continuous(labels = scales::comma) +
  xlab(label = "Chromosome position [bp]")+
  theme_minimal() +
  ylab(label = "Coverage")
library(ggpubr)
plot_chr1 <- ggarrange(no_zoom, zoom, ncol=1, nrow=2)
mean(result$Mean_Coverage)

version






library(dplyr)
library(ggplot2)
library(ggpubr)

# List of data frames
coverage_data_list <- list(coverage_data_1, coverage_data_2, coverage_data_3, coverage_data_4, 
                           coverage_data_5, coverage_data_6, coverage_data_7, coverage_data_8)

# List of chromosome names
chromosomes <- c("chr1","chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8")

# Directory to save the CSV files
output_directory <-  "/netscratch/dep_coupland/grp_fulgione/male/assemblies/SE02_002/RagTag/ragtag_without_correction_updated_ref"
# Change to your desired directory path

# Create the output directory if it doesn't exist
#dir.create(output_directory, showWarnings = FALSE)

# Loop through each data frame
for (i in 1:length(coverage_data_list)) {
  # Get the current coverage data frame
  coverage_data <- coverage_data_list[[i]]
  
  # Define the window size
  window_size <- 1000
  
  # Group the data by Chromosome and calculate mean coverage within each window
  result <- coverage_data %>%
    group_by(Chr, Window = ((Nucleotide - 1) %/% window_size) * window_size + 1) %>%
    summarise(Mean_Coverage = mean(depth)) %>%
    ungroup()
  
  # Create a unique output file name for each chromosome
  output_file_name <- file.path(output_directory, paste("chr_", chromosomes[i], "_coverage_mean.csv", sep=""))
  
  # Save the mean coverage data frame to a CSV file
  write.csv(result, output_file_name, row.names = FALSE)
  
  # Print a message indicating the completion of processing for the current chromosome
  cat(paste("Processed", chromosomes[i], "coverage data and saved to", output_file_name, "\n"))
}


##########################################################################

coverage_data_1 <- result
coverage_data_1 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr1_coverage_mean.csv")
coverage_data_2 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr2_coverage_mean.csv")
coverage_data_3 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr3_coverage_mean.csv")
coverage_data_4 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr4_coverage_mean.csv")
coverage_data_5 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr5_coverage_mean.csv")
coverage_data_6 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr6_coverage_mean.csv")
coverage_data_7 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr7_coverage_mean.csv")
coverage_data_8 <- read.csv("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/chr_chr8_coverage_mean.csv")


# Create an empty list to store the plots
all_plots <- list()

# Loop through each data frame
for (i in 1:length(coverage_data_list)) {
  # Get the current coverage data frame
  coverage_data <- coverage_data_list[[i]] 
  
  # Create a ggplot plot
  no_zoom <- ggplot(coverage_data, aes(x=Window, y=Mean_Coverage)) +
    geom_point(size=0.5)+
    ggtitle(label=paste("Coverage of reads against Scaffolds"))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(labels = scales::comma)+
    xlab(label = "position on chromosome [bp]")
  
  zoom <- ggplot(coverage_data, aes(x=Window, y=Mean_Coverage)) +
    geom_point(size=0.5)+
    ylim(0,200)+
    ggtitle(label=paste("Coverage of reads against Scaffolds (Zoom)"))+
    theme(plot.title = element_text(hjust = 0.5))+
    scale_x_continuous(labels = scales::comma)+
    xlab(label = "position on chromosome [bp]")#+ 
  #geom_text(x = (13332001 + 1000) , y = 200, label = "pericentromere", color = "red", hjust = 0.5, vjust = -0.5)
  
  # Arrange the plots using ggpubr
  plot_chr <- ggarrange(no_zoom, zoom)
  
  # Append the current plot to the list of plots for all chromosomes
  all_plots[[i]] <- plot_chr
}

# Now, all_plots is a list containing one plot for each chromosome.
# You can access individual plots like this:
plot_for_chr1 <- all_plots[[1]] # The plot for chr2
plot_for_chr2 <- all_plots[[2]] # The plot for chr3
plot_for_chr3 <- all_plots[[3]]
plot_for_chr4 <- all_plots[[4]]
plot_for_chr5 <- all_plots[[5]]
plot_for_chr6 <- all_plots[[6]]
plot_for_chr7 <- all_plots[[7]]
plot_for_chr8 <- all_plots[[8]]

chr1_4 <- ggarrange(plot_for_chr1, plot_for_chr2,plot_for_chr3,plot_for_chr4,
          nrow = 4, ncol = 1,
          labels = c("chr1", "chr2", "chr3", "chr4"))
chr5_8 <- ggarrange(plot_for_chr5, plot_for_chr6,plot_for_chr7,plot_for_chr8,
          nrow = 4, ncol = 1,
          labels = c("chr5", "chr6", "chr7", "chr8"))
# Save the plot to a PDF file
pdf_file_name <- paste("chr_", chr, "_coverage_plot.pdf", sep="")
ggsave(pdf_file_name, plot_chr, width = 10, height = 5)

setwd("~/Studium/Köln_Biological Science/Master Thesis/assemblies/SE02_002/")
ggsave("Chromosome1-4_depth.pdf", chr1_4, width = 20, height = 20)
ggsave("Chromosome5-8_depth.pdf", chr5_8, width = 20, height = 20)
