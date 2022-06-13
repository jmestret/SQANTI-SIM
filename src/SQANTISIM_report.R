#######################################
#                                     #
#    SQANTISIM report generation     #
#                                     #
#######################################

# Author: Jorge Mestre Tomas (jormart2@alumni.uv.es)

#######################################
#                                     #
#      PACKAGES AND LIBRARIES         #
#                                     #
#######################################

suppressMessages(library(dplyr))
suppressMessages(library(DT))
suppressMessages(library(fmsb))
suppressMessages(library(ggplot2))
suppressMessages(library(gridExtra))
suppressMessages(library(knitr))
suppressMessages(library(rmarkdown))
suppressMessages(library(tidyr))


#######################################
#                                     #
#      FUNCTIONS  AND CLASSES         #
#                                     #
#######################################

get_performance_metrics <- function(data.query, data.index, MAX_TSS_TTS_DIFF, min_supp=0){
  # Simulated ref transcripts
  if (min_supp > 0){
    idx <- data.index[which(data.index$sim_counts >= min_supp), ]
  } else {
    idx <- data.index
  }
  data.novel <- idx[which(idx$sim_type == 'novel'),]
  data.known <- idx[which(idx$sim_type == 'known'),]
  sim.sc <- unique(data.novel$structural_category)

  # Matches between simulated and reconstructed transcripts:
  # First all splice-junctions must be identical
  # Second, difference between the annotated and reconstructed TSS and TTS must be smaller than MAX_TSS_TTS_DIFF
  matches <- inner_join(data.query, data.index, by=c('junctions', 'chrom')) %>%
    mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
    arrange(difftot) %>%
    distinct(isoform, .keep_all = T)
  matches <- matches[which(matches$sim_type != "absent"),]
  perfect.matches <- matches[which(matches$diffTSS < MAX_TSS_TTS_DIFF & matches$diffTTS < MAX_TSS_TTS_DIFF),]
  cond <- (perfect.matches$exons > 1) | (perfect.matches$strand.x == '+' & perfect.matches$TSS_genomic_coord.x <= perfect.matches$TTS_genomic_coord.y & perfect.matches$TSS_genomic_coord.y <= perfect.matches$TTS_genomic_coord.x) | (perfect.matches$strand.x == '-' & perfect.matches$TTS_genomic_coord.x <= perfect.matches$TSS_genomic_coord.y & perfect.matches$TTS_genomic_coord.y <= perfect.matches$TSS_genomic_coord.x)
  perfect.matches <- perfect.matches[cond,]
  matches <- matches[which(matches$junctions != "" | matches$isoform %in% perfect.matches$isoform),]  # Delete PTP of mono-exons which always match because there is no splice junctions
  
  known.matches <- matches[which(matches$sim_type == "known"),]
  known.perfect.matches <- perfect.matches[which(perfect.matches$sim_type == "known"),]
  novel.matches <- matches[which(matches$sim_type == "novel"),]
  novel.perfect.matches <- perfect.matches[which(perfect.matches$sim_type == "novel"),]
  
  #  Compute metrics
  sqantisim.stats <- data.frame(init=c())
  for (sc in xaxislabelsF1){
    
    if (sc == 'FSM') {
      n.sim <- nrow(data.known)
      TP <- nrow(known.perfect.matches[which(known.perfect.matches$structural_category.x == sc),])
      PTP <- nrow(known.matches[which(known.matches$structural_category.x == sc),]) - TP
      FP <- nrow(data.query[which(data.query$structural_category == sc),]) - TP - PTP
    } else {
      n.sim <- nrow(data.novel[which(data.novel$structural_category == sc),])
      TP <- nrow(novel.perfect.matches[which(novel.perfect.matches$structural_category.y == sc),])
      PTP <- nrow(novel.matches[which(novel.matches$structural_category.y == sc),]) - TP
      FP <- nrow(data.query[which(data.query$structural_category == sc),]) - TP - (nrow(novel.matches[which(novel.matches$structural_category.x == sc),]) - TP)
    }
    FN <- n.sim - TP

    if (sum(n.sim, TP, PTP, FP, FN) > 0){
      sqantisim.stats['Total', sc] <- n.sim
      sqantisim.stats['TP', sc] <- TP
      sqantisim.stats['PTP', sc] <- PTP
      sqantisim.stats['FP', sc] <- FP
      sqantisim.stats['FN', sc] <- FN
      sqantisim.stats['Sensitivity', sc] <- TP/ (TP + FN)
      sqantisim.stats['Precision', sc] <- TP/ (TP + PTP + FP)
      sqantisim.stats['F-score', sc] <- 2*((sqantisim.stats['Sensitivity', sc]*sqantisim.stats['Precision', sc])/(sqantisim.stats['Sensitivity', sc]+sqantisim.stats['Precision', sc]))
      sqantisim.stats['False_Discovery_Rate', sc] <- (FP + PTP) / (FP + PTP +  TP)
      sqantisim.stats['Positive_Detection_Rate', sc] <- (TP + PTP) / n.sim
      sqantisim.stats['False_Detection_Rate', sc] <- (FP) / (FP + PTP + TP)
    }
  }
  
  known.TP <- 0
  known.PTP <- 0
  known.FN <- 0
  novel.TP <- 0
  novel.PTP <- 0
  novel.FN <- 0 
  total.FP <- 0
  for (sc in colnames(sqantisim.stats)){
    if (sc == "FSM"){
      known.TP <- known.TP + sqantisim.stats['TP', sc]
      known.PTP <- known.PTP + sqantisim.stats['PTP', sc]
      total.FP <- total.FP + sqantisim.stats['FP', sc]
    }else{
      novel.TP <- novel.TP + sqantisim.stats['TP', sc]
      novel.PTP <- novel.PTP + sqantisim.stats['PTP', sc]
      total.FP <- total.FP + sqantisim.stats['FP', sc]
    }
  }
  
  sqantisim.stats['Total', 'Known'] <-  nrow(data.known)
  sqantisim.stats['TP', 'Known'] <- known.TP
  sqantisim.stats['PTP', 'Known'] <- known.PTP
  sqantisim.stats['FP', 'Known'] <- total.FP
  sqantisim.stats['FN', 'Known'] <- nrow(data.known) - known.TP
  sqantisim.stats['Precision', 'Known'] <- known.TP / (known.TP + known.PTP + total.FP)
  sqantisim.stats['Sensitivity', 'Known'] <- known.TP / (known.TP + sqantisim.stats['FN', 'Known'])
  sqantisim.stats['F-score', 'Known'] <- 2*((sqantisim.stats['Sensitivity', 'Known']*sqantisim.stats['Precision', 'Known'])/(sqantisim.stats['Sensitivity', 'Known']+sqantisim.stats['Precision', 'Known']))
  sqantisim.stats['False_Discovery_Rate', 'Known'] <- (total.FP + known.PTP) / (total.FP + known.PTP +  known.TP)
  sqantisim.stats['Positive_Detection_Rate', 'Known'] <- (known.TP + known.PTP) / nrow(data.known)
  sqantisim.stats['False_Detection_Rate', 'Known'] <- (total.FP) / (total.FP + known.PTP +  known.TP)
  
  sqantisim.stats['Total', 'Novel'] <-  nrow(data.novel)
  sqantisim.stats['TP', 'Novel'] <- novel.TP
  sqantisim.stats['PTP', 'Novel'] <- novel.PTP
  sqantisim.stats['FP', 'Novel'] <- total.FP
  sqantisim.stats['FN', 'Novel'] <- nrow(data.novel) - novel.TP
  sqantisim.stats['Precision', 'Novel'] <- novel.TP / (novel.TP + novel.PTP + total.FP)
  sqantisim.stats['Sensitivity', 'Novel'] <- novel.TP / (novel.TP + sqantisim.stats['FN', 'Novel'])
  sqantisim.stats['F-score', 'Novel'] <- 2*((sqantisim.stats['Sensitivity', 'Novel']*sqantisim.stats['Precision', 'Novel'])/(sqantisim.stats['Sensitivity', 'Novel']+sqantisim.stats['Precision', 'Novel']))
  sqantisim.stats['False_Discovery_Rate', 'Novel'] <- (total.FP + novel.PTP) / (total.FP + novel.PTP +  novel.TP)
  sqantisim.stats['Positive_Detection_Rate', 'Novel'] <- (novel.TP + novel.PTP) / nrow(data.novel)
  sqantisim.stats['False_Detection_Rate', 'Novel'] <- (total.FP) / (total.FP + novel.PTP +  novel.TP)
  
  col.order <- c("Known", "Novel", "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
  row.order <- c('Total', 'TP', 'PTP', 'FP', 'FN', 'Sensitivity', 'Precision', 'F-score', 'False_Discovery_Rate', 'Positive_Detection_Rate', 'False_Detection_Rate')
  sqantisim.stats <- sqantisim.stats[intersect(row.order, rownames(sqantisim.stats)), intersect(col.order, colnames(sqantisim.stats))]
  
  # Add column with match type for plotting later
  data.known$match_type <- 'FN'
  data.known$match_type[which(data.known$transcript_id %in% known.perfect.matches$transcript_id)] <- 'TP'
  data.known$structural_category <- 'known'
  data.novel$match_type <- 'FN'
  data.novel$match_type[which(data.novel$transcript_id %in% novel.perfect.matches$transcript_id)] <- 'TP'
  
  res <- list(data.known, data.novel, known.matches, novel.matches, known.perfect.matches, novel.perfect.matches, sqantisim.stats)
  names(res) <- c("data.known", "data.novel", "known.matches", "novel.matches", "known.perfect.matches", "novel.perfect.matches", "sqantisim.stats")
  return(res)
}

modify_index_file <- function(index.file, res.full){
  modif.index <- read.table(index.file, header=T, as.is=T, sep="\t")
  modif.index$pipeline_performance <- "FN"
  modif.index$pipeline_performance[which(modif.index$sim_counts <= 0)] <-  "absent"
  modif.index$pipeline_performance[which(modif.index$transcript_id %in% res.full$known.matches$transcript_id | modif.index$transcript_id %in% res.full$novel.matches$transcript_id)] <- "PTP"
  modif.index$pipeline_performance[which(modif.index$transcript_id %in% res.full$known.perfect.matches$transcript_id | modif.index$transcript_id %in% res.full$novel.perfect.matches$transcript_id)] <- "TP"
  write.table(modif.index, file = paste0(index.file, ".eval"), quote = F, sep = "\t", na = "NA",row.names = F)
}

#######################################
#                                     #
#                MAIN                 #
#                                     #
#######################################

# -------------------- Argument parser
args <- commandArgs(trailingOnly = TRUE)
class.file <- args[1] # classification file SQANTI3
junc.file <- args[2] # junctions file SQANTI3
index.file <- args[3] # index file
min.supp <- args[4] # min support reads
src.path <- args[5] # path to src utilities

output_directory <- dirname(dirname(class.file))
output_name <- basename(strsplit(class.file, "_classification.txt")[[1]][1])

# -------------------- Read files
# Read classification file
data.class <- read.table(class.file, header=T, as.is=T, sep="\t")
rownames(data.class) <- data.class$isoform
xaxislevelsF1 <- c("full-splice_match","incomplete-splice_match","novel_in_catalog","novel_not_in_catalog", "genic","antisense","fusion","intergenic","genic_intron");
xaxislabelsF1 <- c("FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
data.class$structural_category = factor(data.class$structural_category,
                                        labels = xaxislabelsF1,
                                        levels = xaxislevelsF1,
                                        ordered=TRUE)

# Read juncs file
data.junction <- read.table(junc.file, header=T, as.is=T, sep="\t")
data.junction <- data.junction %>%
  group_by(isoform) %>%
  summarise(Donors=list(as.numeric(genomic_start_coord)), Acceptors=list(as.numeric(genomic_end_coord)))
data.junction$Donors <- lapply(data.junction$Donors, function(x){
  paste(sort(unlist(x)), collapse=',')
})
data.junction$Acceptors <- lapply(data.junction$Acceptors, function(x){
  paste(sort(unlist(x)), collapse=',')
})
data.junction$junctions <- paste(data.junction$Donors, data.junction$Acceptors, sep=',')

# Combine class and junc
data.query <- full_join(data.class, data.junction, by='isoform')
data.query$junctions[which(is.na(data.query$junctions))] <- ''
data.query <- data.query[,c('isoform', 'chrom', 'strand', 'structural_category', 'junctions', 'TSS_genomic_coord', 'TTS_genomic_coord', 'all_canonical', 'dist_to_CAGE_peak', 'within_CAGE_peak', 'min_cov', 'ratio_TSS')]

# Read index file
data.index <- read.table(index.file, header=T, as.is=T, sep="\t")
data.index$donors <- lapply(data.index$donors, function(x){
  paste(sort(unlist(as.numeric(unlist(strsplit(as.character(x), ",")))+1)), collapse=',')
})
data.index$acceptors <- lapply(data.index$acceptors, function(x){
  paste(sort(unlist(as.numeric(unlist(strsplit(as.character(x), ",")))-1)), collapse=',')
})
data.index$junctions <- paste(data.index$donors, data.index$acceptors, sep=',')
data.index$junctions[which(data.index$junctions == ',')] <- ''
data.index[which(data.index$strand == "-"),"TSS_genomic_coord"]  <- data.index[which(data.index$strand == "-"),"TSS_genomic_coord"] - 1
data.index[which(data.index$strand == "-"),"TTS_genomic_coord"]  <- data.index[which(data.index$strand == "-"),"TTS_genomic_coord"] + 1
data.index$structural_category = factor(data.index$structural_category,
                                      labels = xaxislabelsF1,
                                      levels = xaxislevelsF1,
                                      ordered=TRUE)
data.index$donors <- NULL
data.index$acceptors <- NULL
data.index$sim_type[which(data.index$sim_counts == 0)] <- 'absent' # Ignore not simulated

# -------------------- Performance metrics
# Matched for novel and known
MAX_TSS_TTS_DIFF = 50
res.full <- get_performance_metrics(data.query, data.index, MAX_TSS_TTS_DIFF)
res.min <- get_performance_metrics(data.query, data.index, MAX_TSS_TTS_DIFF, min.supp)

modify_index_file(index.file, res.full)

#######################################
#                                     #
#     TABLE AND PLOT GENERATION       #
#                                     #
#######################################

# -------------------- 
# -------------------- 
# TABLE INDEX
# t1: SQANTISIM metrics
# t2: SQANTISIM metrics above min_supp

# -------------------- 
# -------------------- 
# PLOT INDEX
# p1: simulated expression profile
# p2: structural classification
# p3: novel TP vs FN - mono/multi-exon
# p4: canonical juncs
# p5: within cage peak
# p6: distance to cage peak
# p7: min SJ cov by short reads
# p8: ratio TSS
# px: radar chart from perfomance metrics

print("***Generating plots for the report")

# -------------------- Global plot parameters
# SAME FORMAT AS SQANTI3 REPORT
#myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")
myPalette = c("#3A5A81", "#D31336", "#252131", "#6BAED6","#FC8D59","#78C679","#EE6A50")

cat.palette = c("FSM"="#6BAED6", "ISM"="#FC8D59", "NIC"="#78C679", 
                "NNC"="#EE6A50", "Genic\nGenomic"="#969696", "Antisense"="#66C2A4", "Fusion"="goldenrod1",
                "Intergenic" = "darksalmon", "Genic\nIntron"="#41B6C4")

mytheme <- theme_classic(base_family = "Helvetica") +
  theme(axis.line.x = element_line(color="black", size = 0.4),
        axis.line.y = element_line(color="black", size = 0.4)) +
  theme(axis.title.x = element_text(size=13),
        axis.text.x  = element_text(size=12),
        axis.title.y = element_text(size=13),
        axis.text.y  = element_text(vjust=0.5, size=12) ) +
  theme(legend.text = element_text(size = 11), legend.title = element_text(size=12), legend.key.size = unit(0.5, "cm")) +
  theme(plot.title = element_text(lineheight=.4, size=15, hjust = 0.5)) +
  theme(plot.margin = unit(c(2.5,1,1,1), "cm"))

# -------------------- 
# TABLE 1: SQANTISIM metrics
t1 <- DT::datatable(res.full$sqantisim.stats, class = 'compact', options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(res.full$sqantisim.stats), digits = 3, rows=c(6:11), zero.print = 0)

# TABLE 2: SQANTISIM metrics above min_supp
t2 <- DT::datatable(res.min$sqantisim.stats, class = 'compact', options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(res.min$sqantisim.stats), digits = 3, rows=c(6:11), zero.print = 0)

# -------------------- PLOT FULL
# PLOT 1: simulated expression profile
expr.dist <- data.index[which(data.index$sim_type %in% c('novel', 'known')), c('sim_type', 'sim_counts')]

p1 <- expr.dist %>%
  ggplot(aes(x=sim_counts, fill=sim_type)) +
  geom_histogram(aes(y=stat(count)/sum(count)), color="black", alpha=0.5, position = 'identity', bins = round(sqrt(nrow(expr.dist)))) +
  mytheme +
  scale_fill_manual(values = myPalette[1:2], name='Transcript type') +
  xlab('Number of simulated reads') + 
  ylab('Percentage of transcripts') +
  ggtitle('Simulated expression levels')

# PLOT 2: structural classification
p2 <- data.class %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, fill=structural_category), color="black", size=0.3, width=0.7) +
  scale_x_discrete(drop=FALSE) + 
  xlab('') + 
  ylab('Transcripts %') +
  mytheme +
  geom_blank(aes(y=((..count..)/sum(..count..))), stat = "count") +
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_manual(values = cat.palette, guide='none') +
  ggtitle("Isoform Distribution Across Structural Categories\n\n" ) +
  theme(axis.title.x=element_blank()) +  theme(axis.text.x  = element_text(margin=ggplot2::margin(17,0,0,0), size=12)) +
  theme(legend.justification=c(1,1), legend.position=c(1,1))

# PLOT 3: novel TP vs FN - mono/multi-exon
p3 <- rbind(res.full$data.novel, res.full$data.known) %>%
  mutate(exon_type=ifelse(exons > 1, 'multi-exon', 'mono-exon')) %>%
  group_by(structural_category, match_type, exon_type) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=exon_type), position="fill", stat="identity") +
  scale_fill_manual(values=myPalette[1:2], name='Stats') +
  scale_alpha_manual(values=c(0.5,1), name='Exons') +
  mytheme +
  ylab('Percentage %') +
  xlab('')+
  ggtitle('Single- and Multi-exon identifications')

# PLOT 4: canonical juncs
data.query$match_type <- 'FP'
data.query$match_type[which(data.query$isoform %in% res.full$novel.perfect.matches$isoform)] <- 'TP'
data.query$match_type[which(data.query$isoform %in% res.full$known.perfect.matches$isoform)] <- 'TP'
p4 <- data.query[which(!is.na(data.query$all_canonical)),] %>%
  group_by(structural_category, match_type, all_canonical) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=all_canonical), position="fill", stat="identity") +
  scale_fill_manual(values=myPalette[1:2], name='Stats') +
  scale_alpha_manual(values=c(1, 0.5), name='Junctions') +
  mytheme +
  ylab('Percentage %') +
  xlab('') +
  ggtitle('Canonical or Non Canonical Junctions') +
  theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))

if ('within_CAGE_peak' %in% colnames(data.index)){
  # PLOT 5: within cage peak
  data.query$match_type <- 'FP'
  data.query$match_type[which(data.query$isoform %in% res.full$novel.perfect.matches$isoform)] <- 'novel_TP'
  data.query$match_type[which(data.query$isoform %in% res.full$known.perfect.matches$isoform)] <- 'known_TP'
  p5.known_FN <- data.index[which(data.index$transcript_id %in% res.full$data.known$transcript_id & !(data.index$transcript_id %in% res.full$known.perfect.matches$transcript_id)),]
  p5.known_FN$match_type <- 'known_FN'
  p5.novel_FN <- data.index[which(data.index$transcript_id %in% res.full$data.novel$transcript_id & !(data.index$transcript_id %in% res.full$novel.perfect.matches$transcript_id)),]
  p5.novel_FN$match_type <- 'novel_FN'
  p5.all <- rbind(data.query[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')],
                  p5.known_FN[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')],
                  p5.novel_FN[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')])
  
  p5 <- p5.all[which(!is.na(p5.all$within_CAGE_peak)),] %>%
    group_by(match_type, within_CAGE_peak) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=match_type)) +
    geom_bar(aes(y=value, fill=within_CAGE_peak), position="fill", stat="identity") +
    scale_fill_manual(values=myPalette[1:2], name='CagePeak') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('Within CAGE peak') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
  
  # PLOT 6: distance to cage peak
  p6 <- p5.all[which(!is.na(p5.all$dist_to_CAGE_peak)),] %>%
    ggplot(aes(x=dist_to_CAGE_peak, color=match_type, fill=match_type)) +
    geom_density(alpha=.3) +
    scale_color_manual(values = myPalette) +
    scale_fill_manual(values = myPalette) +
    mytheme +
    ylab('Distance to CAGE peak') +
    xlab('') +
    ggtitle('Distance To Cage Peak') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}

if ('min_cov' %in% colnames(data.index)) {
  # PLOT 7: min SJ cov by short reads
  p7.all <- rbind(data.query[,c('structural_category', 'match_type', 'min_cov')],
                  p5.known_FN[,c('structural_category', 'match_type', 'min_cov')],
                  p5.novel_FN[,c('structural_category', 'match_type', 'min_cov')])
  p7.all$Coverage_SJ <- 'False'
  p7.all[which(p7.all$min_cov>0), 'Coverage_SJ'] <- 'True'
  
  p7 <- p7.all[which(!is.na(p7.all$Coverage_SJ)),] %>%
    group_by(match_type, Coverage_SJ) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=match_type)) +
    geom_bar(aes(y=value, fill=Coverage_SJ), position="fill", stat="identity") +
    scale_fill_manual(values=myPalette[1:2], name='Coverage SJ') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('Splice Junctions Short Reads Coverage') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
  
  # PLOT 8: ratio TSS
  p8.all <- rbind(data.query[,c('structural_category', 'match_type', 'ratio_TSS')],
                  p5.known_FN[,c('structural_category', 'match_type', 'ratio_TSS')],
                  p5.novel_FN[,c('structural_category', 'match_type', 'ratio_TSS')])
  
  p8 <- p8.all[which(!is.na(p8.all$ratio_TSS)),] %>%
    ggplot(aes(x=log(ratio_TSS), color=match_type, fill=match_type)) +
    geom_density(alpha=.3) +
    scale_color_manual(values = myPalette) +
    scale_fill_manual(values = myPalette) +
    mytheme +
    ylab('log TSS ratio') +
    xlab('') +
    ggtitle('Ratio TSS') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}

# PLOT X: Radar chart
# Generated in RMD file

# ----------------- PLOT WITH MIN SUPPORT
p3.min <- rbind(res.min$data.novel, res.min$data.known) %>%
  mutate(exon_type=ifelse(exons > 1, 'multi-exon', 'mono-exon')) %>%
  group_by(structural_category, match_type, exon_type) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=exon_type), position="fill", stat="identity") +
  scale_fill_manual(values=myPalette[1:2], name='Stats') +
  scale_alpha_manual(values=c(0.5,1), name='Exons') +
  mytheme +
  ylab('Percentage %') +
  xlab('')+
  ggtitle('Single- and Multi-exon identifications (min. supp.)')

data.query$match_type <- 'FP'
data.query$match_type[which(data.query$isoform %in% res.min$novel.perfect.matches$isoform)] <- 'TP'
data.query$match_type[which(data.query$isoform %in% res.min$known.perfect.matches$isoform)] <- 'TP'
p4.min <- data.query[which(!is.na(data.query$all_canonical)),] %>%
  group_by(structural_category, match_type, all_canonical) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=all_canonical), position="fill", stat="identity") +
  scale_fill_manual(values=myPalette[1:2], name='Stats') +
  scale_alpha_manual(values=c(1, 0.5), name='Junctions') +
  mytheme +
  ylab('Percentage %') +
  xlab('') +
  ggtitle('Canonical or Non Canonical Junctions (min. supp.)') +
  theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))

if ('within_CAGE_peak' %in% colnames(data.index)){
  # PLOT 5: within cage peak
  data.query$match_type <- 'FP'
  data.query$match_type[which(data.query$isoform %in% res.min$novel.perfect.matches$isoform)] <- 'novel_TP'
  data.query$match_type[which(data.query$isoform %in% res.min$known.perfect.matches$isoform)] <- 'known_TP'
  p5.known_FN <- data.index[which(data.index$transcript_id %in% res.min$data.known$transcript_id & !(data.index$transcript_id %in% res.min$known.perfect.matches$transcript_id)),]
  p5.known_FN$match_type <- 'known_FN'
  p5.novel_FN <- data.index[which(data.index$transcript_id %in% res.min$data.novel$transcript_id & !(data.index$transcript_id %in% res.min$novel.perfect.matches$transcript_id)),]
  p5.novel_FN$match_type <- 'novel_FN'
  p5.all <- rbind(data.query[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')],
                  p5.known_FN[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')],
                  p5.novel_FN[,c('structural_category', 'match_type', 'within_CAGE_peak', 'dist_to_CAGE_peak')])
  
  p5.min <- p5.all[which(!is.na(p5.all$within_CAGE_peak)),] %>%
    group_by(match_type, within_CAGE_peak) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=match_type)) +
    geom_bar(aes(y=value, fill=within_CAGE_peak), position="fill", stat="identity") +
    scale_fill_manual(values=myPalette[1:2], name='CagePeak') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('Within cage peak (min. supp.)') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
  
  # PLOT 6: distance to cage peak
  p6.min <- p5.all[which(!is.na(p5.all$dist_to_CAGE_peak)),] %>%
    ggplot(aes(x=dist_to_CAGE_peak, color=match_type, fill=match_type)) +
    geom_density(alpha=.3) +
    scale_color_manual(values = myPalette) +
    scale_fill_manual(values = myPalette) +
    mytheme +
    ylab('Distance to CAGE peak') +
    xlab('') +
    ggtitle('Distance To Cage Peak (min. supp.)') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}

if ('min_cov' %in% colnames(data.index)) {
  # PLOT 7: min SJ cov by short reads
  p7.all <- rbind(data.query[,c('structural_category', 'match_type', 'min_cov')],
                  p5.known_FN[,c('structural_category', 'match_type', 'min_cov')],
                  p5.novel_FN[,c('structural_category', 'match_type', 'min_cov')])
  p7.all$Coverage_SJ <- 'False'
  p7.all[which(p7.all$min_cov>0), 'Coverage_SJ'] <- 'True'
  
  p7.min <- p7.all[which(!is.na(p7.all$Coverage_SJ)),] %>%
    group_by(match_type, Coverage_SJ) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=match_type)) +
    geom_bar(aes(y=value, fill=Coverage_SJ), position="fill", stat="identity") +
    scale_fill_manual(values=myPalette[1:2], name='Coverage SJ') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('Splice Junctions Short Reads Coverage (min. supp.)') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
  
  p8.all <- rbind(data.query[,c('structural_category', 'match_type', 'ratio_TSS')],
                  p5.known_FN[,c('structural_category', 'match_type', 'ratio_TSS')],
                  p5.novel_FN[,c('structural_category', 'match_type', 'ratio_TSS')])
  
  p8.min <- p8.all[which(!is.na(p8.all$ratio_TSS)),] %>%
    ggplot(aes(x=log(ratio_TSS), color=match_type, fill=match_type)) +
    geom_density(alpha=.3) +
    scale_color_manual(values = myPalette) +
    scale_fill_manual(values = myPalette) +
    mytheme +
    ylab('log TSS ratio') +
    xlab('') +
    ggtitle('Ratio TSS (min. supp.)') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}

# -------------------- Output report
rmarkdown::render(
  input = paste(src.path, 'SQANTISIM_report.Rmd', sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, "_SQANTISIM_report.html")
)

