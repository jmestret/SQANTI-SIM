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

get_performance_metrics <- function(data.query, data.known, MAX_TSS_TTS_DIFF, min_supp=0){
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
  known.matches <- inner_join(data.query, data.known, by=c('junctions', 'chrom')) %>%
    mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
    arrange(difftot) %>%
    distinct(isoform, .keep_all = T)
  known.perfect.matches <- known.matches[which(known.matches$diffTSS < MAX_TSS_TTS_DIFF & known.matches$diffTTS < MAX_TSS_TTS_DIFF),]
  cond <- (known.perfect.matches$exons > 1) | (known.perfect.matches$strand.x == '+' & known.perfect.matches$TSS_genomic_coord.x <= known.perfect.matches$TTS_genomic_coord.y & known.perfect.matches$TSS_genomic_coord.y <= known.perfect.matches$TTS_genomic_coord.x) | (known.perfect.matches$strand.x == '-' & known.perfect.matches$TTS_genomic_coord.x <= known.perfect.matches$TSS_genomic_coord.y & known.perfect.matches$TTS_genomic_coord.y <= known.perfect.matches$TSS_genomic_coord.x)
  known.perfect.matches <- known.perfect.matches[cond,]
  
  novel.matches <- inner_join(data.query, data.novel, by=c('junctions', 'chrom')) %>%
    mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
    arrange(difftot) %>%
    distinct(isoform, .keep_all = T)
  novel.perfect.matches <- novel.matches[which(novel.matches$diffTSS < MAX_TSS_TTS_DIFF & novel.matches$diffTTS < MAX_TSS_TTS_DIFF),]
  cond <- (novel.perfect.matches$exons > 1) | (novel.perfect.matches$strand.x == '+' & novel.perfect.matches$TSS_genomic_coord.x <= novel.perfect.matches$TTS_genomic_coord.y & novel.perfect.matches$TSS_genomic_coord.y <= novel.perfect.matches$TTS_genomic_coord.x) | (novel.perfect.matches$strand.x == '-' & novel.perfect.matches$TTS_genomic_coord.x <= novel.perfect.matches$TSS_genomic_coord.y & novel.perfect.matches$TTS_genomic_coord.y <= novel.perfect.matches$TSS_genomic_coord.x)
  novel.perfect.matches <- novel.perfect.matches[cond,]
  
  #  Compute metrics
  known.metrics <- data.frame(init=c())
  novel.metrics <- data.frame(init=c())
  for (sc in xaxislabelsF1){
    
    if (sc == 'FSM') {
      known.sim <- nrow(data.known)
    } else {
      known.sim <- 0
    }
    
    known.TP <- nrow(known.perfect.matches[which(known.perfect.matches$structural_category.x == sc),])
    known.PTP <- nrow(known.matches[which(known.matches$structural_category.x == sc),]) - known.TP
    known.FN <- known.sim - known.TP
    
    if (sc %in% sim.sc) {
      novel.sim <- nrow(data.novel[which(data.novel$structural_category == sc),])
      novel.TP <- nrow(novel.perfect.matches[which(novel.perfect.matches$structural_category.x == sc),])
      novel.PTP <- nrow(novel.matches[which(novel.matches$structural_category.x == sc),]) - novel.TP
      
      FP <- nrow(data.query[which(data.query$structural_category == sc),]) - known.TP - known.PTP - novel.TP - novel.PTP
      novel.FN <- novel.sim - novel.TP
      
      
      novel.metrics['Total', sc] <- novel.sim
      novel.metrics['TP', sc] <- novel.TP
      novel.metrics['PTP', sc] <- novel.PTP
      novel.metrics['FP', sc] <- FP
      novel.metrics['FN', sc] <- novel.FN
      novel.metrics['Sensitivity', sc] <- novel.TP/ (novel.TP + novel.FN)
      novel.metrics['Precision', sc] <- novel.TP/ (novel.TP + FP)
      novel.metrics['F-score', sc] <- 2*((novel.metrics['Sensitivity', sc]*novel.metrics['Precision', sc])/(novel.metrics['Sensitivity', sc]+novel.metrics['Precision', sc]))
      novel.metrics['False_Discovery_Rate', sc] <- (FP + novel.PTP) / (FP + novel.PTP +  novel.TP)
      novel.metrics['Positive_Detection_Rate', sc] <- (novel.TP + novel.PTP) / novel.sim
      novel.metrics['False_Detection_Rate', sc] <- (FP) / (FP + novel.PTP +  novel.TP)
      
    } else {
      FP <- nrow(data.query[which(data.query$structural_category == sc),]) - known.TP - known.PTP
    }
    
    known.metrics['Total', sc] <- known.sim
    known.metrics['TP', sc] <- known.TP
    known.metrics['PTP', sc] <- known.PTP
    known.metrics['FP', sc] <- FP
    known.metrics['FN', sc] <- known.FN
    known.metrics['Sensitivity', sc] <- known.TP/ (known.TP + known.FN)
    known.metrics['Precision', sc] <- known.TP/ (known.TP + FP)
    known.metrics['F-score', sc] <- 2*((known.metrics['Sensitivity', sc]*known.metrics['Precision', sc])/(known.metrics['Sensitivity', sc]+known.metrics['Precision', sc]))
    known.metrics['False_Discovery_Rate', sc] <- (FP + known.PTP) / (FP + known.PTP +  known.TP)
    known.metrics['Positive_Detection_Rate', sc] <- (known.TP + known.PTP) / known.sim
    known.metrics['False_Detection_Rate', sc] <- (FP) / (FP + known.PTP +  known.TP)
    
  }
  
  known.TP <- 0
  known.PTP <- 0
  known.FP <- 0
  known.FN <- 0
  novel.TP <- 0
  novel.PTP <- 0
  novel.FP <- 0
  novel.FN <- 0 
  for (sc in xaxislabelsF1){
    known.TP <- known.TP + known.metrics['TP', sc]
    known.PTP <- known.PTP + known.metrics['PTP', sc]
    known.FP <- known.FP + known.metrics['FP', sc]
    if (sc %in% sim.sc){
      novel.TP <- novel.TP + novel.metrics['TP', sc]
      novel.PTP <- novel.PTP + novel.metrics['PTP', sc]
      novel.FP <- novel.FP + novel.metrics['FP', sc]
    }
  }
  
  known.metrics['Total', 'global'] <-  nrow(data.known)
  known.metrics['TP', 'global'] <- known.TP
  known.metrics['PTP', 'global'] <- known.PTP
  known.metrics['FP', 'global'] <- known.FP
  known.metrics['FN', 'global'] <- nrow(data.known) - known.TP
  known.metrics['Precision', 'global'] <- known.TP / (known.TP + known.FP)
  known.metrics['Sensitivity', 'global'] <- known.TP / (known.TP + known.metrics['FN', 'global'])
  known.metrics['F-score', 'global'] <- 2*((known.metrics['Sensitivity', 'global']*known.metrics['Precision', 'global'])/(known.metrics['Sensitivity', 'global']+known.metrics['Precision', 'global']))
  known.metrics['False_Discovery_Rate', 'global'] <- (FP + known.PTP) / (FP + known.PTP +  known.TP)
  known.metrics['Positive_Detection_Rate', 'global'] <- (known.TP + known.PTP) / nrow(data.known)
  known.metrics['False_Detection_Rate', 'global'] <- (known.FP) / (known.FP + known.PTP +  known.TP)
  
  novel.metrics['Total', 'global'] <-  nrow(data.novel)
  novel.metrics['TP', 'global'] <- novel.TP
  novel.metrics['PTP', 'global'] <- novel.PTP
  novel.metrics['FP', 'global'] <- novel.FP
  novel.metrics['FN', 'global'] <- nrow(data.novel) - novel.TP
  novel.metrics['Precision', 'global'] <- novel.TP / (novel.TP + novel.FP)
  novel.metrics['Sensitivity', 'global'] <- novel.TP / (novel.TP + novel.metrics['FN', 'global'])
  novel.metrics['F-score', 'global'] <- 2*((novel.metrics['Sensitivity', 'global']*novel.metrics['Precision', 'global'])/(novel.metrics['Sensitivity', 'global']+novel.metrics['Precision', 'global']))
  novel.metrics['False_Discovery_Rate', 'global'] <- (FP + novel.PTP) / (FP + novel.PTP +  novel.TP)
  novel.metrics['Positive_Detection_Rate', 'global'] <- (novel.TP + novel.PTP) / nrow(data.novel)
  novel.metrics['False_Detection_Rate', 'global'] <- (novel.FP) / (novel.FP + novel.PTP +  novel.TP)
  
  col.order <- c("global", "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
  row.order <- c('Total', 'TP', 'PTP', 'FP', 'FN', 'Sensitivity', 'Precision', 'F-score', 'False_Discovery_Rate', 'Positive_Detection_Rate', 'False_Detection_Rate')
  known.metrics <- known.metrics[intersect(row.order, rownames(known.metrics)), intersect(col.order, colnames(known.metrics))]
  novel.metrics <- novel.metrics[intersect(row.order, rownames(novel.metrics)), intersect(col.order, colnames(novel.metrics))]
  
  # Add column with match type for plotting later
  data.known$match_type <- 'FN'
  data.known$match_type[which(data.known$transcript_id %in% known.perfect.matches$transcript_id)] <- 'TP'
  data.known$structural_category <- 'known'
  data.novel$match_type <- 'FN'
  data.novel$match_type[which(data.novel$transcript_id %in% novel.perfect.matches$transcript_id)] <- 'TP'
  
  res <- list(data.known, data.novel, known.matches, novel.matches, known.perfect.matches, novel.perfect.matches, known.metrics, novel.metrics)
  names(res) <- c("data.known", "data.novel", "known.matches", "novel.matches", "known.perfect.matches", "novel.perfect.matches", "known.metrics", "novel.metrics")
  return(res)
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

class.file <- "sqantisim_classification.txt"
junc.file <- "sqantisim_junctions.txt"
index.file <- "sqantisim_index.tsv"
min.supp <- 3

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
#data.index$TSS_genomic_coord  <- data.index$TSS_genomic_coord - 1
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

#######################################
#                                     #
#     TABLE AND PLOT GENERATION       #
#                                     #
#######################################

# -------------------- 
# -------------------- 
# TABLE INDEX
# t1: known metrics

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
# TABLE 1: known metrics
t1 <- DT::datatable(res.full$known.metrics, class = 'compact', options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(res.full$known.metrics), digits = 3, rows=c(6:11), zero.print = 0)

t1.min <- DT::datatable(res.min$known.metrics, class = 'compact', options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(res.min$known.metrics), digits = 3, rows=c(6:11), zero.print = 0)

# TABLE 2: novel metrics
t2 <- DT::datatable(res.full$novel.metrics, class = 'compact',  options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(res.full$novel.metrics), digits = 3, rows=c(6:11), zero.print = 0)

t2.min <- DT::datatable(res.min$novel.metrics, class = 'compact',  options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(res.min$novel.metrics), digits = 3, rows=c(6:11), zero.print = 0)

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

