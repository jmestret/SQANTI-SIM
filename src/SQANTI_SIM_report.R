#######################################
#                                     #
#    SQANTI-SIM report generation     #
#                                     #
#######################################

# Author: Jorge Meste
# Last modified: 03/02/2022 by Jorge Mestre


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
src.path <- args[4] # path to src utilities

output_directory <- dirname(class.file)
output_name <- basename(strsplit(class.file, "_classification.txt")[[1]][1])

class.file <-'talon_cage_classification.txt'
junc.file <-'talon_cage_junctions.txt'
index.file <- 'mix_equal_index.tsv'
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


data.query <- full_join(data.class, data.junction, by='isoform')
data.query$junctions[which(is.na(data.query$junctions))] <- ''
data.query <- data.query[,c('isoform', 'strand', 'structural_category','all_canonical', 'dist_to_cage_peak', 'within_cage_peak', 'junctions', 'TSS_genomic_coord', 'TTS_genomic_coord')]

# Read deleted file
data.index <- read.table(index.file, header=T, as.is=T, sep="\t")
data.index$donors <- lapply(data.index$donors, function(x){
  paste(sort(unlist(as.numeric(unlist(strsplit(as.character(x), ",")))+1)), collapse=',')
})
data.index$acceptors <- lapply(data.index$acceptors, function(x){
  paste(sort(unlist(as.numeric(unlist(strsplit(as.character(x), ",")))-1)), collapse=',')
})
data.index$junctions <- paste(data.index$donors, data.index$acceptors, sep=',')
data.index$junctions[which(data.index$junctions == ',')] <- ''
data.index$TSS_genomic_coord  <- data.index$TSS_genomic_coord - 1
data.index$structural_category = factor(data.index$structural_category,
                                      labels = xaxislabelsF1,
                                      levels = xaxislevelsF1,
                                      ordered=TRUE)
data.index$donors <- NULL
data.index$acceptors <- NULL
data.index$sim_type[which(data.index$sim_counts == 0)] <- 'absent'

# Transcript simulated
data.novel <- data.index[which(data.index$sim_type == 'novel'),]
data.known <- data.index[which(data.index$sim_type == 'known'),]
sim.sc <- unique(data.novel$structural_category)

# Matched for novel and known
known.matches <- inner_join(data.query, data.known, by='junctions') %>%
  mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
  arrange(difftot) %>%
  distinct(isoform, .keep_all = T)
known.perfect.matches <- known.matches[which(known.matches$diffTSS < 50 & known.matches$diffTTS < 50),]
cond <- (known.perfect.matches$exons > 1) | (known.perfect.matches$strand.x == '+' & known.perfect.matches$TSS_genomic_coord.x <= known.perfect.matches$TTS_genomic_coord.y & known.perfect.matches$TSS_genomic_coord.y <= known.perfect.matches$TTS_genomic_coord.x) | (known.perfect.matches$strand.x == '-' & known.perfect.matches$TTS_genomic_coord.x <= known.perfect.matches$TSS_genomic_coord.y & known.perfect.matches$TTS_genomic_coord.y <= known.perfect.matches$TSS_genomic_coord.x)
known.perfect.matches <- known.perfect.matches[cond,]

novel.matches <- inner_join(data.query, data.novel, by='junctions') %>%
  mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
  arrange(difftot) %>%
  distinct(isoform, .keep_all = T)
novel.perfect.matches <- novel.matches[which(novel.matches$diffTSS < 50 & novel.matches$diffTTS < 50),]
cond <- (novel.perfect.matches$exons > 1) | (novel.perfect.matches$strand.x == '+' & novel.perfect.matches$TSS_genomic_coord.x <= novel.perfect.matches$TTS_genomic_coord.y & novel.perfect.matches$TSS_genomic_coord.y <= novel.perfect.matches$TTS_genomic_coord.x) | (novel.perfect.matches$strand.x == '-' & novel.perfect.matches$TTS_genomic_coord.x <= novel.perfect.matches$TSS_genomic_coord.y & novel.perfect.matches$TTS_genomic_coord.y <= novel.perfect.matches$TSS_genomic_coord.x)
novel.perfect.matches <- novel.perfect.matches[cond,]

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
known.metrics['Sensitivity', 'global'] <- known.TP / (known.TP + known.FP)
known.metrics['Precision', 'global'] <- known.TP / (known.TP + known.metrics['FN', 'global'])
known.metrics['F-score', 'global'] <- 2*((known.metrics['Sensitivity', 'global']*known.metrics['Precision', 'global'])/(known.metrics['Sensitivity', 'global']+known.metrics['Precision', 'global']))
known.metrics['False_Discovery_Rate', 'global'] <- (FP + known.PTP) / (FP + known.PTP +  known.TP)
known.metrics['Positive_Detection_Rate', 'global'] <- (known.TP + known.PTP) / nrow(data.known)
known.metrics['False_Detection_Rate', 'global'] <- (known.FP) / (known.FP + known.PTP +  known.TP)

novel.metrics['Total', 'global'] <-  nrow(data.novel)
novel.metrics['TP', 'global'] <- novel.TP
novel.metrics['PTP', 'global'] <- novel.PTP
novel.metrics['FP', 'global'] <- novel.FP
novel.metrics['FN', 'global'] <- nrow(data.novel) - novel.TP
novel.metrics['Sensitivity', 'global'] <- novel.TP / (novel.TP + novel.FP)
novel.metrics['Precision', 'global'] <- novel.TP / (novel.TP + novel.metrics['FN', 'global'])
novel.metrics['F-score', 'global'] <- 2*((novel.metrics['Sensitivity', 'global']*novel.metrics['Precision', 'global'])/(novel.metrics['Sensitivity', 'global']+novel.metrics['Precision', 'global']))
novel.metrics['False_Discovery_Rate', 'global'] <- (FP + novel.PTP) / (FP + novel.PTP +  novel.TP)
novel.metrics['Positive_Detection_Rate', 'global'] <- (novel.TP + novel.PTP) / nrow(data.novel)
novel.metrics['False_Detection_Rate', 'global'] <- (novel.FP) / (novel.FP + novel.PTP +  novel.TP)

col.order <- c("global", "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
row.order <- c('Total', 'TP', 'PTP', 'FP', 'FN', 'Sensitivity', 'Precision', 'F-score', 'False_Discovery_Rate', 'Positive_Detection_Rate', 'False_Detection_Rate')
known.metrics <- known.metrics[intersect(row.order, rownames(known.metrics)), intersect(col.order, colnames(known.metrics))]
novel.metrics <- novel.metrics[intersect(row.order, rownames(novel.metrics)), intersect(col.order, colnames(novel.metrics))]


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
# p3: sensitivity by trans per gene
#   p3.1: discrete
#   p3.2: factor
# p4: novel TP vs FN
# p6: radar chart

print("***Generating plots for the report")

# -------------------- Global plot parameters
# SAME FORMAT AS SQANTI3 REPORT

myPalette = c("#6BAED6","#FC8D59","#78C679","#EE6A50","#969696","#66C2A4", "goldenrod1", "darksalmon", "#41B6C4","tomato3", "#FE9929")

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
t1 <- DT::datatable(known.metrics, class = 'compact', options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(known.metrics), digits = 3, rows=c(6:11), zero.print = 0)

# TABLE 2: novel metrics
t2 <- DT::datatable(novel.metrics, class = 'compact',  options = list(pageLength = 15, dom = 'tip')) %>%
  formatRound(colnames(novel.metrics), digits = 3, rows=c(6:11), zero.print = 0)


# -------------------- 
# PLOT 1: simulated expression profile
expr.dist <- data.index[which(data.index$sim_type %in% c('novel', 'known')), c('sim_type', 'sim_counts')]

p1 <- expr.dist %>%
  ggplot( aes(x=sim_counts, fill=sim_type)) +
  geom_histogram(color='white', alpha=0.5, position = 'identity') +
  mytheme +
  scale_fill_manual(values = c('orange', 'darkcyan'), name='Simulation type') +
  xlab('Counts') + 
  ylab('Number of transcripts') +
  ggtitle('Simulated counts distribution')

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

# PLOT 3: sensitivity by trans per gene
trans.per.gene <- data.index %>%
  group_by(gene_id) %>%
  summarise(n_trans = n())
trans.per.gene <- right_join(trans.per.gene, data.novel, by='gene_id')
trans.per.gene$found <- sapply(trans.per.gene$transcript_id, function(x){
  if (x %in% novel.perfect.matches$transcript_id){
    return('TP')
  } else {return('FN')}
})
trans.per.gene <- trans.per.gene %>%
  group_by(n_trans, found) %>%
  summarise(counts=n()) %>%
  pivot_wider(names_from = found, values_from = counts, values_fill = 0) %>%
  mutate(sensitivity=(TP/(TP+FN)))

p3.1 <- trans.per.gene %>%
  ggplot(aes(x=n_trans, y=sensitivity,)) +
  geom_point(color='orange', size=3) +
  mytheme +
  ylab('Sensitivity') +
  xlab('Number of annotated transcripts per gene')

trans.per.gene$cat <- sapply(trans.per.gene$n_trans, function(x){
  if (x <= 2){
    return('1-2')
  } else if (3 <= x & x <= 10){
    return('3-10')
  } else if (10 < x & x <= 20){
    return('11-20')
  } else {
      return('>20')
    }
})

trans.per.gene$cat <- factor(trans.per.gene$cat, levels=c('1-2', '3-10','11-20', '>20'))
trans.per.gene <- trans.per.gene %>%
  group_by(cat) %>%
  summarise(TP=sum(TP), FN=sum(FN)) %>%
  mutate(sensitivity=(TP/(TP+FN)))

p3.2 <- trans.per.gene %>%
  ggplot(aes(x=cat, y=sensitivity)) +
  geom_bar(stat='identity', color='cyan', fill='darkcyan') +
  mytheme +
  ylab('Sensitivity') +
  xlab('Number of annotated transcripts per gene')

# PLOT 4: novel TP vs FN
data.known$match_type <- 'FN'
data.known$match_type[which(data.known$transcript_id %in% known.perfect.matches$transcript_id)] <- 'TP'
data.known$structural_category <- 'known'
data.novel$match_type <- 'FN'
data.novel$match_type[which(data.novel$transcript_id %in% novel.perfect.matches$transcript_id)] <- 'TP'
p4 <- rbind(data.novel, data.known) %>%
  mutate(exon_type=ifelse(exons > 1, 'multi-exon', 'mono-exon')) %>%
  group_by(structural_category, match_type, exon_type) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=exon_type), position="fill", stat="identity") +
  scale_fill_manual(values=c('orange', 'darkcyan'), name='Stats') +
  scale_alpha_manual(values=c(0.5,1), name='Exons') +
  mytheme +
  ylab('Percentage %') +
  xlab('')+
  ggtitle('TP vs FN - monoexon vs multiexon')

# PLOT 5: canonical juncs

data.query$match_type <- 'FP'
data.query$match_type[which(data.query$isoform %in% novel.perfect.matches$isoform)] <- 'TP'
data.query$match_type[which(data.query$isoform %in% known.perfect.matches$isoform)] <- 'TP'
p5 <- data.query[which(!is.na(data.query$all_canonical)),] %>%
  group_by(structural_category, match_type, all_canonical) %>%
  summarise(value=n()) %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(fill=match_type, y=value, alpha=all_canonical), position="fill", stat="identity") +
  scale_fill_manual(values=c('orange', 'darkcyan'), name='Stats') +
  scale_alpha_manual(values=c(1, 0.5), name='Junctions') +
  mytheme +
  ylab('Percentage %') +
  xlab('') +
  ggtitle('TP vs FP - canonical junctions') +
  theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))

if (any(!is.na(data.class$within_cage_peak))){
  # PLOT 6: within cage peak
  p6 <- data.query[which(!is.na(data.query$within_cage_peak)),] %>%
    group_by(structural_category, match_type, within_cage_peak) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=structural_category)) +
    geom_bar(aes(fill=match_type, y=value, alpha=within_cage_peak), position="fill", stat="identity") +
    scale_fill_manual(values=c('orange', 'darkcyan'), name='Stats') +
    scale_alpha_manual(values=c(0.5, 1), name='Cage Peak') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('TP vs FP - within cage peak') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
  
  data.query$match_type <- 'FP'
  data.query$match_type[which(data.query$isoform %in% novel.perfect.matches$isoform)] <- 'novel_TP'
  data.query$match_type[which(data.query$isoform %in% known.perfect.matches$isoform)] <- 'known_TP'
  p7 <- data.query[which(!is.na(data.query$within_cage_peak)),] %>%
    group_by(match_type, within_cage_peak) %>%
    summarise(value=n()) %>%
    ggplot(aes(x=match_type)) +
    geom_bar(aes(y=value, fill=within_cage_peak), position="fill", stat="identity") +
    scale_fill_manual(values=c('orange', 'darkcyan'), name='CagePeak') +
    mytheme +
    ylab('Percentage %') +
    xlab('') +
    ggtitle('TP vs FP - within cage peak') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
  
  # PLOT 8: distance to cage peak
  p8 <- data.query[which(!is.na(data.query$dist_to_cage_peak)),] %>%
    ggplot(aes(x=dist_to_cage_peak, color=match_type, fill=match_type)) +
    geom_density(alpha=.6) +
    mytheme +
    ylab('Distance to Cage Peak') +
    xlab('') +
    ggtitle('TP vs FP - distance to cage peak') +
    theme(axis.text.x = element_text(angle = 45, margin=ggplot2::margin(17,0,0,0), size=10))
}



# PLOT 6: Radar chart
# In RMD file

# -------------------- Output report
rmarkdown::render(
  input = paste(src.path, 'SQANTI_SIM_report.Rmd', sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, "_SQANTI_SIM_report.html")
)