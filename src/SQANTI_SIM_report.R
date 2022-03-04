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
suppressMessages(library(gridExtra))
suppressMessages(library(knitr))
suppressMessages(library(rmarkdown))
suppressMessages(library(ggplot2))


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

class.file <- 'sqanti_sim_classification.txt'
junc.file <- 'sqanti_sim_junctions.txt'
index.file <- 'mix_index.tsv'


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
data.query <- data.query[,c('isoform', 'strand', 'structural_category', 'junctions', 'TSS_genomic_coord', 'TTS_genomic_coord')]

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

novel.matches <- inner_join(data.query, data.novel, by='junctions') %>%
  mutate(diffTSS = abs(TSS_genomic_coord.x - TSS_genomic_coord.y), diffTTS = abs(TTS_genomic_coord.x - TTS_genomic_coord.y), difftot = diffTSS+diffTTS) %>%
  arrange(difftot) %>%
  distinct(isoform, .keep_all = T)
novel.perfect.matches <- novel.matches[which(novel.matches$diffTSS < 50 & novel.matches$diffTTS < 50),]

known.metrics <- data.frame(init=c())
novel.metrics <- data.frame(init=c())
for (sc in xaxislabelsF1){
  known.TP <- nrow(known.perfect.matches[which(known.perfect.matches$structural_category.x == sc),])
  known.PTP <- nrow(known.matches[which(known.matches$structural_category.x == sc),]) - known.TP
  
  if (sc %in% sim.sc) {
    novel.TP <- nrow(novel.perfect.matches[which(novel.perfect.matches$structural_category.x == sc),])
    novel.PTP <- nrow(novel.matches[which(novel.matches$structural_category.x == sc),]) - novel.TP
    
    FP <- nrow(data.query[which(data.query$structural_category == sc),]) - known.TP - novel.TP
    novel.FN <- nrow(data.novel[which(data.novel$structural_category == sc),]) - novel.TP
    
    
    novel.metrics['TP', sc] <- novel.TP
    novel.metrics['PTP', sc] <- novel.PTP
    novel.metrics['FP', sc] <- FP
    novel.metrics['FN', sc] <- novel.FN
    novel.metrics['Sensitivity', sc] <- novel.TP/ (novel.TP + novel.FN)
    novel.metrics['Precision', sc] <- novel.TP/ (novel.TP + FP)
    
  } else {
    FP <- nrow(data.query[which(data.query$structural_category == sc),]) - known.TP
  }
  
  known.metrics['TP', sc] <- known.TP
  known.metrics['PTP', sc] <- known.PTP
  known.metrics['FP', sc] <- FP
  known.metrics['Precision', sc] <- known.TP/ (known.TP + FP)
  
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



known.metrics['TP', 'global'] <- known.TP
known.metrics['PTP', 'global'] <- known.PTP
known.metrics['FP', 'global'] <- known.FP
known.metrics['FN', 'global'] <- nrow(data.known) - known.TP
known.metrics['Sensitivity', 'global'] <- known.TP / (known.TP + known.FP)
known.metrics['Precision', 'global'] <- known.TP / (known.TP + known.metrics['FN', 'global'])
known.metrics['F-score', 'global'] <- 2*((known.metrics['Sensitivity', 'global']*known.metrics['Precision', 'global'])/(known.metrics['Sensitivity', 'global']+known.metrics['Precision', 'global']))

novel.metrics['TP', 'global'] <- novel.TP
novel.metrics['PTP', 'global'] <- novel.PTP
novel.metrics['FP', 'global'] <- novel.FP
novel.metrics['FN', 'global'] <- nrow(data.novel) - novel.TP
novel.metrics['Sensitivity', 'global'] <- novel.TP / (novel.TP + novel.FP)
novel.metrics['Precision', 'global'] <- novel.TP / (novel.TP + novel.metrics['FN', 'global'])
novel.metrics['F-score', 'global'] <- 2*((novel.metrics['Sensitivity', 'global']*novel.metrics['Precision', 'global'])/(novel.metrics['Sensitivity', 'global']+novel.metrics['Precision', 'global']))


col.order <- c("global", "FSM", "ISM", "NIC", "NNC", "Genic\nGenomic",  "Antisense", "Fusion","Intergenic", "Genic\nIntron")
row.order <- c('TP', 'PTP', 'FP', 'FN', 'Sensitivity', 'Precision', 'F-score')
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
# p3: 

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
t1 <- DT::datatable(known.metrics) %>%
  formatRound(colnames(known.metrics), digits = 3, rows=c(5,6,7), zero.print = 0)

# TABLE 2: novel metrics
t2 <- DT::datatable(novel.metrics) %>%
  formatRound(colnames(novel.metrics), digits = 3, rows=c(5,6,7), zero.print = 0)


# -------------------- 
# PLOT 1: simulated expression profile
expr.dist <- data.index[which(data.index$sim_type %in% c('novel', 'known')), c('sim_type', 'sim_counts')]

p1 <- expr.dist %>%
  ggplot( aes(x=sim_counts, fill=sim_type)) +
  geom_histogram(color='white', alpha=0.6, position = 'identity') +
  mytheme +
  scale_fill_manual(values = c('orange', 'darkcyan')) +
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

# PLOT 3: novel precision vs sensitivity



# -------------------- Output report
rmarkdown::render(
  input = paste(src.path, 'SQANTI_SIM_report.Rmd', sep = "/"),
  output_dir = output_directory,
  output_file = paste0(output_name, "_SQANTI_SIM_report.html")
)

