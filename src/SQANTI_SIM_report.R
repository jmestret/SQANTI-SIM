#######################################
#                                     #
#    SQANTI-SIM report generation     #
#                                     #
#######################################

# Author: Jorge Meste
# Last modified: 28/02/2022 by Jorge Mestre


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

setClass('QueryIsoform', slots = list(
  id='character', gene_id='character', str_class='character',
  junctions='vector', start='numeric', end='numeric'
))

calc_known_stats <- function(sc, data.query, data.known){
  matches <- apply(data.query, 1, function(x, data.known){
    for (row in 1:nrow(data.known)){
      ref_donor <- data.known[row,'Donors']
      ref_acceptor <- data.known[row,'Acceptors']
      ref_start <- data.known[row,'TSS']
      ref_end <- data.known[row,'TTS']
      
      if (setequal(x$Donors, ref_donor) && setequal(x$Acceptors, ref_acceptor) && abs(start - x$TSS) < 50 && abs(end - x$TTS) < 50){
        return(TRUE)
      }
    }
    return(FALSE)
  })
  TP <- sum(matches==T)
  FP <- sum(matches==F)
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
del.file <- args[3] # deleted transcripts file preparatory step
cat.file <- args[4] # categories files
expr.file <- args[5] # sim counts file sim step


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

data.query <- inner_join(data.class(), data.junction)
data.query <- data.query[,c('isoform', 'structural_category', 'Donors', 'Acceptors', 'TSS_genomic_coord', 'TTS_genomic_coord')]

# Read deleted file
data.del <- read.table(del.file, header=T, as.is=T, sep="\t")
data.del$Donors <- lapply(data.del$Donors, function(x){
  as.numeric(unlist(strsplit(as.character(x), ",")))+1
})
data.del$Acceptors <- lapply(data.del$Acceptors, function(x){
  as.numeric(unlist(strsplit(as.character(x), ",")))+1
})
data.del$structural_category = factor(data.del$SC,
                                      labels = xaxislabelsF1,
                                      levels = xaxislevelsF1,
                                      ordered=TRUE)
data.del$SC <- NULL

# Read cat file
data.cat <- read.table(cat.file, header=T, as.is=T, sep="\t")
# Sum 1 to donors and deduct 1 to acceptors to match with SQANTI3 format
data.cat$Donors <- lapply(data.cat$Donors, function(x){
  as.numeric(unlist(strsplit(as.character(x), ",")))+1
})
data.cat$Acceptors <- lapply(data.cat$Acceptors, function(x){
  as.numeric(unlist(strsplit(as.character(x), ",")))+1
})
data.cat$structural_category = factor(data.cat$SC,
                                        labels = xaxislabelsF1,
                                        levels = xaxislevelsF1,
                                        ordered=TRUE)
data.cat$SC <- NULL

# Read sim expr file
data.expr <- read.table(expr.file, header=F, as.is = T, sep="\t")
colnames(data.expr) <- c('TransID', 'counts')

# Transcript simulated with reference
data.known <- data.cat[which(data.cat$TransID %in% data.expr$TransID & !(data.cat$TransID %in% data.del$TransID)),]

known.metrics <- list()
for (sc in xaxislabelsF1){
  list[sc] <- calc_known_stats(sc, data.query, data.known)
}

#######################################
#                                     #
#     TABLE AND PLOT GENERATION       #
#                                     #
#######################################

# -------------------- 
# -------------------- 
# TABLE INDEX
# t1:

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
# PLOT 1: simulated expression profile
expr.dist <- data.expr
expr.dist$type <- lapply(expr.dist$TransID, function(x){
  ifelse(x %in% del.file$TransID, 'novel', 'known')
})

p1 <- expr.dist %>%
  ggplot( aes(x=counts, fill=type)) +
  geom_histogram(alpha=0.6, position = 'identity') +
  mytheme +
  scale_fill_manual(values = c('orange', 'darkcyan'), guide='none') +
  ggtitle('Simulated counts distribution')

# PLOT 2: structural classification
p2 <- data.class %>%
  ggplot(aes(x=structural_category)) +
  geom_bar(aes(y = (..count..)/sum(..count..)*100, alpha=0.6, fill=structural_category), color="black", size=0.3, width=0.7) +
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













