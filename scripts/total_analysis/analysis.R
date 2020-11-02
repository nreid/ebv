#################################################
# Load necessary libraries
#################################################

library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggrepel)

#################################################
# Set some directories for input/output files
#################################################

list.files()

# Script assumes kallisto counts are located in this directory
count_dir <- "../../results/total_analysis/counts/" # or appropriate path to the counts 


#################################################
# Read in count data
#################################################



# create a empty dataframe called co to merge the data into
co <- data.frame()

# using for loop read all the count files in the count_dir path
# iterate over each "abundance.tsv" file
# merge each file, i, into a single data frame

count_files <- list.files(path = count_dir,recursive=TRUE,pattern="abundance.tsv",full.names=TRUE)

for (i in count_files) {

  # print the file that is being loaded
  print(paste0("reading file: ", i))

  # read file i as a data frame
  f <- read.table(i, sep = "\t", header = TRUE)[,c(1,4)]

  # rename the columns
  sam <- str_extract(i, regex("[^\\/]+(?=/abundance.tsv)"))
  colnames(f) <- c("gene_id", sam)

  #if the counts is empty just copy the f to m
  if(length(co) == 0){
    co <- f
    
  } else 
  {
    #if the dataframe is not empty then merge the data
    co <- merge(co, f, by.x = "gene_id", by.y = "gene_id")
  }
  rm(f)
}

#grab the rows from the 1st column and use it as the row-names in the dataframe
rownames(co) <- co[,1]

# table containing only counts
m <- co[,-1]

# create a vector containing transcript lengths
# read in length
df_length <- read.table(count_files[1], sep = "\t", header = TRUE)[,1:2]

# create named vector
mylength <- df_length$length
names(mylength) <- df_length[,1]


#################################################
# Create metadata table
#################################################

# here we create the table
Sample = str_extract(count_files, regex("[^\\/]+(?=/abundance.tsv)"))
CellLine = substr(Sample,1,2)
myfactors <- data.frame(Sample, CellLine)
myfactors

# the order and identity of samples in data frame m and myfactor must be the same
# ensure all the columns in dataframe m and myfactor is present
all(colnames(m) %in% myfactors$Sample)

# Then check the order is the same in both objects
all(colnames(m) == myfactors$Sample)

# If are not, so order them according to the sample names
m <- m[, myfactors$Sample]

# Now check the order is correct after sorting 
all(colnames(m) == myfactors$Sample)


#################################################
# Preliminary data exploration and filtering
#################################################

# how many transcripts did we quantify expression for?
dim(m)

# how many (non-contaminant) transcripts have expression for n samples?
prevalence <- rowSums(m > 0)
table(prevalence)

# what proportion of expression data maps to low prevalence transcripts?
colSums(m[prevalence<2,])/colSums(m)

# what is the distribution of total expression across samples?
rowSums(m) %>% log(.,10) %>% hist(.,100)

subg <- rowSums(m) > 20
table(subg)

m <- m[subg,]
mylength <- mylength[subg]

# make an MDS 


# downstream analysis shows that miA_1_S5 is a very problematic sample. exclude it.
	# in the PCA, it is way off by itself on PC1, which explains 43% of the variance of the rld transformed counts

# m <- m[,!(colnames(m) %in% "miA_1_S5")]
# myfactors <- myfactors[myfactors[,1]!="miA_1_S5",]

#################################################
#  Analyze the data using DESeq2
#################################################

# we've already created the necessary component objects
  # use them to create a "DESeqDataSet" object
  # we are rounding 'm' because DESeq2 expects integer counts,
  # but Kallisto estimates the counts, resulting in fractional numbers

ddsKAL <- DESeqDataSetFromMatrix(
  countData = round(m),
  colData = myfactors,
  design = ~ CellLine
  )

######################################################
# Reset treatment factors
######################################################

# To see the levels as they are now:
ddsKAL$CellLine

# To replace the order with one of your choosing, create a vector with the order you want:
treatments <- c("RP","mi","sh","RE")

# Then reset the factor levels:
ddsKAL$CellLine <- factor(ddsKAL$CellLine, levels = treatments)

# verify the order
ddsKAL$CellLine

######################################################
# Run the statistical analysis
######################################################

dds <- DESeq(ddsKAL)

######################################################
# Data visualization
######################################################

# set a directory to save figures in
figdir <- "../../doc/figs/"
dir.create(figdir)

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

dat <- plotPCA(vsd, intgroup="CellLine",returnData=TRUE)
percentVar <- round(100 * attr(dat, "percentVar"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group)) +
	geom_point() + 
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	geom_text_repel(aes(label=name))
p

# to save the plot
ggsave(filename="1_PCAv1.png",plot=p,device="png",path=figdir)


######################################################
# RE-Run the statistical analysis
######################################################

# exclude sample miA_1_S5 as an outlier.
	# it is in index position 1

dds <- DESeq(ddsKAL[,-1])


######################################################
# Back to data visualization
######################################################

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

dat <- plotPCA(vsd, intgroup="CellLine",returnData=TRUE)
percentVar <- round(100 * attr(dat, "percentVar"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group)) +
	geom_point() + 
	xlab(paste0("PC1: ",percentVar[1],"% variance")) +
	ylab(paste0("PC2: ",percentVar[2],"% variance")) +
	geom_text_repel(aes(label=name))
p

# to save the plot
ggsave(filename="2_PCAv2.png",plot=p,device="png",path=figdir)


######################################################
# Get a table of results
######################################################

# list coefficient names

resultsNames(dds)

# get results tables for all 6 pairwise contrasts

miRP <- results(dds, name="CellLine_mi_vs_RP")
shRP <- results(dds, name="CellLine_sh_vs_RP")
RERP <- results(dds, name="CellLine_RE_vs_RP")
REmi <- results(dds, contrast=c("CellLine","RE","mi"))
shmi <- results(dds, contrast=c("CellLine","sh","mi"))
REsh <- results(dds, contrast=c("CellLine","sh","RE"))

# get a quick summary of the table
summary(shRP)


######################################################
# Get a table of shrunken log2 fold changes
######################################################

# see coefficient names:
resultsNames(dds)

# get shrunken log fold changes, specifying the coefficient or contrast

miRPs <- lfcShrink(dds, coef="CellLine_mi_vs_RP",type="ashr")
shRPs <- lfcShrink(dds, coef="CellLine_sh_vs_RP",type="ashr")
RERPs <- lfcShrink(dds, coef="CellLine_RE_vs_RP",type="ashr")
REmis <- lfcShrink(dds, contrast=c("CellLine","RE","mi"),type="ashr")
shmis <- lfcShrink(dds, contrast=c("CellLine","sh","mi"),type="ashr")
REshs <- lfcShrink(dds, contrast=c("CellLine","sh","RE"),type="ashr")
	
######################################################
# Back to data visualization
######################################################

# MA plot
plotMA(RERPs, ylim=c(-4,4))

plotMA(shRPs, ylim=c(-4,4))

# distribution of log2 fold changes:
  # there should be a peak at 0
  # some weird spikes are accounted for by low expression genes filtered by independent filtering
hist(RERP$log2FoldChange,breaks=2000,xlim=c(-4,4))
abline(v=0,col="red",lwd=2)

##############

#Volcano plot

# negative log-scaled adjusted p-values
log_padj <- -log(res_shrink$padj,10)
log_padj[log_padj > 100] <- 100

# plot
plot(x=res_shrink$log2FoldChange,
     y=log_padj,
     pch=20,
     cex=.2,
     col=(log_padj > 3)+1, # color padj < 0.1 red
     ylab="negative log-scaled adjusted p-value",
     xlab="shrunken log2 fold changes")

#############



##############

# heatmap of DE genes

# regularized log transformation of counts
rld <- rlog(dds, blind=FALSE)

# first for RERP (RPE vs RE)
# get top 50 log fold change genes (excluding cook's cutoff outliers)
top50 <- data.frame(RERPs) %>%
  filter(!is.na(padj)) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() %>% 
  head(.,n=50)

df <- data.frame(colData(dds)[,"CellLine"])
  rownames(df) <- colnames(dds)
  colnames(df) <- "CellLine"

pheatmap(
  assay(rld)[top50,c(6,7,8,1,2,9,10,11,3,4,5)], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  annotation_col=df
  )

# to save the plot
png(filename=paste0(figdir,"heatmap_REvsRPE.png"), width=800,height=800)
pheatmap(
  assay(rld)[top50,c(6,7,8,1,2,9,10,11,3,4,5)], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  annotation_col=df
  )
dev.off()

# now for sh vs RPE
# get top 50 log fold change genes (excluding cook's cutoff outliers)
top50 <- data.frame(shRPs) %>%
  filter(!is.na(padj)) %>% 
  arrange(-abs(log2FoldChange)) %>% 
  rownames() %>% 
  head(.,n=50)

df <- data.frame(colData(dds)[,"CellLine"])
  rownames(df) <- colnames(dds)
  colnames(df) <- "CellLine"


pheatmap(
  assay(rld)[top50,c(6,7,8,1,2,9,10,11,3,4,5)], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  annotation_col=df
  )

# to save the plot
png(filename=paste0(figdir,"heatmap_shvsRPE.png"), width=800,height=800)
pheatmap(
  assay(rld)[top50,c(6,7,8,1,2,9,10,11,3,4,5)], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=FALSE,
  annotation_col=df
  )
dev.off()

