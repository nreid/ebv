library(tidyverse)
library(DESeq2)
library(pheatmap)
library(ggrepel)


#################################
# check basic raw sequence stats
#################################

# read in read summary data from multiqc before/after trimming---------------------------------------------------------
mqcb <- read.table("../../results/04a_multiqc_small/multiqc_data/multiqc_fastqc.txt",sep="\t",header=TRUE)
	mqcb <- mqcb[,c(1,5,8,9,10)]
	colnames(mqcb) <- c("sample","total_raw","GC_raw","total_deduplicated_raw","meanlength_raw")
mqca <- read.table("../../results/07a_multiqc_small_trimmed/multiqc_data/multiqc_fastqc.txt",sep="\t",header=TRUE)
	mqca <- mqca[,c(1,5,8,9,10)]
	colnames(mqca) <- c("sample","total_trim","GC_trim","total_deduplicated_trim","meanlength_trim")
	mqca[,1] <- gsub(".trim","",mqca[,1])

meta <- full_join(mqcb,mqca) %>% 
	separate(sample,"_",into=c("line","kd","id"),remove=FALSE) %>%
	mutate(kd=gsub("^a","A",kd) %>% gsub("^d","D",.) %>% gsub("^l","L",.))
rownames(meta) <- meta[,1]

# make some plots---------------------------------------------------------------------------------------------------------

# characteristics of read data 

ggplot(meta, aes(x=line, y=total_trim)) + geom_boxplot() + geom_jitter(width=0.1)

# before/after trimming

ggplot(meta,aes(x=total_raw,y=total_trim,color=line)) + geom_point()

ggplot(meta,aes(x=GC_raw,y=GC_trim,color=line)) + geom_point() + geom_jitter()

ggplot(meta,aes(x=total_deduplicated_raw,y=total_deduplicated_trim,color=line)) + geom_point()


###############################################
# analyze results from ShortStack
###############################################

# read in count data from ShortStack-------------------------------------------------------------------------------

sscounts <- read.table("../../results/small_analysis/shortstack_results/Counts.txt", header=TRUE)
	colnames(sscounts) <- gsub(".trim","",colnames(sscounts))

# read in Results.txt

ssres <- read.table("../../results/small_analysis/shortstack_results/Results.txt",comment.char="",header=TRUE)

# exploratory plots----------------------------------------------------------------------------------------------------------------

# what does the correlation among samples look like?

cor(log(ssmat+1)) %>% image()


# put count data into DESeq2 object------------------------------------------------------------------------------------

# count matrix
ssmat <- as.matrix(sscounts[,-c(1:3)])
	rownames(ssmat) <- sscounts$Locus

# sample data table
coldata <- meta[,2:4]
	rownames(coldata) <- meta[,1]
	# match order to ssmat
	coldata <- coldata[colnames(ssmat),]
	# test to make sure it's right:
	all(rownames(coldata) == colnames(ssmat))

# make deseq object
	# design matrix is cell line + knockdown. 
	# unclear how to set this up. certainly can't do interactions terms with zero replication within cell lines?
	# maybe this whole thing is only suitable for preliminary analysis?
dds <- DESeqDataSetFromMatrix(countData = ssmat,
                              colData = coldata,
                              design = ~ line + kd)

# relevel the knockdown factors so CONTROL is the base level
	# not going to relevel cell lines. 
facs <- c("CONTROL", "Ago1KD", "Ago2KD", "Ago3KD", "Ago4KD", "DicerKD", "DroshaKD", "LaKD")

dds$kd <- factor(dds$kd,levels=facs)


#############################
# do analysis
#############################

# fit model---------------------------------------------------------------------
dds <- DESeq(dds)

# some visualizations----------------------------------------------------------------------

figdir <- "../../doc/figs/"

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- vst(dds, blind=FALSE)

plotPCA(vsd, intgroup="line")

# just using ebv
ebv <- grep("EBV",rownames(vsd))
plotPCA(vsd[ebv,], intgroup="line")

# alternatively, using ggplot

dat <- plotPCA(vsd, intgroup="line",returnData=TRUE)
percentVar <- round(100 * attr(dat, "percentVar"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group)) +
  geom_point() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel(aes(label=name))
p

ggsave(filename="shortstack_PCA.png",plot=p,device="png",path=figdir)


# heatmap


png(filename=paste0(figdir,"shortstack_heatmap_ebv_only.png"), width=800,height=800)

df <- coldata[,-3]

pheatmap(
  assay(vsd)[ebv,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df
  )

dev.off()

###############################################
# analyze results from miRDeep2
###############################################

# load in results files and fix them up-----------------------------------------------------------------------------------------------------------------------------

# get the config file, which maps miRDeep2 sample numbers to our sample IDs
cf <- read.table("../../results/small_analysis/mirdeep2_all/config.txt") %>%
	mutate(V1=str_extract(V1,regex("[^/\\.]+(?=\\.)")))

# get the table of basic data on mirBase known miRNAs that were found in our data
	# ignoring novel ones because there weren't many good ones, and none on EBV
mdr <- read.table("../../results/small_analysis/mirdeep2_all/mature_miRNA_found.txt",skip=1,header=TRUE,sep="\t",comment.char="")

# get the table of read counts for known miRNAs, rename the header lines according to config.txt
	# some mature miRNAs are identical and have identical counts, but come from different precursors
mdc <- read.table("../../results/small_analysis/mirdeep2_all/quant/miRNAs_expressed_all_samples_1616442575.csv",comment.char="",header=TRUE,sep="\t")
	colnames(mdc)[5:36] <- cf[,1]
	colnames(mdc)[37:68] <- paste(cf[,1],"norm",sep="_")

# filter it down to detected miRNAs
	# note that one EBV miRNA has non-zero read counts but was not detected by miRDeep2
	# overall only 1% of counted reads are discarded by this filter
mdc <- filter(mdc, X.miRNA %in% mdr$mature.miRBase.miRNA)

# extract only the raw counts, as a matrix
mdraw <- as.matrix(mdc[,5:36])
	rownames(mdraw) <- mdc[,1]


# put count data into DESeq2 object------------------------------------------------------------------------------------


# sample data table
coldata <- meta[,2:4]
	rownames(coldata) <- meta[,1]
	# match order to mdraw
	coldata <- coldata[colnames(mdraw),]
	# test to make sure it's right:
	all(rownames(coldata) == colnames(mdraw))

# make deseq object
	# design matrix is cell line + knockdown. 
	# unclear how to set this up. certainly can't do interactions terms with zero replication within cell lines?
	# maybe this whole thing is only suitable for preliminary analysis?
	# NOTE: some mature miRNAs are identical and have identical counts, but come from different precursors

mdd <- DESeqDataSetFromMatrix(countData = mdraw,
                              colData = coldata,
                              design = ~ line + kd)

# relevel the knockdown factors so CONTROL is the base level
	# not going to relevel cell lines. 
facs <- c("CONTROL", "Ago1KD", "Ago2KD", "Ago3KD", "Ago4KD", "DicerKD", "DroshaKD", "LaKD")

mdd$kd <- factor(mdd$kd,levels=facs)


# fit model---------------------------------------------------------------------
mdd <- DESeq(mdd)

# some visualizations----------------------------------------------------------------------

figdir <- "../../doc/figs/"

# PCA plot

# normalized, variance-stabilized transformed counts for visualization
vsd <- varianceStabilizingTransformation(mdd)

plotPCA(vsd, intgroup="line")

# just using ebv
ebv <- grep("ebv",rownames(vsd))
plotPCA(vsd[ebv,], intgroup="line")

# alternatively, using ggplot

dat <- plotPCA(vsd, intgroup="line",returnData=TRUE)
percentVar <- round(100 * attr(dat, "percentVar"))

p <- ggplot(dat,aes(x=PC1,y=PC2,col=group)) +
  geom_point() + 
  xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  ylab(paste0("PC2: ",percentVar[2],"% variance")) +
  geom_text_repel(aes(label=name))
p

ggsave(filename="mirdeep2_PCA.png",plot=p,device="png",path=figdir)

# heatmap

png(filename=paste0(figdir,"mirdeep2_heatmap_ebv_only.png"), width=800,height=800)

df <- coldata[,-3]

pheatmap(
  assay(vsd)[ebv,], 
  cluster_rows=TRUE, 
  show_rownames=TRUE,
  cluster_cols=TRUE,
  annotation_col=df
  )

dev.off()


#######################
# table output
#######################


reads <- data.frame(total=meta[colnames(ssmat),9],shortstack=colSums(ssmat),mirdeep2=colSums(mdraw)) %>%
	mutate(shortstack_pct=shortstack/total,mirdeep2_pct=mirdeep2/total)

ebvreads <- data.frame(total=meta[colnames(ssmat),9],shortstack=colSums(ssmat[grep("EBV",rownames(ssmat)),]),mirdeep2=colSums(mdraw[grep("ebv",rownames(mdraw)),])) %>%
	mutate(shortstack_pct=shortstack/total,mirdeep2_pct=mirdeep2/total)

write.csv(reads,"../../doc/reads.csv")
write.csv(ebvreads,"../../doc/ebvreads.csv")