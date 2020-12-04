#######################################
### Generic script to make boxplots ###
#######################################

### load libraries
library(ggplot2)
library(RColorBrewer)
library(tidyr)
library(ggpubr)
library(dplyr)

#############################
### plot a single boxplot ###
#############################
mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/mmrf_select-genes_tmm_log2.rds")
group <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/factors/translocations_factors.rds")
## choose genes outside of 24 targets
mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/mmrf_tmm_log2.rds")
g <- c("FCRL1", "FCRL2", "FCRL3", "FCRL4", "FCRL5", "FCRL6")
g <- c("MDM2", "MDM4", "TP53", "CDKN1A")
g <- c("FUT3", "FGFR3", "WHSC1")
g <- c("ICAM1", "MAF", "NFKB1", "FUT3", "WHSC1", "FGFR3")
g <- c("ZFP36L1", "ANGEL1")
g <- c("FUT3", "FGFR3", "WHSC1", "BCMA", "CD138", "SLAMF7", "CD56", "CD200")
g <- read.csv("~/Documents/osu/Perna_lab/surface_proteins/ms_mmrf_overlap_155.csv", stringsAsFactors=TRUE)
g <- g$UniProt.gene
mat <- mat[which(rownames(mat) %in% c("MAF", "ICAM1", "NFKB1")), ] ## ICAM1 with MAF
mat <- mat[which(rownames(mat) %in% c("FUT3", "FGFR3")), ] ## FUT3 with FGFR3
which(rownames(mat) %in% g)
mat <- mat[which(rownames(mat) %in% g), ]
rowMeans(mat)
cor(t(mat), method="spearman")

## re-level as needed
levels(group)
# group <- factor(group, levels=levels(group)[c(2, 1)])

# ## Rename BI_TP53 factors:
# levels(group)[levels(group)==2] <- "Normal"
# levels(group)[levels(group)==1] <- "Impacted"
# levels(group)[levels(group)==0] <- "Inactive"

##############################
######## module genes ########
group <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/factors/translocations_factors.rds")
levels(group)
### set up to get module genes:
## read in modules
modules <- readRDS("~/Documents/osu/Perna_lab/network_mining/mmrf/rds/gamma60_log2-mat_spearman.rds")
## full mmrf expression matrix
mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/mmrf_tmm_log2.rds")
## make matrix
mat <- mat[which(rownames(mat) %in% c(modules@clusters.names[[23]], "ICAM1", "MAF", "CBFA2T3", "IRF8", "AP1G2", "ARG2")), ]
mat <- mat[which(rownames(mat) %in% c(modules@clusters.names[[47]], "FUT3", "ZFP36L1")), ]
mat <- mat[which(rownames(mat) %in% c(modules@clusters.names[[9]], "CD19")), ]
mat <- mat[which(rownames(mat) %in% c(modules@clusters.names[[20]], "BCMA")), ]
mat <- mat[which(rownames(mat) %in% c(modules@clusters.names[[11]], "IGF1R", "CD70", "MYB")), ]
mat <- mat[which(rownames(mat) %in% c(modules@clusters.names[[44]], "LY9")), ]
mat <- mat[which(rownames(mat) %in% c(modules@clusters.names[[27]], "BST2", "IRF9")), ]
##############################

# ##################### only for BCMA Arkansas data #####################
# ###### detour: get Arkansas baseline or relapse group ######
# mat <- mat[, which(colnames(mat) %in% names(group)[group=="baseline"])]
# group <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/factors/arkansas_baseline_bcma_factors.rds")
# levels(group)
# ### match up sample names
# all.equal(colnames(mat), names(group)) # TRUE
# #######################################################################

### specify names for plotting:
filepath <- "~/Downloads/Rshiny_cell-markers/Rplots/additional_analysis/t4-14_upreg_genes_translocation-pts"
legend.name <- ""
plot.title <- ""
fact <- TRUE # Plotting by factors? TRUE or FALSE
gene <- "BCMA"
gene <- c("MDM2", "MDM4", "CDKN1A", "TP53", "BCMA", "GUCY2F")
upreg.genes <- c("FGFR3", "MMSET") # make sure to change rownames of WHSC1 to MMSET
## change gene names
rownames(mat)[rownames(mat)=="WHSC1"] <- "MMSET"

############# automated process #############
### match up sample names
all.equal(colnames(mat), names(group)) # FALSE
## if FALSE and dimensions of mat is larger than group:
if(!(all.equal(colnames(mat), names(group)) == TRUE)[1]) {
  mat <- mat[, which(colnames(mat) %in% names(group))]
  ## check again
  all.equal(colnames(mat), names(group)) # TRUE
}

# ### calculate correlations in specific patient subsets
# mat.subset <- mat[, which(colnames(mat) %in% names(group)[group=="t(11;14)"])]
# cor(t(mat.subset))
cor(t(mat))
write.csv(cor(t(mat), method="spearman"), file="~/Downloads/Rshiny_cell-markers/Rplots/additional_analysis/t4-14_upreg_genes_translocation-pts.csv")

### wrangle data
d <- data.frame(Genes=rownames(mat), mat)
# head(d[, 1:5])
dd <- d %>% gather(Factor, Expression, -Genes)
# head(dd)

### output gene order for surface protein selection
d.rowMeans <- rowMeans(mat)
head(d.rowMeans)
d.rowMeans <- d.rowMeans[order(d.rowMeans, decreasing=TRUE)]
head(d.rowMeans)
write.csv(d.rowMeans, "~/Documents/osu/Perna_lab/surface_proteins/ms_mmrf_overlap_115_genes_RNAseq_ordered.csv")

## plot heatmap?
library(reshape2)
dd <- melt(d.rowMeans)
dd$Gene <- rownames(dd)
head(dd)
ncols <- 10
mycolors <- colorRampPalette(brewer.pal(7, "RdBu"))(ncols)
p <- ggplot(dd, aes(x=Gene, y=NULL, fill=value)) +
  geom_tile() +
  # theme_pubr() + # doesn't look good when used with coord_fixed()
  coord_fixed() +
  scale_fill_gradientn(colours = rev(mycolors), name="log2(expr)") +
  theme(axis.text.x = element_text(angle=90),
        text=element_text(size=10)) +
  xlab("") + 
  ylab("")
p

######### If selecting genes: #########
if(!is.null(gene)) {
  dd <- dd[which(dd$Genes %in% gene), ]
}

### re-format for faceted boxplot
dd$Upreg <- factor(ifelse(dd$Genes %in% upreg.genes, 1, 0))

######### If there is a factor: #########
if(fact) {
  ## Factor right now are sample names, so change to actual factor names
  for(fac in levels(group)) {
    dd$Factor[which(dd$Factor %in% names(group)[group==fac])] <- fac
  }
  dd$Factor <- factor(dd$Factor)
  levels(dd$Factor)
  ## Re-level the factor according to order given in group
  dd$Factor <- factor(dd$Factor, levels=levels(group))
}
levels(dd$Factor)

# ### Save plot
tiff(filename=paste0(filepath, ".tiff"), res = 300, width = 10, height = 5, units = "in")
pdf(file=paste0(filepath, ".pdf"), width=6.5, height=4.5)
ncols <- ifelse(is.factor(dd$Factor), length(levels(dd$Factor)), 1) # returns num factors or 1
mycolors <- colorRampPalette(brewer.pal(3, "Greys"))(ncols)
# ## provide gene order from make_boxplots
# genes.order <- readRDS("~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_order.rds")
# dd$Genes <- factor(dd$Genes, levels=genes.order)
dd$Genes <- reorder(factor(dd$Genes), -dd$Expression) # re-order by expression
### plot - CHANGE fill as needed
ggplot(dd, aes(x=Genes, y=Expression, fill=Factor)) +
  geom_boxplot(notch=FALSE, outlier.shape=NA, show.legend=TRUE) +
  scale_fill_manual(values=mycolors, name=legend.name) +
  # ylim(-10, 20) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, size=10, hjust=0.95,vjust=0.5)) +
  theme(axis.text.x=element_text(face="italic")) +
  # theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  # theme(axis.title = element_text(size=12)) +
  theme(legend.position="right") +
  ggtitle(plot.title) +
  ylab("log2(Expression)") +
  xlab("") +
  facet_grid(~dd$Upreg, scales="free", space="free")
dev.off()

###################################################
################ plot using ggpubr ################
genes.order <- readRDS("~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_order.rds")
dd$Genes.relevel <- factor(dd$Genes, levels=genes.order)
head(dd)
## make data frame with statistical values for manual option
# results <- data.frame(gene=c("CD70", "IGF1R"),
#                       group1=c("1", "1"),
#                       group2=c("3", "3"),
#                       p.signif=c("**", "*"),
#                       y.position=c(16, 16))
results <- data.frame(read.csv(file="~/Downloads/Rshiny_cell-markers/Rtables/updated_stats/Translocations/Translocations_upreg-genes.csv"))
## choose colors
ncols <- ifelse(is.factor(dd$Factor), length(levels(dd$Factor)), 1) # returns num factors or 1
mycolors <- colorRampPalette(brewer.pal(5, "Blues"))(ncols)
## plot with options
tiff(filename=paste0(filepath, ".tiff"), res = 300, width = 10, height = 6, units = "in")
p <- ggboxplot(dd, x="Genes.relevel", y="Expression", fill="Factor",
          palette=mycolors, outlier.shape=NA, notch=TRUE) + 
  stat_pvalue_manual(data=results, x="gene", label="p.signif", 
                     position=position_identity())
# label="p.signif"
# label="p = {p}"
ggpar(p, main=plot.title, xlab="", ylab="log2(Expression)",
      legend="top", legend.title=legend.name)
# legend.title=legend.name
dev.off()

#############################
### plot combined boxplot ###
#############################

### load data
# mmrf_select-genes_tmm_log2_37.rds
mat1 <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/arkansas_select-probes_expr_log2_37.rds")
mat2 <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/boston_select-probes_expr_log2_37.rds")
mat3 <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/mmrf_select-genes_tmm_log2_37.rds")

keep.genes <- read.csv("~/Downloads/Rshiny_cell-markers/Rdata/genes_keep.csv", stringsAsFactors=FALSE)

mat1 <- mat1[match(keep.genes$Genes, rownames(mat1)), ]
mat2 <- mat2[match(keep.genes$Genes, rownames(mat2)), ]
mat3 <- mat3[match(keep.genes$Genes, rownames(mat3)), ]

# saveRDS(mat1, file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/genes24_arkansas_v2.rds")
# saveRDS(mat2, file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/genes24_boston_v2.rds")
# saveRDS(mat3, file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/genes24_mmrf_v2.rds")

# mat1 <- mat1[which(rownames(mat1) %in% keep.genes$Genes), ]
# mat2 <- mat2[which(rownames(mat2) %in% keep.genes$Genes), ]
# mat3 <- mat3[which(rownames(mat3) %in% keep.genes$Genes), ]

### combine all three datasets into one boxplot
### calculate z-score for each gene
### given dataframe d, calculate z-score in a col-wise manner

## dataset 1
d <- data.frame(Genes=rownames(mat1), mat1)
d1.zscore <- sapply(d[-c(1)], function(x) scale(x))
d1.zscore[1:5, 1:5]
d1.zscore <- data.frame(Genes=d$Genes, d1.zscore)

## dataset 2
d <- data.frame(Genes=rownames(mat2), mat2)
d2.zscore <- sapply(d[-c(1)], function(x) scale(x))
d2.zscore[1:5, 1:4]
d2.zscore <- data.frame(Genes=d$Genes, d2.zscore)

## dataset 3
d <- data.frame(Genes=rownames(mat3), mat3)
d3.zscore <- sapply(d[-c(1)], function(x) scale(x))
d3.zscore[1:5, 1:5]
d3.zscore <- data.frame(Genes=d$Genes, d3.zscore)

######### order genes by MMRF expression ######### 
# mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/genes24_mmrf.rds")
# ### wrangle data
# d <- data.frame(Genes=rownames(mat), mat)
# head(d[, 1:5])
# dd <- d %>% gather(Factor, Expression, -Genes)
# head(dd)
# genes.order <- reorder(dd$Genes, -dd$Expression)
# genes.order <- factor(levels(genes.order), levels=levels(genes.order))
# saveRDS(genes.order, file="~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_order.rds")

# keep.genes <- read.csv("~/Downloads/Rshiny_cell-markers/Rdata/genes_keep.csv", stringsAsFactors=FALSE)
# 
# d1.zscore <- d1.zscore[which(d1.zscore$Genes %in% keep.genes$Genes), ]
# d2.zscore <- d2.zscore[which(d2.zscore$Genes %in% keep.genes$Genes), ]
# d3.zscore <- d3.zscore[which(d3.zscore$Genes %in% keep.genes$Genes), ]

### consolidate datasets into one
all.equal(d1.zscore$Genes, d2.zscore$Genes)
all.equal(d1.zscore$Genes, d3.zscore$Genes)
combined.zscore <- cbind(d1.zscore, d2.zscore[-c(1)], d3.zscore[-c(1)])
dim(combined.zscore) # genes in row
combined.zscore[1:5, 1:5]

ttg <- gather(combined.zscore, 
              key="Dataset",
              value="Expression",
              -Genes)
head(ttg)

### differentiate sample names by dataset
ttg$Dataset[which(ttg$Dataset %in% colnames(mat1))] <- "GSE31161"
ttg$Dataset[which(ttg$Dataset %in% colnames(mat2))] <- "GSE9782"
ttg$Dataset[which(ttg$Dataset %in% colnames(mat3))] <- "CoMMpass"
table(ttg$Dataset)
ttg$Dataset <- factor(ttg$Dataset)
levels(ttg$Dataset)
ttg$Dataset <- factor(ttg$Dataset, levels(ttg$Dataset)[c(2, 3, 1)])
saveRDS(ttg, file="~/Downloads/Rshiny_cell-markers/Rdata/plotting/combined_24_genes_v2.rds")
control <- filter(ttg, Genes=="CLOCK")
mean(control$Expression)

### Save plot
tiff(file=paste0(filepath, ".tiff"), res = 300, width = 10, height = 7, units = "in")
ncols <- 3
mycolors <- colorRampPalette(brewer.pal(5, "Greys"))(ncols)
## provide gene order from make_boxplots
# genes.order <- readRDS("~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_order.rds")
# ttg$Genes.relevel <- factor(ttg$Genes, levels=genes.order)
# order for 25 genes:
ttg$Genes.relevel <- reorder(ttg$Genes, -ttg$Expression)
saveRDS(levels(ttg$Genes.relevel), "~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_order.rds")
## plot - CHANGE fill as needed
ggplot(ttg, aes(x=Genes.relevel, y=Expression, fill=Dataset)) +
  geom_boxplot(notch=TRUE, outlier.shape=NA, show.legend=TRUE) +
  scale_fill_manual(values=mycolors) +
  theme_pubr() +
  theme(axis.text.x = element_text(angle=90, size=10)) +
  theme(axis.title = element_text(size=12)) +
  ggtitle("Gene expression of cell surface markers") + 
  # ylab("log2(Expression)") +
  ylab("Relative expression (z-score)") +
  xlab("Genes")
dev.off()


#### ggpubr 
p <- ggboxplot(dd, x="Genes.relevel", y="Expression", fill="Factor",
               palette=mycolors, outlier.shape=NA, notch=TRUE) + 
  stat_pvalue_manual(data=results, x="gene", label="p.signif", 
                     position=position_identity())
# label="p.signif"
# label="p = {p}"
ggpar(p, main=plot.title, xlab="", ylab="log2(Expression)",
      legend="top", legend.title=legend.name)
