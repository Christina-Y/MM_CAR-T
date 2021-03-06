### genes20 dataframes are generated by gse_data_processing.Rmd and mmrf_cell_markers.Rmd
genes20.boston <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_boston.rds")
genes20.arkansas <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_arkansas.rds")
genes20.mmrf <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes25_mmrf.rds")

### plotting boxplot
library(ggplot2)
library(RColorBrewer)
library(tidyr)

### plot single figure with boxplots of all genes
d <- data.frame(Genes=rownames(genes20.mmrf), 
                genes20.mmrf)
head(d[, 1:5])
dd <- d %>% gather(patient, expression, -Genes)
head(dd)

# Define the number of colors you want
ncols <- 29
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(ncols)
ggplot(dd, aes(x=reorder(Genes, -expression), y=expression, fill=Genes)) + 
  geom_boxplot(notch=TRUE, outlier.size=1, show.legend=FALSE) +
  scale_fill_manual(values=mycolors) +
  theme(axis.text.x = element_text(angle=90, size=10)) +
  theme(axis.title = element_text(size=12)) + 
  ggtitle("Gene Expression in MMRF Data") + ylab("log2(Expression)")


### combine all three datasets into one boxplot
## calculate z-score for each gene
## given dataframe d, calculate z-score in a col-wise manner
## boston
d <- data.frame(Genes=rownames(genes20.boston), 
                genes20.boston)
# write.csv(d, file="~/Downloads/Rshiny_cell-markers/boston_to_zscore.csv", row.names=FALSE)
boston.zscore <- sapply(d[-c(1)], function(x) scale(x))
boston.zscore[1:5, 1:5]
boston.zscore <- data.frame(Genes=d$Genes, boston.zscore)

## arkansas
d <- data.frame(Genes=rownames(genes20.arkansas), 
                genes20.arkansas)
arkansas.zscore <- sapply(d[-c(1)], function(x) scale(x))
arkansas.zscore <- data.frame(Genes=d$Genes, arkansas.zscore)
arkansas.zscore[1:5, 1:5]

## mmrf
d <- data.frame(Genes=rownames(genes20.mmrf), 
                genes20.mmrf)
mmrf.zscore <- sapply(d[-c(1)], function(x) scale(x))
mmrf.zscore <- data.frame(Genes=d$Genes, mmrf.zscore)
mmrf.zscore[1:5, 1:5]

### plot single dataset of z-score to check
## reshape d in order for plotting
dd <- boston.zscore %>% gather(patient, expression, -Genes)
head(dd)

# define colors for 24 genes
ncols <- 24
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(ncols)
ggplot(dd, aes(x=Genes, y=expression, fill=Genes)) +
  geom_boxplot(notch=TRUE, outlier.size=1, show.legend=FALSE) +
  scale_fill_manual(values=mycolors) +
  theme(axis.text.x = element_text(angle=90, size=10)) +
  theme(axis.title = element_text(size=12)) +
  ggtitle("Gene Expression in Boston Data") + 
  ylab("Relative expression (z-score)")

### how to consolidate all three datasets into one?
## Boston patient names all start with GSM2
## Arkansas patient names start with GSM7 or GSM9
## MMRF patient names all start with MMRF
combined.zscore <- cbind(mmrf.zscore, boston.zscore[-c(1)], arkansas.zscore[-c(1)])
combined.zscore <- cbind(mmrf.zscore, arkansas.zscore[-c(1)])
dim(combined.zscore) # 25 x 1898
combined.zscore[1:5, 1:5]

ttg <- gather(combined.zscore, 
              key="Dataset",
              value="Expression",
              -Genes)
head(ttg)

## change Dataset: from patient names to the dataset name
# MMRF patient names all start with MMRF
# Boston patient names all start with GSM2
# Arkansas patient names start with GSM7 or GSM9
ttg$Dataset[grepl("MMRF", ttg$Dataset)] <- "MMRF"
ttg$Dataset[grepl("GSM2", ttg$Dataset)] <- "Boston"
ttg$Dataset[grepl("GSM7|GSM9", ttg$Dataset)] <- "Arkansas"

# ### genes are ordered alphabetically
# ncols <- 3
# mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(ncols)
# ggplot(ttg, aes(x=Genes, y=Expression, fill=Dataset)) + 
#   geom_boxplot(notch=TRUE, outlier.size=1, show.legend=TRUE) +
#   scale_fill_manual(values=mycolors) +
#   theme(axis.text.x = element_text(angle=90, size=10)) +
#   theme(axis.title = element_text(size=12)) + 
#   ggtitle("Gene Expression") + ylab("Gene expression (z-score)")

### order by gene expression, highest values first
ncols <- 2
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(ncols)
ggplot(ttg, aes(x=reorder(Genes, -Expression), y=Expression, fill=Dataset)) + 
  geom_boxplot(notch=TRUE, outlier.shape=NA, show.legend=TRUE) +
  scale_fill_manual(values=mycolors) +
  theme(axis.text.x = element_text(angle=90, size=10)) +
  theme(axis.title = element_text(size=12)) + 
  ggtitle("Gene Expression of Cell Surface Markers") + 
  ylab("Relative expression (z-score)") + xlab("Genes")

### NOTE: will need this gene order for all plots
genes.order <- reorder(ttg$Genes, -ttg$Expression)
genes.order <- factor(levels(genes.order), levels=levels(genes.order))
saveRDS(genes.order, file="~/Downloads/Rshiny_cell-markers/Rdata/genes21/genes24_order.rds")

### Save plot
tiff("~/Downloads/Rshiny_cell-markers/Rplots/figure1/combined_boxplot_genes24.tiff", res = 300, width = 10, height = 7, units = "in")
ncols <- 3
mycolors <- colorRampPalette(brewer.pal(9, "Blues"))(ncols)
### order by gene expression, highest values first
ggplot(ttg, aes(x=reorder(Genes, -Expression), y=Expression, fill=Dataset)) + 
  geom_boxplot(notch=TRUE, outlier.shape=NA, show.legend=TRUE) +
  scale_fill_manual(values=mycolors) +
  theme(axis.text.x = element_text(angle=90, size=10)) +
  theme(axis.title = element_text(size=12)) + 
  ggtitle("Gene Expression of Cell Surface Markers") + 
  ylab("Relative expression (z-score)") + xlab("Genes")
dev.off()

#####
### Is there significant difference between the datasets?
### aov
## Randomized Block Design (B is the blocking factor)?
fit <- aov(Expression ~ Genes + Dataset, data=ttg)
summary(fit)
## Tukey post-hoc test
TukeyHSD(fit, which="Dataset") ## no sig diff between datasets

