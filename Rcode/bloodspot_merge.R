library(tidyverse)
library(purrr)
library(RColorBrewer)

directory <- "~/Documents/osu/Perna_lab/bloodspot/HemaExplorer/"

RowToNames <- function(df, row.num = 1) {
  
  col.names <- df[row.num,]
  
  df %<>% set_names(col.names)
  
  df[-row.num,]
}

blood.pooled.sample <-
	list.files(directory, full.names = TRUE, pattern = '*csv') %>%
	map( ~ {
		base <-
			read.delim(.x, sep = ',', stringsAsFactors = FALSE, header = FALSE) %>%
			t %>%
			as_data_frame %>%
			RowToNames

		vals <-
			base %>%
			select(-1) %>%
			mutate_each(funs(as.numeric(.)))

		hgnc <-
			names(base)[[1]] %>%
			str_split(' ') %>%
			unlist %>%
			head(1)

		blood.master <-
			base %>%
			select(1) %>%
			set_names('cell') %>%
			mutate(cell = toupper(cell)) %>%
			mutate(hgnc = hgnc) %>%
			bind_cols(vals) %>%
			arrange(hgnc, cell) %>%
			mutate(cell = ifelse(is.na(cell), 'nan', cell)) %>%
			gather(probe, exp, -cell, -hgnc) %>%
			as_data_frame
	})

str(blood.pooled.sample) # list of 24 genes with expression by different probes

### take the probe with largest mean - function
probe.max <- function(pooled.probes) {
  probe.mean <-  
    pooled.probes %>%
    group_by(probe) %>%
    summarize(mean.exp=mean(exp))
  
  ## find probe name with max mean expression
  max <- probe.mean$probe[which(probe.mean$mean.exp==max(probe.mean$mean.exp))]
  
  # return(max)
  
  ## filter pooled dataset by probe name
  pooled.max <- filter(pooled.probes, probe==max)
    # pooled.probes %>%
    # group_by(probe) %>%
    # filter(probe==max)

  return(pooled.max)
}

### implement - returns a list
# blood.pooled.probes <- vector(mode = "list", length = 25)
blood.pooled.probes <- map(blood.pooled.sample, probe.max)
# lapply(blood.pooled.sample, function(x) probe.max(x))
# map(blood.pooled.sample, probe.max)

### merge list
blood.pooled.probes.merge <- blood.pooled.probes %>%
  data.table::rbindlist(., fill = TRUE) %>%
  unique %>%
  as_tibble()
filter(blood.pooled.probes.merge, cell=="B CELLS")

## save probe list
temp <- blood.pooled.probes.merge %>% select(hgnc, probe) %>% unique
write.csv(temp, file="~/Documents/osu/Perna_lab/bloodspot/HemaExplorer/probe_list.csv", row.names=FALSE)

### this merges and calculates the mean for each gene
blood.pooled.probes.mean <-
  blood.pooled.probes %>%
  map( ~ {
    .x %>%
      group_by(cell, hgnc, probe) %>%
      mutate(exp = mean(exp)) %>%
      unique
  }) %>%
  data.table::rbindlist(., fill = TRUE) %>%
  unique %>%
  tbl_df

### what's the mean expression of CD4+ T CELLS, CD8+ T CELLS, NK CELLS?
#### and which genes are lower?
cd4.cells <- filter(blood.pooled.probes.mean, cell=="CD4+ T CELLS")
a <- filter(cd4.cells, exp < mean(exp))
cd8.cells <- filter(blood.pooled.probes.mean, cell=="CD8+ T CELLS")
b <- filter(cd8.cells, exp < mean(exp))
nk.cells <- filter(blood.pooled.probes.mean, cell=="NK CELLS")
c <- filter(nk.cells, exp < mean(exp))
intersect(intersect(a$hgnc, b$hgnc), c$hgnc) %>% unique
# [1] "CD19"      "CD1D"      "CD200"     "CD70"      "FUT3"      "GPRC5D"    "SDC1"     
# [8] "TNFRSF13B" "TNFRSF17" 

### genes with expression lower than mean exp in HSCs
hsc.cells <- filter(blood.pooled.probes.mean, cell=="HSC_BM")
d <- filter(hsc.cells, exp < mean(exp))
unique(d$hgnc)
# [1] "CD19"      "CD1D"      "CD38"      "CD70"      "CD86"      "FUT3"      "GPRC5D"   
# [8] "ICAM1"     "IGKC"      "ITGB7"     "LY9"       "NCAM1"     "SDC1"      "SLAMF7"   
# [15] "TNFRSF13B" "TNFRSF17" 

dd <- blood.pooled.probes.merge
dd$cell <- factor(dd$cell) # make into factor

### change gene names
dd$hgnc <- factor(dd$hgnc)
levels(dd$hgnc)[levels(dd$hgnc)=="TNFRSF17"] <- "BCMA"
levels(dd$hgnc)[levels(dd$hgnc)=="NCAM1"] <- "CD56"
levels(dd$hgnc)[levels(dd$hgnc)=="SDC1"] <- "CD138"
levels(dd$hgnc)[levels(dd$hgnc)=="TNFRSF13B"] <- "TACI"

### gene order by expression
dd$hgnc <- reorder(dd$hgnc, dd$exp)

###> re-name cell types - TBD
levels(dd$cell)
levels(dd$cell) <- c("B CELLS", "CD14+ MONOCYTES", "CD4+ T CELLS",
                     "CD8+ T CELLS", "CMP", "EARLY HPC GMP", "GMP",
                     "HSC BM", "mDC", "MEP", "MY BM", "NK CELLS",
                     "pDC", "PROMYELOCYTES", "PMN BM", "PMN PB")

saveRDS(dd, "~/Downloads/Rshiny_cell-markers/Rdata/plotting/bloodspot_normal-immune_dd.rds")

### heatmap
ncols <- 10
mycolors <- colorRampPalette(brewer.pal(7, "RdBu"))(ncols)
p <- ggplot(dd, aes(y=hgnc, x=cell, fill=exp)) +
  geom_tile() +
  # theme_pubr() + # doesn't look good when used with coord_fixed()
  coord_fixed() +
  scale_fill_gradientn(colours = rev(mycolors), name="log2(expr)") +
  theme(axis.text.x = element_text(angle=90),
        text=element_text(size=10)) +
  theme(axis.text.y = element_text(face="italic")) +
  xlab("") + 
  ylab("")
p

### save plot
tiff(filename="~/Downloads/Rshiny_cell-markers/Rplots/final_fig3/immune_cells_heatmap_v2.tiff", 
     res = 300, width = 8, height = 10, units = "in")
p
dev.off()


### heatmap version with B cells, T cells, and HSCs
cart.cells <- c("B CELLS", "CD4+ T CELLS", "CD8+ T CELLS", "HSC BM")
dd.plot <- filter(dd, cell == cart.cells)
levels(dd.plot$hgnc)[levels(dd.plot$hgnc)=="CD44"] <- "*CD44"

ncols <- 10
mycolors <- colorRampPalette(brewer.pal(7, "RdBu"))(ncols)
p <- ggplot(dd.plot, aes(y=hgnc, x=cell, fill=exp)) +
  geom_tile() +
  # theme_pubr() + # doesn't look good when used with coord_fixed()
  coord_fixed() +
  scale_fill_gradientn(colours = rev(mycolors), name="log2(expr)") +
  theme(axis.text.x = element_text(angle=90, hjust=0.95, vjust=0.5)) + #, face="italic")) +
  theme(axis.text.y = element_text(face="italic")) +
  xlab("") + 
  ylab("")
p

# p <- ggplot(dd.plot, aes(x=hgnc, y=exp, fill=cell)) +
#   geom_col(position="dodge") +
#   # scale_fill_gradientn(colours = rev(mycolors), name="log2(expr)") +
#   theme(axis.text.x = element_text(angle=90),
#         text=element_text(size=10)) +
#   theme(axis.text.x = element_text(face="italic")) +
#   xlab("") + 
#   ylab("log2(expr)")


### save plot
tiff(filename="~/Documents/osu/Perna_lab/Rshiny_cell-markers/Rplots/final_fig3/selected_immune_cells_heatmap_v1.tiff", 
     res = 300, width = 4, height = 5, units = "in")
p
dev.off()

#############
### boxplot - doesn't look great: too much going on
## choose colors
# ncols <- ifelse(is.factor(dd$cell), length(levels(dd$cell)), 1) # returns num factors or 1
# mycolors <- colorRampPalette(brewer.pal(9, "Greys"))(ncols)
# p <- ggplot(dd, aes(x=hgnc, y=exp, fill=cell)) +
#   # ylim(-10, 20) +
#   geom_boxplot(notch=FALSE, outlier.shape=NA, show.legend=TRUE) +
#   scale_fill_manual(values=mycolors, name=legend.name) +
#   theme_pubr() +
#   theme(legend.position="right") +
#   ylab("log2(Expression)") +
#   xlab("") +
#   theme(axis.text.x = element_text(angle=90, size=10)) #+
#   # facet_grid(~dd$cell, scales="free", space="free") #+
#   # theme(strip.text = element_text(colour = c("#CA0020")))
# p




# ###############################
# blood.pooled.sample.mean <-
#   blood.pooled.sample %>%
#   map( ~ {
#     .x %>%
#       group_by(cell, hgnc, probe) %>%
#       mutate(exp = mean(exp)) %>%
#       unique
#   }) %>%
#   data.table::rbindlist(., fill = TRUE) %>%
#   unique %>%
#   tbl_df 
# 
# blood.pooled.sample.hsc.mean <-
# 	blood.pooled.sample.mean %>%
# 	filter(cell == 'HSC_BM') %>%
# 	rename(hsc.exp = exp) %>%
# 	select(-cell)
# 
# blood.pooled.sample.mean <-
# 	blood.pooled.sample.mean %>%
# 	filter(cell != 'HSC_BM') %>%
# 	full_join(blood.pooled.sample.hsc.mean) %>%
# 	filter(!is.na(hsc.exp))
# 
# blood.pooled.ratio <-
# 	blood.pooled.sample.mean %>%
# 	rowwise %>%
# 	mutate(ratio = exp / hsc.exp) %>%
# 	ungroup
# 
# 
# mutate(gain.amp.fisher = fisher.test( data.frame( a = c(gain.lt.a, gain.gt.eq.a),
# 												  b = c(gain.lt.b, gain.gt.eq.b) ) )$p.value ) %>%
# # p adjustments
# mutate(amp.pos.fisher.adj  = p.adjust(amp.pos.fisher,  'BH')) %>%
# # 0.05 cutoffs
# mutate(amp.pos.fisher.adj.cut  = ifelse(amp.pos.fisher.adj  < 0.05, amp.pos.fisher.adj,  NA)) %>%
# 
# # log
# mutate(amp.pos.fisher.adj.cut.log  = log10(amp.pos.fisher.adj.cut)) %>%
# 
# 
# library(Hmisc)
# 
# cor(mtcars, use="pairwise.complete.obs", method="pearson") 
# 
# blood.pooled.ratio %>%
# split(cell) %>%
# map( ~ {
# 	x = .x$
# 	rcorr(x=1:10, y=11:20, type="pearson")
# })
# 
# 
# #-------------------------------------------------
# 
# list.files('~/Desktop/diff_exp')
# 
# library(limma)
# targets <- readTargets(path = '~/Desktop/diff_exp')
# 
# # 17.3.3 The expression profiles
# 
# x <- read.ilmn(files="~/Desktop/diff_exp/probe profile.txt", ctrlfiles="~/Desktop/diff_exp/control probe profile.txt", other.columns="Detection")
# 
# table(x$genes$Status)
# 
# 
# # 17.3.4 How many probes are truly expressed?
# 
# pe <- propexpr(x)
# dim(pe) <- c(4,3)
# dimnames(pe) <- list(CellType=c("MS","Stroma","ML","LP"),Donor=c(1,2,3))
# 
# 
# # 17.3.5 Normalization and filtering
# 
# y <- neqc(x)
# expressed <- rowSums(y$other$Detection < 0.05) >= 3
# 
# y <- y[expressed,]
# plotMDS(y, labels=targets$CellType)
# 
# 
# # 17.3.6 Within-patient correlations
# 
# ct <- factor(targets$Type)
# 
# design <- model.matrix(~0+ct)
# 
# colnames(design) <- levels(ct)
# 
# dupcor <- duplicateCorrelation(y, design, block = targets$Donor)
# 
# dupcor$consensus.correlation
# 
# 
# # 17.3.7 Differential expression between cell types
# 
# fit <- lmFit(y, design, block = targets$Donor, correlation = dupcor$consensus.correlation)
# contrasts <- makeContrasts(mL-MS, pL-MS, mL-pL, levels = design)
# fit2 <- contrasts.fit(fit, contrasts)
# fit2 <- eBayes(fit2, trend = TRUE)
# summary(decideTests(fit2, method = "global"))
# 
# topTable(fit2, coef=1)
# 
# 
# # 17.3.8 Signature genes for luminal progenitor cells
# 
# ct <- relevel(ct, ref = "pL")
# design <- model.matrix(~ct)
# fit <- lmFit(y, design, block = targets$Donor, correlation = dupcor$consensus.correlation)
# fit2 <- fit[,c("ctMS", "ctmL")]
# fit2 <- eBayes(fit2, trend = TRUE)
# 
# results <- decideTests(fit2, lfc = 1)
# vennDiagram(results, include=c("up", "down"))

