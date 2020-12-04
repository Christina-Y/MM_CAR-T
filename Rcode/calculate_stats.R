####################################
### Run ANOVA and Dunnett's Test ###
####################################

### load libraries
library(DescTools)
library(tidyr)
`%notin%` <- Negate(`%in%`)

### load data
mat <- readRDS("~/Documents/osu/Perna_lab/Rshiny_cell-markers/Rdata/mm_pts/genes24_mmrf_v2.rds")
group <- readRDS("~/Documents/osu/Perna_lab/Rshiny_cell-markers/Rdata/factors/mmrf_bcma_factors.rds")

mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/genes24_mmrf.rds")
group <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/factors/hrd_factors_all.rds")
group2 <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/factors/translocations_factors.rds")

### load data GSE6477
mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/normals/genes24_gse6477_log2_full.rds")
group <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/factors/gse6477_norm_all-mm_factors.rds")
levels(group) <- c("Normal donors", "MM patients")
table(group)
mat <- mat[, which(colnames(mat) %in% names(group))]

### choose genes outside of 24 targets
mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/mmrf_tmm_log2.rds")
mat <- mat[which(rownames(mat) %in% c("MAF", "MAFB")), ] ## MAF
g <- c("FCRL1", "FCRL2", "FCRL3", "FCRL4", "FCRL5", "FCRL6", "MDM2", "MDM4", "TP53", "CDKN1A")
g <- c("FCRL1", "FCRL2", "FCRL3", "FCRL4", "FCRL5", "FCRL6")
g <- c("FUT3", "FGFR3", "ICAM1", "MAF", "WHSC1")
mat <- mat[which(rownames(mat) %in% g), ]
## re-level as needed
levels(group)
# group <- factor(group, levels=levels(group)[c(3, 2, 1)])

### HRD data: remove specific translocations
dim(mat)
head(group2)
length(group)
group.rm4.14 <- group[which(names(group) %notin% names(group2)[group2 %in% c("t(4;14)")])]
length(group.rm4.14)
group.rm14.16 <- group[which(names(group) %notin% names(group2)[group2 %in% c("t(14;16)")])]
length(group.rm14.16)
group.rm14.16.rm4.14 <- group[which(names(group) %notin% names(group2)[group2 %in% c("t(14;16)", "t(4;14)")])]
length(group.rm14.16.rm4.14)
group.rm14.16.rm14.20 <- group[which(names(group) %notin% names(group2)[group2 %in% c("t(14;16)", "t(14;20)")])]
length(group.rm14.16.rm14.20)
group.rm14.16.rm14.20.rm4.14 <- group[which(names(group) %notin% names(group2)[group2 %in% c("t(14;16)", "t(14;20)", "t(4;14)")])]
length(group.rm14.16.rm14.20.rm4.14)



### specify file path to save results in:
filepath <- "~/Documents/osu/Perna_lab/Rshiny_cell-markers/Rtables/updated_stats/BCMA_expression_mmrf_v2/"

### match up sample names
mat <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/mm_pts/genes24_mmrf.rds")
group <- group.rm14.16.rm14.20
all.equal(colnames(mat), names(group)) # FALSE
## if FALSE and dimensions of mat is larger than group:
mat <- mat[, which(colnames(mat) %in% names(group))]
## check again
all.equal(colnames(mat), names(group)) # TRUE

######################## only for Arkansas data ########################
###### detour: get Arkansas baseline or relapse group ######
mat <- mat[, which(colnames(mat) %in% names(group)[group=="relapse"])]
group <- readRDS(file="~/Downloads/Rshiny_cell-markers/Rdata/factors/arkansas_relapse_bcma_factors.rds")
levels(group)
### match up sample names
all.equal(colnames(mat), names(group)) # TRUE
########################################################################

############# automated process #############
### wrangle data
d <- data.frame(Genes=rownames(mat), mat)
head(d[, 1:5])
dd <- d %>% gather(Factor, Expression, -Genes)
head(dd)
## Factor right now are sample names, so change to actual factor names
for(fac in levels(group)) {
  dd$Factor[which(dd$Factor %in% names(group)[group==fac])] <- fac
}
dd$Factor <- factor(dd$Factor)
levels(dd$Factor)
## Re-level the factor according to order given in group
dd$Factor <- factor(dd$Factor, levels=levels(group))

### ANOVA or t-test
ttest.pvals <- NULL
for(gene in unique(dd$Genes)){
  dd.gene <- dd[dd$Genes==gene, ]
  
  ## specify statistical test to run based on number of factors:
  if(length(levels(dd.gene$Factor))==2) {
    #> t-test
    test.t <- t.test(dd.gene$Expression ~ dd.gene$Factor)
    test.t$p.value
    
    ttest.pvals <- rbind(ttest.pvals, c(gene, test.t$p.value))
    
  } 
  else {
    #> ANOVA
    test.aov <- aov(dd.gene$Expression ~ dd.gene$Factor)
    aov.pval <- summary(aov(test.aov))[[1]][["Pr(>F)"]][[1]]
    
    # if significant, then Dunnett's test
    # pval is after adjustment for the multiple comparisons
    if(aov.pval <= 0.05) {
      print(paste("yes", gene))
      
      test.dunnett <- DunnettTest(x=dd.gene$Expression, 
                                  g=dd.gene$Factor, 
                                  control=levels(dd.gene$Factor)[1])
      
      write.csv(test.dunnett[[1]], file=paste0(filepath, gene, "_dunnett.csv"))
    }
    
  }
}
### process t-test pval table
if(!is.null(ttest.pvals)) {
  ttest.pvals <- data.frame(Genes=ttest.pvals[, 1],
                            Pval=as.numeric(ttest.pvals[, 2]))
  ttest.pvals$fdr <- p.adjust(ttest.pvals$Pval, method="fdr")
  ttest.pvals
  ttest.pvals$Genes[which(ttest.pvals$fdr < 0.05)]
  
  write.csv(ttest.pvals, file=paste0(filepath, "t_test.csv"), row.names=FALSE)
}

