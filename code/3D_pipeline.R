##save to object
DDD.data <- list()
################################################################################
##----->> Set folders to read and save results
data.folder <- file.path(getwd(),'data') # for .RData format results
result.folder <- file.path(getwd(),'result') # for 3D analysis results in csv files
figure.folder <- file.path(getwd(),'figure')# for figures
report.folder <- file.path(getwd(),'report')

DDD.data$data.folder <- data.folder
DDD.data$result.folder <- result.folder
DDD.data$figure.folder <- figure.folder
DDD.data$report.folder <- report.folder

if(!file.exists(data.folder))
  dir.create(path = data.folder,recursive = T)
if(!file.exists(result.folder))
  dir.create(path = result.folder,recursive = T)
if(!file.exists(figure.folder))
  dir.create(path = figure.folder,recursive = T)
if(!file.exists(report.folder))
  dir.create(path = report.folder,recursive = T)

### Set the input data folder
##----->> folder of input files
input.folder <- 'salmon_SRP220383'
quant.folder <- 'salmon_SRP220383/quant'

################################################################################
##----->> parameters of tximport to generate read counts and TPMs
quant_method <- 'salmon' # abundance generator
tximport_method <- 'lengthScaledTPM' # method to generate expression in tximport

################################################################################
##----->> parameters for data pre-processing
### has sequencign replicates?
has_srep <- T

### parameter for low expression filters
cpm_cut <- 1
cpm_samples_n <- 3

### parameter for batch effect estimation
has_batcheffect <- T
RUVseq_method <- 'RUVr' # RUVseq_method is one of 'RUVr', 'RUVs' and 'RUVg'

### data normalisation parameter
norm_method <- 'TMM' ## norm_method is one of 'TMM','RLE' and 'upperquartile'

################################################################################
##----->> parameters for 3D analysis
pval_adj_method <- 'BH'
pval_cut <- 0.01
l2fc_cut <- 1
DE_pipeline <- 'limma'
deltaPS_cut <- 0.1
DAS_pval_method <- 'F-test'

################################################################################
##----->> heatmap
dist_method <- 'euclidean'
cluster_method <- 'ward.D'
cluster_number <- 10

################################################################################
##----->> TSIS
TSISorisokTSP <- 'isokTSP'
TSIS_method_intersection <- method_intersection <- 'mean'
TSIS_spline_df <- spline_df <- NULL
TSIS_prob_cut <- 0.5
TSIS_diff_cut <- 1
TSIS_adj_pval_cut <- 0.05
TSIS_time_point_cut <- 1
TSIS_cor_cut <- 0

################################################################################
##----->> Meta table includes sample information, e.g. conditions, bio-reps, seq-reps, abundance paths, etc.
metatable <- read.csv(file.path(getwd(),'metatable.csv'))
##select the columns of experimental design
factor_col <- c('time','day')
brep_col <- 'bio_rep'
srep_col <- 'seq_rep'
quant_col <- 'quant_files'

##arrange information in the metatable
metatable$label <- as.vector(interaction(metatable[,factor_col]))
metatable$sample.name <- as.vector(interaction(metatable[,c(factor_col,brep_col,srep_col)]))
metatable$quant.folder <-  file.path(quant.folder,metatable$quant_files,
                                     ifelse(quant_method=='salmon','quant.sf','abundance.h5'))

##----->> Transcript-gene association mapping
mapping <-read.csv(file.path(getwd(),'mapping.csv'))
mapping <- data.frame(as.matrix(mapping),stringsAsFactors = F)
rownames(mapping) <- mapping$TXNAME

################################################################################
##----->> Generate gene expression
##
txi_genes <- tximport(metatable$quant.folder,dropInfReps = T,
                      type = quant_method, tx2gene = mapping,
                      countsFromAbundance = tximport_method)

## give colunames to the datasets
colnames(txi_genes$counts) <-
  colnames(txi_genes$abundance) <-
  colnames(txi_genes$length) <-metatable$sample.name

## save the data
write.csv(txi_genes$counts,file=paste0(result.folder,'/counts_genes.csv'))
write.csv(txi_genes$abundance,file=paste0(result.folder,'/TPM_genes.csv'))
save(txi_genes,file=paste0(data.folder,'/txi_genes.RData'))

################################################################################
##----->> Generate transcripts expression
txi_trans<- tximport(metatable$quant.folder, 
                     type = quant_method, tx2gene = NULL,
                     countsFromAbundance = tximport_method,
                     txOut = T,dropInfReps = T)

## give colunames to the datasets
colnames(txi_trans$counts) <- 
  colnames(txi_trans$abundance) <-
  colnames(txi_trans$length) <-metatable$sample.name

## save the data
write.csv(txi_trans$counts,file=paste0(result.folder,'/counts_trans.csv'))
write.csv(txi_trans$abundance,file=paste0(result.folder,'/TPM_trans.csv'))
save(txi_trans,file=paste0(data.folder,'/txi_trans.RData'))

################################################################################
##extract gene and transcript read counts
genes_counts <- txi_genes$counts
trans_counts <- txi_trans$counts
trans_TPM <- txi_trans$abundance

##If no sequencing replicates, genes_counts and trans_counts remain the same by 
if(has_srep){
  idx <- paste0(metatable$label,'.',metatable[,brep_col])
  genes_counts <- sumarrays(genes_counts,group = idx)
  trans_counts <- sumarrays(trans_counts,group = idx)
  metatable_new <- metatable[metatable[,srep_col]==metatable[,srep_col][1],]
} else {
  metatable_new <- metatable
}

################################################################################
##----->> Do the filters
target_high <- low.expression.filter(abundance = trans_counts,
                                     mapping = mapping,
                                     abundance.cut = cpm_cut,
                                     sample.n = cpm_samples_n,
                                     unit = 'counts',
                                     Log=F)
##save expressed genes and transcripts
save(target_high,file=paste0(data.folder,'/target_high.RData'))

################################################################################
##----->> Mean-variance plot
## transcript level

counts.raw = trans_counts[rowSums(trans_counts>0)>0,]
counts.filtered = trans_counts[target_high$trans_high,]
mv.trans <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = metatable_new$label)
### make plot
fit.raw <- mv.trans$fit.raw
fit.filtered <- mv.trans$fit.filtered
mv.trans.plot <- function(){
  par(mfrow=c(1,2))
  plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                   l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
  title('\n\nRaw counts (transcript level)')
  plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                   l = fit.filtered$l,lwd=2,col='black')
  title('\n\nFiltered counts (transcript level)')
  lines(fit.raw$l, col = "gold",lty=4,lwd=2)
  legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
         legend = c('low-exp removed','low-exp kept'))
}
mv.trans.plot()

### save to figure folder
png(filename = paste0(figure.folder,'/Transcript mean-variance trend.png'),
    width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
mv.trans.plot()
dev.off()

pdf(file = paste0(figure.folder,'/Transcript mean-variance trend.pdf'),
    width = 25/2.54,height = 12/2.54)
mv.trans.plot()
dev.off()

################################################################################
## gene level
counts.raw = genes_counts[rowSums(genes_counts>0)>0,]
counts.filtered = genes_counts[target_high$genes_high,]
mv.genes <- check.mean.variance(counts.raw = counts.raw,
                                counts.filtered = counts.filtered,
                                condition = metatable_new$label)
### make plot
fit.raw <- mv.genes$fit.raw
fit.filtered <- mv.genes$fit.filtered
mv.genes.plot <- function(){
  par(mfrow=c(1,2))
  plotMeanVariance(x = fit.raw$sx,y = fit.raw$sy,
                   l = fit.raw$l,lwd=2,fit.line.col ='gold',col='black')
  title('\n\nRaw counts (gene level)')
  plotMeanVariance(x = fit.filtered$sx,y = fit.filtered$sy,
                   l = fit.filtered$l,lwd=2,col='black')
  title('\n\nFiltered counts (gene level)')
  lines(fit.raw$l, col = "gold",lty=4,lwd=2)
  legend('topright',col = c('red','gold'),lty=c(1,4),lwd=3,
         legend = c('low-exp removed','low-exp kept'))
}
mv.genes.plot()

### save to figure folder
png(filename = paste0(figure.folder,'/Gene mean-variance trend.png'),
    width = 25/2.54,height = 12/2.54,units = 'in',res = 300)
mv.genes.plot()
dev.off()

pdf(file = paste0(figure.folder,'/Gene mean-variance trend.pdf'),
    width = 25/2.54,height = 12/2.54)
mv.genes.plot()
dev.off()

################################################################################
##----->> trans level
data2pca <- trans_counts[target_high$trans_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots
groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                groups = groups,plot.title = 'Transcript PCA: bio-reps',
                ellipse.type = ellipse.type,
                add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                groups = groups,plot.title = 'Transcript PCA: average expression',
                ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()


################################################################################
##----->> genes level
data2pca <- genes_counts[target_high$genes_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots

groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                groups = groups,plot.title = 'genescript PCA: bio-reps',
                ellipse.type = ellipse.type,
                add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                groups = groups,plot.title = 'genescript PCA: average expression',
                ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

design <- condition2design(condition = metatable_new$label,
                           batch.effect = NULL)

################################################################################
##----->> trans level
trans_batch <- remove.batch(read.counts = trans_counts[target_high$trans_high,],
                            condition = metatable_new$label,
                            design = design,
                            contrast=NULL,
                            group = metatable_new$label,
                            method = RUVseq_method)
save(trans_batch,file=paste0(data.folder,'/trans_batch.RData')) 

################################################################################
##----->> genes level
genes_batch <- remove.batch(read.counts = genes_counts[target_high$genes_high,],
                            condition = metatable_new$label,
                            design = design,
                            contrast=NULL,
                            group = metatable_new$label,
                            method = RUVseq_method)
save(genes_batch,file=paste0(data.folder,'/genes_batch.RData')) 


################################################################################
## DO the PCA again
################################################################################

##----->> trans level
data2pca <- trans_batch$normalizedCounts[target_high$trans_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots
groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                groups = groups,plot.title = 'Transcript PCA: bio-reps',
                ellipse.type = ellipse.type,
                add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA batch effect removed Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA batch effect removed Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                groups = groups,plot.title = 'Transcript PCA: average expression',
                ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Transcript PCA batch effect removed average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript PCA batch effect removed average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()


################################################################################
##----->> genes level
data2pca <- genes_batch$normalizedCounts[target_high$genes_high,]
dge <- DGEList(counts=data2pca) 
dge <- calcNormFactors(dge)
data2pca <- t(counts2CPM(obj = dge,Log = T))
dim1 <- 'PC1'
dim2 <- 'PC2'
ellipse.type <- 'polygon' #ellipse.type=c('none','ellipse','polygon')

##--All Bio-reps plots
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col] ## colour on biological replicates
# groups <- metatable_new$label ## colour on condtions
g <- plotPCAind(data2pca = data2pca,dim1 = dim1,dim2 = dim2,
                groups = groups,plot.title = 'genescript PCA: bio-reps',
                ellipse.type = ellipse.type,
                add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA batch effect removed Bio-reps.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA batch effect removed Bio-reps.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

##################################################
##--average expression plot
rownames(data2pca) <- gsub('_','.',rownames(data2pca))
groups <- metatable_new[,brep_col]
data2pca.ave <- rowmean(data2pca,metatable_new$label,reorder = F)
groups <- unique(metatable_new$label)
g <- plotPCAind(data2pca = data2pca.ave,dim1 = 'PC1',dim2 = 'PC2',
                groups = groups,plot.title = 'genescript PCA: average expression',
                ellipse.type = 'none',add.label = T,adj.label = F)

g

### save to figure
png(filename = paste0(figure.folder,'/Gene PCA batch effect removed average expression.png'),
    width = 15/2.54,height = 13/2.54,units = 'in',res = 300)
print(g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene PCA batch effect removed average expression.pdf'),
    width = 15/2.54,height = 13/2.54)
print(g)
dev.off()

################################################################################
##----->> trans level
dge <- DGEList(counts=trans_counts[target_high$trans_high,],
               group = metatable_new$label,
               genes = mapping[target_high$trans_high,])
trans_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(trans_dge,file=paste0(data.folder,'/trans_dge.RData'))

################################################################################
##----->> genes level
dge <- DGEList(counts=genes_counts[target_high$genes_high,],
               group = metatable_new$label)
genes_dge <- suppressWarnings(calcNormFactors(dge,method = norm_method))
save(genes_dge,file=paste0(data.folder,'/genes_dge.RData'))

################################################################################
##----->> distribution plot
sample.name <- paste0(metatable_new$label,'.',metatable_new[,brep_col])
condition <- metatable_new$label

###--- trans level
data.before <- trans_counts[target_high$trans_high,]
data.after <- counts2CPM(obj = trans_dge,Log = T)
g <- boxplotNormalised(data.before = data.before,
                       data.after = data.after,
                       condition = condition,
                       sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Transcript expression distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

pdf(file = paste0(figure.folder,'/Transcript expression distribution.pdf'),
    width = 20/2.54,height = 20/2.54)
do.call(grid.arrange,g)
dev.off()

###--- genes level
data.before <- genes_counts[target_high$genes_high,]
data.after <- counts2CPM(obj = genes_dge,Log = T)
g <- boxplotNormalised(data.before = data.before,
                       data.after = data.after,
                       condition = condition,
                       sample.name = sample.name)
do.call(grid.arrange,g)

### save to figure
png(filename = paste0(figure.folder,'/Gene expression distribution.png'),
    width = 20/2.54,height = 20/2.54,units = 'in',res = 300)
do.call(grid.arrange,g)
dev.off()

pdf(file = paste0(figure.folder,'/Gene expression distribution.pdf'),
    width = 20/2.54,height = 20/2.54)
do.call(grid.arrange,g)
dev.off()

RNAseq_info <- data.frame(
  Description=c('Raw transcripts',
                'Raw genes',
                'Samples',
                'Samples after merging seq-reps',
                'Condition of interest',
                'CPM cut-off',
                'Min samples to CPM cut-off',
                'Expressed transcripts',
                'Expressed genes'),
  Number=c(length(mapping$TXNAME),
           length(unique(mapping$GENEID)),
           nrow(metatable),
           nrow(metatable_new),
           length(unique(metatable$label)),
           cpm_cut,
           cpm_samples_n,
           length(target_high$trans_high),
           length(target_high$genes_high))
)
DDD.data$RNAseq_info <- RNAseq_info

RNAseq_info

################################################################################
##----->> pair-wise contrast groups
contrast_pw <- c('T2.Day1-T2.Day0','T2.Day4-T2.Day0','T3.Day1-T3.Day0','T3.Day4-T3.Day0')

##----->> group mean contrast groups
contrast_mean <- c('(T2.Day1+T3.Day1)/2-(T2.Day0+T3.Day0)/2','(T2.Day4+T3.Day4)/2-(T2.Day0+T3.Day0)/2')

##----->> group differences contrast groups
contrast_pgdiff <- c('(T2.Day1-T3.Day1)-(T2.Day0-T3.Day0)','(T2.Day4-T3.Day4)-(T2.Day0-T3.Day0)')

##----->> put together
contrast <- unique(c(contrast_pw,contrast_mean,contrast_pgdiff))

DDD.data$contrast_pw <- contrast_pw
DDD.data$contrast_mean <- contrast_mean
DDD.data$contrast_pgdiff <- contrast_pgdiff
DDD.data$contrast <- contrast
