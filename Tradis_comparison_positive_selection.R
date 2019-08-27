#!/usr/bin/env Rscript

# PODNAME: tradis_comparison.R
# ABSTRACT: tradis_comparison.R

library("getopt")
options(width=80)
opt = getopt(matrix( c('help', 'h', 0, "logical", 
                       'verbose', 'v', 0, "integer",
                       'controls', 'c', 1, "character",
                       'conditions', 'm', 1, "character",
                       'output', 'o', 1, "character",
                       'plot', 'p', 1, "character",
                       'logFC', 'l', 2, "numeric",
                       'qval', 'q', 2, "numeric",
                       'ins_ratio', 'i', 2, "numeric",
                       'filter', 'f', 0, "logical",
                       'mincount', 't', 1, "integer"
), ncol=4, byrow=TRUE ) );

if(! is.null(opt$help) || is.null(opt$controls )  || is.null(opt$conditions ) || is.null(opt$output) || is.null(opt$plot))
{
  cat(paste("Usage: tradis_comparison.R [-h] [-f] [-t read cutoff] [-l logFC cutoff] [-q q-value cutoff] [-i ins_index ratio cutoff] --controls controls.txt --conditions conditions.txt -o outputfile.csv -p outputplot.pdf \n\n"));
  writeLines(c(strwrap("Compares two experimental conditions using the method of Dembek et al. mBio 2015. Read counts per gene are compared using edgeR. Insertion index ratios between the condition and control are calculated to identify read count changes resulting from a reduced number of insertion sites (likely due to secondary mutations combined with positive selection). This analysis requires experimental replicates."),
               "\n\nRequired Arguments:\n",
               strwrap("--controls : Text file listing the gene-wise insert site csv files for pre-selection (ie 'control') libraries."),
               strwrap("--conditions : Text file listing the gene-wise insert site csv files for libraries exposed to selection."),
               strwrap("-o : output filename"), 
               strwrap("-p : output filename for diagnostic plots"), 
               "\nOptional Arguments:\n",
               strwrap("-l : threshold for log2FC values (default = 2)"),
               strwrap("-q : threshold q-value for significant hits (default = 0.01)"),
               strwrap("-i : threshold insertion index ratio between condition and control for significant hits (default = 0.8)"),
               strwrap("-f : enable filtering on minimum read counts"),
               strwrap("-t : if --filter is enabled, sets minimum read count necessary in one condition for a gene to be included in the comparison."),"\n"))
  q(status=1);
}

if( is.null(opt$filter)) {opt$filter=FALSE}
if( is.null(opt$mincount)) {opt$mincount = 0}
if( is.null(opt$qval)) {opt$qval = 0.01}
if( is.null(opt$logFC)) {opt$logFC = 2}
if( is.null(opt$ins_ratio)) {opt$ins_ratio = 0.8}

#load required packages
library("edgeR")
library("dplyr")

# parse contols and conditions files to lists
control_files <- scan(opt$controls, what="", sep="\n")
condition_files <- scan(opt$conditions, what="", sep="\n")

if(length(control_files) < 2 || length(condition_files) < 2){
  print("2 or more controls/conditions must be provided")
}
if(length(control_files) != length(condition_files)){
  print("Unequal number of conditions and controls provided")
}

control_list = list()
for(i in 1:length(control_files)){
  control_list[[i]] <- read.table(control_files[i], sep=",",header=TRUE, quote="\"", stringsAsFactors=F)
}
condition_list = list()
for(i in 1:length(condition_files)){
  condition_list[[i]] <- read.table(condition_files[i], sep=",",header=TRUE, quote="\"", stringsAsFactors=F)
}

#only look at genes with counts > 0 (or input alternative) in some condition
all_list <- c(control_list, condition_list)

# make list of rows where read count > 0 (or input alternative) in all controls and conditions
read_counts = do.call(cbind, lapply(all_list, function(x){ x$read_count }))

#old case for only 0.
if(! opt$filter){
  zeros = apply( apply(read_counts, 1, ">", 0), 2, any )
} else {
  zeros_cont = apply( apply(read_counts[,1:length(control_files)], 1, ">", opt$mincount), 2, all )
  zeros_cond = apply( apply(read_counts[,(length(control_files) + 1):(length(control_files) + length(condition_files))], 1, ">", opt$mincount), 2, all )
  zeros = (zeros_cont | zeros_cond)
}

# remove these rows
noness_list = lapply(all_list, function(x){ x[zeros,] } )

#build count matrix
count_mat <- do.call(cbind, lapply(noness_list, function(x){x[,7]}))
conds = c()
for(i in 1:length(control_files)){
  conds <- c(conds, "ctrl")
}
for(i in 1:length(condition_files)){
  conds <- c(conds,"cond")
}
conds <- as.factor(conds)

#open plot file
opt$plot = paste(opt$plot,"_plot.pdf", sep = "")
pdf( opt$plot )

#edgeR
d <- DGEList(counts = count_mat, group=conds)
if(length(condition_files) == 3){
  plotMDS.DGEList(d, labels=c("ctrl1","ctrl2","ctrl3","cond1","cond2","cond3"))
} else if(length(condition_files) == 2){
  plotMDS.DGEList(d, labels=c("ctrl1","ctrl2","cond1","cond2"))
} else if(length(condition_files) == 4){
  plotMDS.DGEList(d, labels=c("ctrl1","ctrl2","ctrl3", "ctrl4","cond1","cond2","cond3","cond4"))
}

d <- calcNormFactors(d)
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
de.tgw <- exactTest(d,pair=c("ctrl","cond"))

ctrl1_noness <- noness_list[[1]]
diff <- cbind(ctrl1_noness[,1:2],ctrl1_noness[,11],de.tgw$table,q.value=p.adjust(de.tgw$table$PValue,"BH"))

#volcano plot
plot(diff$logFC, -log(diff$q.value, base=2), xlim=range(c(-6,6)),xlab="Log2 Fold-Change, cond - Ctrl",ylab="-Log2 Q-value", cex = .8, pch=20)
abline(h=-log(opt$qval), col="red")
abline(v=-(opt$logFC), col="red")
abline(v=(opt$logFC), col="red")

#calculate insertion index ratio
ins_index_control = as.data.frame(do.call(cbind, lapply(control_list, function(x){ x$ins_index })))
ins_index_condition = as.data.frame(do.call(cbind, lapply(condition_list, function(x){ x$ins_index })))
ins_index_ratio <- rowMeans(ins_index_condition)/rowMeans(ins_index_control)

#add zero-read genes to results table, add insertion index ratio data and call hits
fullgenelist <- read.csv(control_files[1], sep=",", header = TRUE, stringsAsFactors = FALSE)
fullgenelist <- fullgenelist[,c(1:2,11)]
colnames(fullgenelist)[3] <- "fcn"
colnames(diff)[3] <- "fcn"
results <- right_join(diff, fullgenelist, by = c("locus_tag","gene_name","fcn"))
results$ins_index_ratio <- ins_index_ratio
results <- mutate(results, hit = (logFC <= -(opt$logFC) | logFC >= opt$logFC) & q.value < opt$qval & ins_index_ratio > opt$ins_ratio)

#write results
write.csv(results, file=paste(opt$output,".csv", sep = ""), quote=TRUE, row.names=FALSE)