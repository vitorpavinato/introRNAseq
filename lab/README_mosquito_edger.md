# Introduction to RNAseq analysis in R with edgeR

For this exercise we are going to analyze the RNAseq data from mosquitoes fed with blood and compare to the expression of mosquitoes fed with sugar. The idea is to find differently expressed genes when mosquitoes are fed with blood. 

The count data was already normalized.

Load packages
```{r load_packages, echo=FALSE, warning=FALSE, message=FALSE}
library(edgeR)
library(gplots)
library(tidyverse)
```

Import the data that is located at `data/` into R by running this code below (copy and paste the code into R console):
```{r}
infile = "lab/data/normalized_counts.csv"
counts = read.table(infile, header=TRUE, sep="\t", row.names=1 )
```

**Question 1**: How many features the dataset contain?

**Question 2**: Check the first 5 rows. How many columns the dataset has? Could you see the count data for each replication?

Read in Annotation
```{r read_annotation}
anno <- read.delim("lab/data/ensembl_mg_52.txt",as.is=T)
dim(anno) #211488     17
head(anno)
tail(anno)
# Because the features selected to represent the expression data was exons (AAEL000004-RA-E1), we need to have a column with exons IDs to compare.
any(duplicated(anno$Exon.stable.ID))
```

Run the code below to merge the normalized read count with the annotation table
```{r}
counts_anno <- merge(counts, anno, by.x=1, by.y=15, all.x = TRUE, all.y)
counts_anno <- counts_anno[,c(1:6, 10)]

counts_anno.aggr <- counts_anno %>%
  group_by(Gene.stable.ID) %>%
  summarise(
    Sugar_1 = mean(sugar_1, na.rm=T),
    Sugar_2 = mean(Sugar_2, na.rm=T),
    Blood_1 = mean(Blood_1, na.rm=T),
    Blood_2 = mean(Blood_2, na.rm=T)
  )

Gene_IDs <- counts_anno.aggr$Gene.stable.ID
counts_anno.aggr <- as.matrix(counts_anno.aggr[,-c(1)])
row.names(counts_anno.aggr) <- Gene_IDs
head(counts_anno.aggr)

counts_anno.aggr <- counts_anno.aggr[complete.cases(counts_anno.aggr), ]
counts_anno.aggr <- counts_anno.aggr[!is.na(Gene_IDs), ]
head(counts_anno.aggr)
```

Run the code below to define the design matrix for the analysis.
For this dataset, the design matrix should be 2x2 (2 replicates of "blood" and 2 replicates of "sugar")
```{r}
reps1 = as.integer(2)
reps2 = as.integer(2)

# Set up the conditions based on the experimental setup.
sugar = rep("sugar", reps1); print(sugar)
blood = rep("blood", reps2); print(blood)
```
Run the code below to keep only the count data from the input (remove the additional information):
```{r}
# Create the groups.
group=c(sugar, blood)

head(counts_annot.aggr, n = 5)
```

**Question 3**: Take a look at `group`. How does the data is organized? What does it represent?


Creates a DGEList object from a table of counts and group.
```{r}
dge0 <- DGEList(counts=counts_annot.aggr, group=group)
```

Here we will filter low-expressed genes, remove any row (gene) whose max value (for the row) is less than the cutoff (2).
```{r filter}
cutoff <- 2
drop <- which(apply(cpm(dge0), 1, max) < cutoff)
dge <- dge0[-drop,]
dim(dge) # number of genes left
```

Visualizing the data with a Multidimensional scaling (MDS) plot (**PCA Plot**)
```{r mds, fig.width=6}
plotMDS(dge, col = c(rep("red",2), rep("blue",2)), cex=1)
```

Extract the edgeR results
```{r}
# Maximizes the negative binomial conditional common likelihood to estimate a common dispersion value across all genes.
dis <- estimateCommonDisp(dge)

# Estimates tagwise dispersion values by an empirical Bayes method based on weighted conditional maximum likelihood.
tag <- estimateTagwiseDisp(dis)

# Compute genewise exact tests for differences in the means between the groups.
etx <- exactTest(tag)

# Extracts the most differentially expressed genes.
etp <- topTags(etx, n=nrow(counts))
```

Run the code below to get the **Volcano Plot**
```{r}
alpha <- 0.01 # Threshold on the adjusted p-value
cols <- densCols(etp$table$logFC, -log10(etp$table$PValue))
plot(etp$table$logFC, -log10(etp$table$FDR), col=cols, panel.first=grid(),
     main="Volcano plot", xlab="Effect size: log2(fold-change)", ylab="-log10(adjusted p-value)",
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-2,2), col="brown")
abline(h=-log10(alpha), col="brown")
```
Run the code below to get the normalized counts.
```{r}
# Get the scale of the data
scale = dge$samples$lib.size * dge$samples$norm.factors

# Get the normalized counts
normed = round(t(t(counts_annot.aggr)/scale) * mean(scale))
```

**Question 4**: What does it mean normalized data? Why does we care about it? Could we use the raw counts to make comparisons? Is different than the normalization already applied to the data?

Run the code below to turn the edgeR results into data frames.
```{r}
# Select the edger result dataframe.
data = etp$table

# Create column placeholders.
data$baseMean = 1
data$baseMeanA = 1
data$baseMeanB = 1
data$foldChange = 2 ^ data$logFC

# Rename the column.
names(data)[names(data)=="logFC"] <-"log2FoldChange"

# Compute the adjusted p-value
data$PAdj = p.adjust(data$PValue, method="hochberg")

# Reorganize the columns. Make errPct the last column before normalized data.
data = data[c(5,6,7,8,1,2,3,9,4)]

# Sort the data by PValue to compute false discovery counts.
data = data[with(data, order(PValue, -foldChange)), ]

# Compute the false discovery counts on the sorted table.
data$falsePos = 1:nrow(data) * data$FDR
```

Merge the two datasets by row names
```{r}
# Create a merged output that contains the normalized counts.
total <- merge(data, normed, by='row.names')

# Sort the data for the output. 
total = total[with(total, order(PValue, -foldChange)), ]

# Rename columns for consistency with other methods.
colnames(total)[1] <- "name"
```

Run the code to get the final output
```{r}
# Start index for normalized data for condition 1.
start1 = 12

# End index for normalized data for condition 1.
end1 = start1 + reps1 - 1

# Start index for normalized data for condition 2.
start2 = end1 + 1

# End index for normalized data for condition 2.
end2 = start2 + reps2 - 1

total$baseMean = rowMeans(total[start1:end2])
total$baseMeanA = rowMeans(total[,start1: end1])
total$baseMeanB = rowMeans(total[,start2: end2])

# Round the numbers
total$foldChange = round(total$foldChange, 3)
total$FDR = round(total$FDR, 4)
total$PAdj = round(total$PAdj, 4)
total$logCPM = round(total$logCPM, 1)
total$log2FoldChange = round(total$log2FoldChange, 1)
total$baseMean = round(total$baseMean, 1)
total$baseMeanA = round(total$baseMeanA, 1)
total$baseMeanB = round(total$baseMeanB, 1)
total$falsePos = round(total$falsePos, 0)
```

Run the code to export the output to `lab/results` folder.
```{r}
outfile = "lab/results/results_mosquito_edger_genes.csv"
write.csv(total, file=outfile, row.names=FALSE, quote=FALSE)
```

**Heatmap for edgeR results**

Run the code to draw a heatmap from the output that contains a normalized matrix
```{r}
# Set the plot dimensions.
WIDTH = 12
HEIGHT = 13

# Default FDR cutoff of 5%
LIMIT = 5

# Input values should be in percent!
LIMIT = LIMIT/100

# Set the margins
MARGINS = c(9, 12)

# Relative heights of the rows in the plot.
LHEI = c(1, 5)

# Read normalized counts from the standard input
normMatrix = "lab/results/results_mosquito_edger_genes.csv"
data = read.csv(normMatrix, header=T, as.is=TRUE)

# Subset data for values under a threshold.
data = subset(data, data$FDR <= LIMIT)

# The heatmap row names will be taken from the first column.
row_names = data[, 1]

# The normalized data starts past the rightmost of these columns.
idx = which(colnames(data) == "falsePos") + 1

# The normalized counts are on the right size.
counts = data[, idx : ncol(data)]

# Load the data from the second column on.
values = as.matrix(counts)

# Adds a little noise to each element to avoid the clustering function failure on zero variance rows.
values = jitter(values, factor = 1, amount = 0.00001)

# Normalize each row to a z-score
zscores = NULL
for (i in 1 : nrow(values)) {
    row = values[i,]
    zrow = (row - mean(row)) / sd(row)
    zscores = rbind(zscores, zrow)
}

# Set the row names on the zscores.
row.names(zscores) = row_names

# Turn the data into a matrix for heatmap2.
zscores = as.matrix(zscores)
```

Plot
```{r}
pdf("lab/results/heatmap_mosquito_edger_genes.pdf", width = WIDTH, height = HEIGHT)
# Set the color palette.
col = greenred

# Draw the heatmap.
heatmap.2(zscores, col=col, density.info="none", Colv=NULL,
          dendrogram="row", trace="none", margins=MARGINS, lhei=LHEI
          )
# Turn off the device.
dev.off()
```

Take a look at the heatmap (open the pdf at `lab/results/`)

**Final questions**:

1. How many genes are **ACTUALLY** differently expressed between the two conditions? We are going to assume that differently expressed genes (DEG) are those with FDR < 0.05. This is not the only possible way to look at DEG. You can, for example, use the values in the adjusted pvalue column (PAdj).

You can run this code in R to find out:

```{r}
write.table(total[which(total$FDR < 0.05), "name"], file="lab/results/results_mosquito_edger_genes_signf_names.txt", 
          row.names=FALSE, col.names = FALSE, quote=FALSE)
```
How about only significant DEG but upper-regulated in sugar-fed mosquito?
```{r}
write.table(total[which(total$FDR < 0.05 & total$log2FoldChange > 0), "name"], file="lab/results/results_mosquito_edger_genes_signf_upper_names.txt", 
          row.names=FALSE, col.names = FALSE, quote=FALSE)
```

How about only significant DEG but down-regulated in sugar-fed mosquito?
```{r}
write.table(total[which(total$FDR < 0.05 & total$log2FoldChange < 0), "name"], file="lab/results/results_mosquito_edger_genes_signf_down_names.txt", 
          row.names=FALSE, col.names = FALSE, quote=FALSE)
```

2. Can you find where the differently exons came from (identify the parental gene)? **Hint:** Use the annotation file.

```{r}
counts_anno_geneIDs <- counts_anno[,c(6,7)]
counts_anno_geneIDs <- counts_anno_geneIDs %>% distinct()

total_anno <- merge(total, counts_anno_geneIDs, by.x='name', by.y="Gene.stable.ID", sort = FALSE)



write.csv(total_anno[which(total_anno$FDR < 0.05), ], 
          file="lab/results/results_mosquito_edger_genes_signf_table.csv", 
          row.names=FALSE, quote=FALSE)
```
How about the information of upper-regulated genes in mosquito fed with sugar?
```{r}
write.csv(total_anno[which(total_anno$FDR < 0.05 & total_anno$log2FoldChange > 0), ], 
          file="lab/results/results_mosquito_edger_genes_signf_upper_table.csv", 
          row.names=FALSE, quote=FALSE)
```

And about the information of down-regulated genes in mosquito fed with sugar?
```{r}
write.csv(total_anno[which(total_anno$FDR < 0.05 & total_anno$log2FoldChange < 0), ], 
          file="lab/results/results_mosquito_edger_genes_signf_down_table.csv", 
          row.names=FALSE, quote=FALSE)
```

```{r}
sessionInfo()
```