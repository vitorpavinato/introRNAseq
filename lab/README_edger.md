# Introduction to RNAseq analysis in R with edgeR

The code for edgeR and heatmaps were modified from [Biostar handbook: RNAseq by example](http://www.biostarhandbook.com).

This is an introduction to the analysis of RNAseq data using R. We have assumed that the count data was already produced by a pipeline. Examples on how produce count data from raw fastq reads from whole transcriptome sequence will be show in lab. There are two main ways to get count data: 1) by mapping reads to a reference genome and 2) by mapping reads to a transcriptome. 

At some specific points, there are some questions to be answered. At the end, will should be able to identify the deferentially expressed genes between the two conditions. 

The example was extracted from the [Biostar handbook: RNAseq by example](http://www.biostarhandbook.com). There it was used to show how different methods for the analysis of RNAseq data can sometimes give you different results. Don't worry, it might happen most of the times in bioinformatics. Different methods were developed with different assumptions in mind. To be a good bioinformatician is to know that is going behind the scene and make reasonable decision about the results obtained with different methods.

We are going to analyze the data from the **Golden Snidget**. "The Golden Snidget, is a magical golden bird with fully rotational wings; it will help you understand how RNA-Seq data analysis works, from its strengths to its limitations."

"The Golden Snidget is not a toy model! You will see that the same code and approach used to analyze the Golden Snidget data can be used almost unchanged in many other realistic data analysis scenarios."


### **What do we know about the Golden Snidget?**

Partly mechanical and fully fictitious, the Golden Snidget turns out to be the perfect candidate for RNA-Seq data analysis.

![image 1](images/golden.png)

What makes it so well suited? Two peculiar features. First, the Golden Snidget can only be in two moods:
1. BORED when it flies nice and straight.
2. EXCITED when it constantly changes direction.

Additionally, the Golden Snidget is so magically regulated that, in each state, we know precisely how many transcripts are expressed. The organism is so perfectly controlled that even genes could be named in a way that describes their changes for both states. For example, a gene might be called:

```bash
ABC-1000-UP-4
```

The gene was named as such because in the BORED state, at every moment, each cell has 1000 copies of the transcript corresponding to gene ABC. More- over, the name also conveys that in the EXCITED state, this same gene is be up-regulated 4 fold compared to the BORED state. We call state where the mRNA concentrations do not change ( the synthesis and the degrada- tion rates are identical) as ???steady state mRNA???. In the Golden Snitch all transcripts are present in steady state.

Another way to explain the naming scheme is that in the BORED state, there will be 1000 copies of the ABC transcript, whereas in the EXCITED state, each cell will contain 4000 copies of the ABC transcript.

Thus, just by looking at a gene name, we already know how much of it is inside the cells. We can create a quantification matrix for the real transcript abundance inside the cell from the name of the genes. Our data would look like this:

```bash
|feature name  |BORED    |EXCITED    |foldChange|
|--------------|---------|-----------|----------|
|ABC-1000-UP-4 |     1000|       4000|         4|
|ABD-40-DOWN-2 |       40|         20|         2|
|ABE-90-SAME-1 |       90|         90|         1|
 ...
```

Across the two states each gene can only be either UP, DOWN regulated or stay the SAME. Note how the quantification matrix above represents what ???really??? happens in the cell.

The whole point of the RNA-Seq analysis is to find the counts matrix. As it happens, for the Golden Snidget we already know what the count matrix ought to be. It remains to be seen how well results from an RNA-Seq data can reproduce the expected numbers or ratios. The unique gene naming convention will allow you to spot check any intermediate result and evaluate how well the analysis is going even before you get to the final results.

Exercise extracted from: [Biostar handbook: RNAseq by example](http://www.biostarhandbook.com).

Load libraries
```{r}
library(edgeR)
library(gplots)
```

Import the data that is located at `data/` into R by running this code below (copy and paste the code into R console):
```{r}
infile = "lab/data/counts.txt"
counts = read.table(infile, header=TRUE, sep="\t", row.names=1 )
```

**Question 1**: How many features the dataset contain?

**Question 2**: Check the first 5 rows. How many columns the dataset has? Could you see the count data for each replication?

**Question 3**: Are there transcripts from more than one strand?

Run the code below to define the design matrix for the analysis.
For this dataset, the design matrix should be 3x3 (3 replicates of condition 1 and 3 replicates of condition 2)
```{r}
reps1 = as.integer(3)
reps2 = as.integer(3)

# Set up the conditions based on the experimental setup.
cond1 = rep("cond1", reps1); print(cond1)
cond2 = rep("cond2", reps2); print(cond2)
```

Run the code below to keep only the count data from the input (remove the additional information):
```{r}
# Assume the last N + M columns are the count matrix.
idx = ncol(counts) - (reps1 + reps2)

# Cut out the valid columns.
counts = counts[-c(1:idx)]

# Create the groups.
group=c(cond1, cond2)

head(counts, n = 5)
```

**Question 4**: Take a look at `group`. How does the data is organized? What does it represent?

Creates a DGEList object from a table of counts and group.
```{r}
dge <- DGEList(counts=counts, group=group)
```

Visualizing the data with a Multidimensional scaling (MDS) plot (**PCA Plot**)
```{r mds, fig.width=6}
plotMDS(dge, col = c(rep("red",3), rep("blue",3)), cex=1)
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
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
```

Run the code below to get the normalized counts.
```{r}
# Get the scale of the data
scale = dge$samples$lib.size * dge$samples$norm.factors

# Get the normalized counts
normed = round(t(t(counts)/scale) * mean(scale))
```

**Question 5**: What does it mean normalized data? Why does we care about it? Could we use the raw counts to make comparisons?

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
outfile = "lab/results/results_edger.csv"
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
normMatrix = "lab/results/results_edger.csv"
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
pdf("lab/results/heatmap_edger.pdf", width = WIDTH, height = HEIGHT)
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

1. How many genes are **EXPECTED** to be  differently expressed between the two conditions?

Hint: Try to count how many genes in `results/results_edger.csv` has `UP` or `DOWN` in its name. You can open the file in excel and count, or you can go the the terminal and type:

```bash
cat results_edger.csv | egrep "UP|DOWN" | wc -l
```

2. How many genes are **ACTUALLY** differently expressed between the two conditions? We are going to assume that differently expressed genes (DEG) are those with FDR < 0.05. This is not the only possible way to look at DEG. You can, for example, use the values in the adjusted pvalue column (PAdj).

3. Does all genes were correctly identified by edgeR as differently expressed? 

3. How many edgeR missed? 

4. Why did it happens?

5. How many false positives? A false positive here is a gene that was erroneously identified as DEG. Here, it is easy to spot a false positive gene because we can count how many of genes with the tag `_SAME_` in its name on the list of genes with FDR < 0.05. 