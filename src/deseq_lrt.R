library(DESeq2)

BiocParallel::register(BiocParallel::MulticoreParam(8))

# load DSeq object
dds <- readRDS("output/deseq2/dds.Rds")

# subset samples
cd_lrt <- colData(dds)[colData(dds)$timepoint %in% c("12hr", "24hr", "48hr"),]
counts_lrt <- counts(dds)[,rownames(cd_lrt)]

# filter read counts
filtered_counts <- subset(counts_lrt,
                          rowMeans(counts_lrt) > 10 |
                              rowMax(counts_lrt) > 10)

# regenerate object
dds_lrt <- DESeqDataSetFromMatrix(
    countData = filtered_counts,
    colData = cd_lrt, 
    design = ~ timepoint + treatment + timepoint:treatment)

# run likelihood ratio test
dds_lrt <- DESeq(dds_lrt,
                 test= "LRT",
                 reduced = ~ timepoint + treatment,
                 parallel = TRUE)

# get results
res <- results(dds_lrt, alpha = 0.001)
sig_genes <- rownames(subset(res, padj < 0.001))

# output results
write.csv(subset(res, padj <0.001),
          "output/deseq2/lrt_sig.csv")
write.csv(counts(dds_lrt, normalized = TRUE)[sig_genes, ],
          "output/deseq2/lrt_sig_counts.csv")

# visualize
subset(res, padj <0.001)
vst <- varianceStabilizingTransformation(dds_lrt)
plotPCA(vst, intgroup = c("timepoint", "treatment"))
res[order(res$padj), ]
plotCounts(dds_lrt, "FBgn0014076", intgroup = c("timepoint", "treatment"))
