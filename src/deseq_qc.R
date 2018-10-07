library(DESeq2)

BiocParallel::register(BiocParallel::MulticoreParam(8))

# load DSeq object
dds <- readRDS("output/deseq2/dds.Rds")

# normalise for library size
dds <- estimateSizeFactors(dds)

# transform read counts
vst <- varianceStabilizingTransformation(dds)   # use blind = FALSE for QC

# get transformed counts
vst_assay <- assay(vst)
vst_vars <- rowVars(vst_assay)
vst_sorted <- vst_assay[order(vst_vars),]

# write the sorted VST table
write.csv(vst_sorted, "output/deseq2/vst_sorted.csv")

# does this make sense?
dds$timepoint <- factor(dds$timepoint)
plotCounts(dds, "FBgn0034519", "timepoint")

# plot deseq's PCA
plotPCA(vst, intgroup = c("timepoint", "treatment"), ntop = Inf)


