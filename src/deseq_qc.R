library(DESeq2)

BiocParallel::register(BiocParallel::MulticoreParam(8))

# load DSeq object
dds <- readRDS("output/deseq2/dds.Rds")

# normalise for library size
dds <- estimateSizeFactors(dds)

# transform read counts
vst <- varianceStabilizingTransformation(dds)   # use blind = FALSE for QC

# plot deseq's PCA
plotPCA(vst, intgroup = c("timepoint", "treatment"), ntop = Inf)
