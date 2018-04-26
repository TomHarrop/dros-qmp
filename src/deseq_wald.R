library(data.table)
library(DESeq2)

BiocParallel::register(BiocParallel::MulticoreParam(8))

# load DSeq object
dds <- readRDS("output/deseq2/dds.Rds")

# subset samples
cd_wald <- colData(dds)[colData(dds)$timepoint == "4872hr",]
counts_wald <- counts(dds)[,rownames(cd_wald)]

# generate batch variable
cd_wald_dt <- data.table(data.frame(cd_wald), keep.rownames = TRUE)
cd_wald_dt[as.numeric(gsub(".*([[:digit:]])$", "\\1", rn)) > 2, batch := "2"]
cd_wald_dt[is.na(batch), batch := "1"]

# regenerate the dds object
dds_wald <- DESeqDataSetFromMatrix(
    countData = counts_wald,
    colData = data.frame(cd_wald_dt, row.names = "rn"),
    design = ~ batch + treatment
)

# run DE analysis
dds_wald <- DESeq(dds_wald,
                  parallel = TRUE)

# extract results
alpha <- 0.1
lfc <- log(1.5, 2)
res_qmp_etoh <- results(dds_wald,
                        contrast = c("treatment", "QMP", "EtOH"),
                        alpha = alpha,
                        lfcThreshold = lfc)
res_stay_etoh <- results(dds_wald,
                         contrast = c("treatment", "stay", "EtOH"),
                         alpha = alpha,
                         lfcThreshold = lfc)
res_stay_qmp <- results(dds_wald,
                        contrast = c("treatment", "stay", "QMP"),
                        alpha = alpha,
                        lfcThreshold = lfc)

# convert to a giant data.table
ConvertResultsToDatatable <- function(x) {
    my_dt <- data.table(data.frame(x), keep.rownames = TRUE)
    setnames(my_dt, "rn", "gene")
    return(my_dt)
}

results_list <- list(qmp_etoh = res_qmp_etoh,
                     stay_etoh = res_stay_etoh,
                     stay_qmp = res_stay_qmp)
results_tables <- lapply(results_list, ConvertResultsToDatatable)
combined_results <- rbindlist(results_tables, idcol = "contrast")

# playtime
combined_results[padj < alpha]
all_sig_genes <- combined_results[padj < alpha, unique(gene)]
counts(dds_wald, normalized = TRUE)[all_sig_genes, ]

# write output
fwrite(combined_results,
       "output/deseq2/wald_results_combined.csv")
fwrite(combined_results[padj < alpha],
       "output/deseq2/wald_results_combined_sig.csv")


