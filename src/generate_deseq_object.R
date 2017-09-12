#!/usr/bin/env Rscript

library(data.table)
library(DESeq2)
library(rtracklayer)
library(tximport)

gff <- "data/transcriptome/dmel-all-r6.17.gff"

# read gff using GenomicRanges
gr <- import.gff3(gff,
                  feature.type = c("exons", "CDS", "mRNA", "gene"))

# extract a data.frame of tx to gene
mrnas <- gr[gr$type == "mRNA",]
mrna_dt <- as.data.table(mcols(mrnas))
tx2gene <- data.frame(mrna_dt[, .(TXNAME = ID,
                                  GENEID = as.character(Parent))])

# get file list
file_list <- list.files("output/salmon",
                        full.names = TRUE,
                        recursive = TRUE,
                        pattern = "quant.sf")
names(file_list) <- gsub(".*/", "", dirname(file_list))

# import quant files
txi <- tximport(file_list, type = "salmon", tx2gene = data.frame(tx2gene))

# generate col_data
col_data <- data.table(samplename = names(file_list))
col_data[, timepoint := paste0(gsub("^([[:digit:]]+).*", "\\1", samplename),
                               "hr")]
col_data[, treatment := gsub("^[[:digit:]]+([[:alpha:]]+).*", "\\1",
                             samplename)]

# generate DESeq object
dds <- DESeqDataSetFromTximport(
    txi,
    colData = data.frame(col_data, row.names = 'samplename'),
    design = ~ 1
)

# filter read counts
deseq_counts <- counts(dds)
filtered_counts <- deseq_counts[rowMeans(deseq_counts) > 10, ]
dds_filtered <- DESeqDataSetFromMatrix(countData = filtered_counts,
                                       colData = colData(dds),
                                       design = ~ 1
)


# save DESEQ object
saveRDS(dds_filtered, "output/deseq2/dds.Rds")
