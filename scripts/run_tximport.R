library(tximport)
library(tibble)
library(tidyr)
library(dplyr)
library(readr)
library(GenomicFeatures)

print(sessionInfo())
options(warn = 1, keep.source = TRUE, error = quote({
      # Debugging in R
      #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/index.shtml
      #
      # Post-mortem debugging
      #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/pmd.shtml
      #
      # Relation functions:
      #   dump.frames
      #   recover
      # >>limitedLabels  (formatting of the dump with source/line numbers)
      #   sys.frame (and associated)
      #   traceback
      #   geterrmessage
      #
      # Output based on the debugger function definition.

      # TODO: setup option for dumping to a file (?)
      # Set `to.file` argument to write this to a file for post-mortem debugging    
      dump.frames()  # writes to last.dump
        n <- length(last.dump)
        if (n > 0) {
                calls <- names(last.dump)
            cat("Environment:\n", file = stderr())
                cat(paste0("  ", seq_len(n), ": ", calls), sep = "\n", file = stderr())
                cat("\n", file = stderr())
                  }

          if (!interactive()) q()
}))

gtf <- "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf"
files <- snakemake@input[['quants']]
sample_ids <- snakemake@params[['sample_ids']]

# Prepare the transcript -> gene map
txdb <- makeTxDbFromGFF(gtf)
tx2gene <- transcripts(txdb,columns=c("TXNAME", "GENEID"))
tx2gene$GENEID <- unlist(tx2gene$GENEID)
tx2gene_df <- as_tibble(tx2gene)[,c("TXNAME", "GENEID")]

# Run tximport
# with both settings of countsFromAbundance
txi_scaled <- tximport(files, type = "salmon", tx2gene = tx2gene_df, ignoreTxVersion=TRUE, countsFromAbundance = "scaledTPM")
txi_lengthScaled <- tximport(files, type = "salmon", tx2gene = tx2gene_df, ignoreTxVersion=TRUE, countsFromAbundance = "lengthScaledTPM")

# Add the sample ids back
num_reads_scaled <- txi_scaled$counts
colnames(num_reads_scaled) <- sample_ids
num_reads_scaled <- as_tibble(num_reads_scaled) |>
    mutate(gene_id = rownames(txi_scaled$counts)) |>
    pivot_longer(!gene_id, names_to="sample_id", values_to="count_scaled")

num_reads_lengthScaled <- txi_lengthScaled$counts
colnames(num_reads_lengthScaled) <- sample_ids
num_reads_lengthScaled <- as_tibble(num_reads_lengthScaled) |>
    mutate(gene_id = rownames(txi_lengthScaled$counts)) |>
    pivot_longer(!gene_id, names_to="sample_id", values_to="count_lengthScaled")

num_reads <- inner_join(num_reads_scaled, num_reads_lengthScaled, by=c("gene_id", "sample_id"))

# Gather the results
write_tsv(num_reads, snakemake@output[['quant_file']])
