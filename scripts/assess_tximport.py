import pathlib
import matplotlib
import pylab
import gtfparse
import polars as pl

df = pl.read_parquet("results/tximport/all_values.parquet")

# Load input (true) transcript counts - same for all BEERS2 runs
raw_input_quant = []
for i in range(1,9):
    quants = pl.read_csv(f"data/unbiased/sample{i}/input_quant.txt", sep="\t") \
            .with_columns(pl.lit(i).alias("sample"))
    raw_input_quant.append(quants)
raw_input_quant = pl.concat(raw_input_quant) \
        .with_columns(pl.col('BeersTranscriptID').apply(lambda x: x.split("_")[1]).alias("transcript_id"))
input_quant = raw_input_quant.groupby(["sample", "transcript_id"])\
        .agg(pl.col('count').sum())

# Load the true read counts - different for each BEERS2 run
raw_read_counts = []
for rundir in pathlib.Path("data/").glob("*"):
    if not rundir.is_dir():
        continue
    run = rundir.name
    for i in range(1,9):
        quants = pl.read_csv(rundir / f"sample{i}/output_quant.txt", sep="\t") \
                .with_columns(
                    pl.lit(i).alias("sample"),
                    pl.lit(run).alias("run"),
                )
        raw_read_counts.append(quants)
raw_read_counts = pl.concat(raw_read_counts) \
        .with_columns(pl.col('BeersTranscriptID').apply(lambda x: x.split("_")[1]).alias("transcript_id"))
read_counts = raw_read_counts.groupby(['run', 'sample', 'transcript_id']) \
        .agg(pl.col('count').sum())

# Read the transcript-to-gene mapping
gtf = gtfparse.read_gtf("/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/Mus_musculus.GRCm38.93.gtf.gz")

tx_to_gene = gtf.filter(pl.col('transcript_id') != '') \
        .groupby("transcript_id") \
        .agg(pl.col('gene_id').first().alias("gene_id"))

# Compute the gene-level true read counts
gene_read_counts = read_counts.join(
        tx_to_gene,
        on="transcript_id",
        how="left") \
    .groupby(['run', 'sample', "gene_id"]) \
    .agg(pl.col('count').sum().alias('true_count')) \
    .with_columns(pl.col("sample").cast(pl.Int64).alias("sample_id"))

data = df.join(
        gene_read_counts \
                .select(["run", "sample_id", "gene_id", "true_count"]),
        on=["run", "sample_id", "gene_id"],
        how="left"
    )



# Load the transcript-level quants
salmon_quants = pl.read_parquet("results/salmon_quants.parquet") \
        .rename({"sample": "sample_id", "GeneID": "gene_id"})
# and convert to gene-level by addition
salmon_gene_quants = salmon_quants.groupby(["run", "sample_id", "GC_correct", "Pos_correct", "Seq_correct", "gene_id"]) \
        .agg(pl.col("NumReads").sum()) \
        .join(gene_read_counts, on=["run", "sample_id", "gene_id"])


# Compare to the true values
scaled_corrs = data.groupby(['run', 'sample_id', 'Pos_correct', 'GC_correct', 'Seq_correct']) \
        .agg(pl.spearman_rank_corr("count_scaled", "true_count").alias("corr_scaled"))
lengthScaled_corrs = data.groupby(['run', 'sample_id', 'Pos_correct', 'GC_correct', 'Seq_correct']) \
        .agg(pl.spearman_rank_corr("count_lengthScaled", "true_count").alias("corr_lengthScaled"))
naive_corrs = salmon_gene_quants.groupby(['run', 'sample_id', 'GC_correct', 'Pos_correct', 'Seq_correct']) \
        .agg(pl.spearman_rank_corr("NumReads", "true_count").alias("corr_naive")) \
        .with_columns(
                pl.col("GC_correct").map_dict({True: "yes", False: "no"}),
                pl.col("Pos_correct").map_dict({True: "yes", False: "no"}),
                pl.col("Seq_correct").map_dict({True: "yes", False: "no"}),
        )

corrs = scaled_corrs.join(
        lengthScaled_corrs,
        on = ["run", "sample_id", "Pos_correct", "GC_correct", "Seq_correct"],
        how = "outer",
    ).join(
        naive_corrs,
        on=["run", "sample_id", "Pos_correct", "GC_correct", "Seq_correct"],
        how="outer",
    )

## NOTE:
# after running this comparison
print(corrs)
