import polars as pl
import gtfparse
import seaborn as sns

sample_ids = list(range(1,9))

bias_order = ["none", "med", "high"]
#run_order = run_configs.keys()
run_order = [
    'unbiased',
    'GC_bias_med',
    'GC_bias_high',
    'pos_3prime_bias_med',
    'pos_3prime_bias_high',
    'primer_bias_high',
    'all_bias',
]
MIN_NUM_READS = 100

out_dir = pathlib.Path("results/compare_gene_accuracy")
out_dir.mkdir(exist_ok=True)

# Load BEERS data and salmon vlaues
input_files = dict(
        salmon = "results/salmon_quants.parquet",
        beers_in = "results/beers_input_quants.parquet",
)
salmon = pl.read_parquet(input_files['salmon']) \
            .with_columns(
                    pl.when(pl.col("GeneID").is_null())
                        .then(pl.col("TranscriptID"))
                        .otherwise(pl.col("GeneID"))
                        .alias("GeneID")
            )
beers_in = pl.read_parquet(input_files['beers_in'])

# Read the transcript-to-gene mapping
gtf = gtfparse.read_gtf("/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/Mus_musculus.GRCm38.93.gtf.gz")

tx_to_gene = gtf.filter(pl.col('transcript_id') != '') \
        .groupby("transcript_id") \
        .agg(pl.col('gene_id').first().alias("gene_id"))

# Map transcript IDs according to salmon's collapsing of identical transcripts
# Otherwise, comparison is unfair - can't possibly tell these apart
# DuplicateRef -> RetainedRef
salmon_collapse = pl.read_csv("/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/index/duplicate_clusters.tsv", sep="\t") \
        .with_columns(
                pl.col("DuplicateRef").str.split(".").arr.get(0).alias("TranscriptID"),
                pl.col("RetainedRef").str.split(".").arr.get(0).alias("RetainedRef"),
        ) \
        .select(["TranscriptID", "RetainedRef"])
beers_in_collapsed = beers_in.join(
        salmon_collapse,
        on = "TranscriptID",
        how = "left"
        ).select(
            "run",
            "sample",
            pl.when(pl.col("RetainedRef").is_not_null()).then(pl.col("RetainedRef")).otherwise(pl.col("TranscriptID")).alias("TranscriptID"),
            "count",
        )


# Aggregate gene-level across all the transcripts (and alleles)
gene_beers_in = (beers_in_collapsed.join(
            tx_to_gene.rename({"transcript_id": "TranscriptID", "gene_id": "GeneID"}),
            on = "TranscriptID",
            how = "left",
        )
        .with_columns( pl.col("GeneID").fill_null(pl.col("TranscriptID")) )
        #.select(pl.exclude("BeersTranscriptID"))
        .groupby(['GeneID', 'run', 'sample'])
        .agg([
            pl.col("count").sum().alias("true_count")
        ]))

# input molecule counts give the underlying true TPM, after normalizing
true_tpm = gene_beers_in.with_columns([
    (1e6*pl.col('true_count') / pl.col('true_count').sum()).over(['run', 'sample']).alias('true_tpm'),
])

salmon_gene = salmon.groupby(["GeneID", "run", "sample", "GC_correct", "Pos_correct", "Seq_correct"]) \
        .agg([
            pl.col("NumReads").sum(),
            pl.col("TPM").sum(),
        ])

data = salmon_gene.join(
        true_tpm,
        on = ["GeneID", "run", "sample"],
        how = "outer"
    )
# A small number of genes do not appear in Salmon output, but the quants are low
# so we drop them (they have nulls) and
# genes with no expression do not appear in BEERS quant files, so fill with zeros
data = data.drop_nulls("TPM") \
        .with_columns(
            pl.col("true_tpm").fill_null(0).keep_name()
        )

corr = data.groupby(["run", "sample", "GC_correct", "Pos_correct", "Seq_correct"])\
        .agg(pl.spearman_rank_corr("true_tpm", "TPM").alias("corr"))
def MARD(col1, col2):
    A = pl.col(col1)
    B = pl.col(col2)
    return ((A - B).abs() /(A + B).abs()).drop_nans().mean()
mard = data.groupby(["run", "sample", "GC_correct", "Pos_correct", "Seq_correct"])\
        .agg(MARD("true_tpm", "TPM").alias("MARD"))

## Draw plots
## Compare correlation no-correct vs all-correct
comparison = pl.concat([
    corr.filter(pl.col("GC_correct") & pl.col("Pos_correct") & pl.col("Seq_correct"))
        .select("run", "sample", "corr", pl.lit(True).alias("corrected")),
    corr.filter(~pl.col("GC_correct") & ~pl.col("Pos_correct") & ~pl.col("Seq_correct"))
        .select("run", "sample", "corr", pl.lit(False).alias("corrected")),
    ])

fig = sns.catplot(
    comparison.to_pandas(),
    y = "run",
    x = "corr",
    hue = "corrected",
    order = run_order,
    kind = "box",
)
sns.swarmplot(comparison.to_pandas(), y="run", x="corr", hue="corrected", order=run_order, palette="dark:k", size=3, dodge=True, legend=False)
fig.savefig(out_dir / "corrected_vs_uncorrected.spearman_correlation.png", dpi=300)

## Compare MARD no-correct vs all-correct
mard_comparison = pl.concat([
    mard.filter(pl.col("GC_correct") & pl.col("Pos_correct") & pl.col("Seq_correct"))
        .select("run", "sample", "MARD", pl.lit(True).alias("corrected")),
    mard.filter(~pl.col("GC_correct") & ~pl.col("Pos_correct") & ~pl.col("Seq_correct"))
        .select("run", "sample", "MARD", pl.lit(False).alias("corrected")),
    ])

fig = sns.catplot(
    mard_comparison.to_pandas(),
    y = "run",
    x = "MARD",
    hue = "corrected",
    order = run_order,
    kind = "box",
)
sns.swarmplot(mard_comparison.to_pandas(), y="run", x="MARD", hue="corrected", order=run_order, palette="dark:k", size=3, dodge=True, legend=False)
fig.savefig(out_dir / "corrected_vs_uncorrected.MARD.png", dpi=300)

