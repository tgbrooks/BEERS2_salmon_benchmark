import pathlib
import pandas
import seaborn as sns

outdir = pathlib.Path(snakemake.output[0])
outdir.mkdir(exist_ok=True)

summary = pandas.read_csv("results/gc_summary.txt", sep="\t")
long_summary = summary.melt(
        id_vars = ['GC_bias', 'pos_3prime_bias', 'sample'],
        var_name = "GC pct",
        value_name = "density",
).reset_index()
long_summary['GC pct'] = long_summary['GC pct'].astype(float)


fig = sns.relplot(
        y = "density",
        x = "GC pct",
        hue = "GC_bias",
        col = "sample",
        row = "pos_3prime_bias",
        data = long_summary,
        kind = "line",
)
fig.savefig(outdir / "summary.png")
