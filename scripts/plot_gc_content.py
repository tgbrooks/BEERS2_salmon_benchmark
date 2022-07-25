import pathlib
import pandas
import seaborn as sns

outdir = pathlib.Path(snakemake.output[0])
outdir.mkdir(exist_ok=True)

summary = pandas.read_csv(snakemake.input.gc_summary, sep="\t")
total_count = summary.groupby(['run', 'sample'])['count'].sum()
total_count.name = 'total_count'
summary = pandas.merge(summary, total_count.reset_index(), on=['run', 'sample'])
summary['normalized_count'] = summary['count'] / summary['total_count']


fig = sns.relplot(
        y = "count",
        x = "GC_content",
        hue = "run",
        col = "sample",
        data = summary,
        kind = "line",
)
fig.savefig(outdir / "summary.png", dpi=300)

fig = sns.relplot(
        y = "normalized_count",
        x = "GC_content",
        hue = "run",
        col = "sample",
        data = summary,
        kind = "line",
)
fig.savefig(outdir / "summary.normalized.png", dpi=300)
