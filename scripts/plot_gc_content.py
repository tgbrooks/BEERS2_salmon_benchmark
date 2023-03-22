import pathlib
import pandas
import seaborn as sns

outdir = pathlib.Path(snakemake.output[0])
outdir.mkdir(exist_ok=True)

run_configs = snakemake.config['run_configs']
run_order = run_configs.keys()

summary = pandas.read_csv(snakemake.input.gc_summary, sep="\t")
total_count = summary.groupby(['run', 'sample'])['count'].sum()
total_count.name = 'total_count'
summary = pandas.merge(summary, total_count.reset_index(), on=['run', 'sample'])
summary['normalized_count'] = summary['count'] / summary['total_count']

fig = sns.relplot(
        y = "count",
        x = "GC_content",
        hue = "run",
        hue_order = run_order,
        col = "sample",
        data = summary,
        kind = "line",
)
fig.savefig(outdir / "summary.png", dpi=300)

fig = sns.relplot(
        y = "normalized_count",
        x = "GC_content",
        hue = "run",
        hue_order = run_order,
        col = "sample",
        data = summary,
        kind = "line",
)
fig.savefig(outdir / "summary.normalized.png", dpi=300)

fig = sns.relplot(
        y = "count",
        x = "GC_content",
        hue_order = run_order,
        hue = "run",
        data = summary.query("sample == 1"),
        kind = "line",
)
fig.savefig(outdir / "sample1.summary.png", dpi=300)

fig = sns.relplot(
        y = "normalized_count",
        x = "GC_content",
        hue = "run",
        hue_order = run_order,
        data = summary.query("sample == 1"),
        kind = "line",
)
fig.savefig(outdir / "sample1.summary.normalized.png", dpi=300)
