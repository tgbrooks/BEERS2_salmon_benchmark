import pathlib
import pandas
import seaborn as sns

outdir = pathlib.Path(snakemake.output[0])
outdir.mkdir(exist_ok=True)

summary = pandas.read_csv(snakemake.input.gc_summary, sep="\t")


fig = sns.relplot(
        y = "count",
        x = "GC_content",
        hue = "run",
        col = "sample",
        data = summary,
        kind = "line",
)
fig.savefig(outdir / "summary.png")
