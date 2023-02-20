import pathlib
import pandas
import numpy as np
import matplotlib
import pylab

outdir = pathlib.Path("results/pos_cov/")
outdir.mkdir(exist_ok=True)
genedir = outdir / "geneplots"
genedir.mkdir(exist_ok=True)

cov = pandas.read_csv("results/positional_coverage.txt.gz", sep="\t")

# Position within the gene
t = np.arange(0,1, 0.025)
t_str = t.astype(str)

for gene in ['ENSMUST00000082409', 'ENSMUST00000082407', 'ENSMUST00000184437', 'ENSMUST00000084013', 'ENSMUST00000188299', 'ENSMUST00000082418', 'ENSMUST00000200021', 'ENSMUST00000198477']:
    genecov = cov[(cov.transcript_id == gene) & (cov['sample'] == 1)]
    fig, ax = pylab.subplots(layout="constrained", figsize=(6,3))
    for run, runcov in genecov.groupby("run"):
        ax.plot(
            t,
            runcov[t_str].to_numpy().flatten(),
            label=run,
        )
    ax.set_title(gene)
    ax.set_xlim(0,1)
    ax.set_xlabel("location in transcript")
    ax.set_ylabel("read depth")
    fig.legend(loc="outside right upper", title="run")
    fig.savefig(genedir/f"{gene}.png", dpi=300)
