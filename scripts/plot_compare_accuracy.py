import pathlib

import pandas
import seaborn as sns

#run_configs = snakemake.config['run_configs']
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

outdir = pathlib.Path(snakemake.output.outdir)
outdir.mkdir(exist_ok=True)

stats = pandas.read_csv(snakemake.input.stats_file, sep="\t")

stats_keys = [c for c in stats.columns
                if c not in ['sample', 'run', 'correction type']]

correction_type_order = ['baseline', 'Seq', 'Pos', 'Pos Seq', 'GC', 'GC Seq', 'GC Pos', 'GC Pos Seq']

#### PLOTS ####
for stat in stats_keys:
    # Plot correlation summary stats
    fig = sns.catplot(
        data=stats,
        x="correction type",
        y=stat,
        order = correction_type_order,
        col = "run",
        col_order = run_order,
        col_wrap = 3,
        sharey=False,
    )
    [tick.set_rotation(90) for ax in fig.axes.flatten() for tick in ax.get_xticklabels()]
    fig.savefig(outdir / f"{stat}.png", dpi=300)
