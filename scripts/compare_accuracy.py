import pathlib

import numpy
import scipy.stats
import pandas
import seaborn as sns

out_dir = pathlib.Path(snakemake.output.dir)
out_dir.mkdir(exist_ok=True)

salmon = pandas.read_csv(snakemake.input.salmon, sep="\t")
beers_in = pandas.read_csv(snakemake.input.beers_in, sep="\t")
beers_out = pandas.read_csv(snakemake.input.beers_out, sep="\t")

# Merge the two alleles and pre/post mRNA together
# BEERS transcript ids are {sample}_{trancriptId}_{allele} and then maybe additional comments (particularly '_pre_mRNA')
#beers_in['allele'] = beers_in.TranscriptID.apply(lambda x: x.split('_')[2])
beers_in['TranscriptID'] = beers_in.TranscriptID.apply(lambda x: x.split('_')[1])
beers_in = beers_in.groupby(['TranscriptID', 'sample', 'GC_bias', 'pos_3prime_bias']).sum().reset_index()

beers_out['TranscriptID'] = beers_out.TranscriptID.apply(lambda x: x.split('_')[1])
beers_out = beers_out.groupby(['TranscriptID', 'sample', 'GC_bias', 'pos_3prime_bias']).sum().reset_index()

# input molecule counts give the underlying true TPM, after normalizing
true_tpm = beers_in.set_index("TranscriptID").groupby(['sample', 'GC_bias', 'pos_3prime_bias']).apply(lambda x: 1e6*x['count'] / x['count'].sum()).reset_index()
true_tpm.rename(columns={"count": "true_tpm"}, inplace=True)

# output molecules from BEERS give the true fragment counts
true_counts = beers_out.copy()
true_counts.rename(columns={"count": "true_count"}, inplace=True)

# Join with the Salmon quantified values
data = pandas.merge(salmon, true_tpm, on=['sample', 'GC_bias', 'pos_3prime_bias', 'TranscriptID'])
data = pandas.merge(data, true_counts, on=['sample', 'GC_bias', 'pos_3prime_bias', 'TranscriptID'])


# Summary stats
# Spearman correlations
count_corr = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: scipy.stats.spearmanr(x.NumReads, x.true_count, nan_policy='omit')[0])
tpm_corr = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: scipy.stats.spearmanr(x.TPM, x.true_tpm, nan_policy='omit')[0])
# Absolute median deviation - log2 scale
count_abs_med_dev = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: numpy.log2((x['NumReads'] + 1) / (x['true_count']+1)).median())
tpm_abs_med_dev = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: numpy.log2((x['TPM'] + 1) / (x['true_tpm']+1)).median())

stats = pandas.DataFrame({
    "count_corr": count_corr,
    "tpm_corr": tpm_corr,
    "count_abs_med_dev": count_abs_med_dev,
    "tpm_abs_med_dev": tpm_abs_med_dev,
}).reset_index()
stats.to_csv(out_dir / "summary_stats.txt", sep="\t", index=None)

# Plot summary stats
to_int_map = {"none": 0, "med": 1, "high": 2}
stats['GC_bias_int'] = stats['GC_bias'].map(to_int_map)
fig = sns.relplot(
    data=stats,
    x="GC_bias",
    y="count_corr",
    hue="GC_correct", 
    style="Pos_correct",
    col="sample",
    row="pos_3prime_bias",
    kind="line",
)
fig.savefig(out_dir / "corr.png", dpi=300)

# Generate figures
import pylab
import math
for measure in ['tpm', 'count']:
    truth_var = {"count": "true_count", "tpm": "true_tpm"}[measure]
    quant_var = {"count": "NumReads", "tpm": "TPM"}[measure]
    for sample in data['sample'].unique(): 
        by_study = data[data['sample'] == sample].set_index(["GC_bias", "pos_3prime_bias", "GC_correct", "Pos_correct"])
        conditions = by_study.index.unique()
        N_COLS = 4
        N_ROWS = math.ceil(len(conditions) / N_COLS)
        fig, axes = pylab.subplots(nrows=N_ROWS, ncols=N_COLS, figsize=(20,20), sharex=True, sharey=True)
        for condition, ax in zip(conditions, axes.flatten()):
            cond_data = by_study.loc[condition]
            ax.loglog(cond_data[truth_var], cond_data[quant_var], color="k", linestyle='', marker='o')
            #ax.set_aspect(1.0)
            ax.set_title(condition)
        #fig.tight_layout()
        fig.savefig(f"{snakemake.output.dir}/{measure}.scatter.sample{sample}.png", dpi=150)
