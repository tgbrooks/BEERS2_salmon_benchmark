import pathlib

import numpy
import scipy.stats
import pandas
import seaborn as sns
import pylab

bias_order = ["none", "med", "high"]
correction_type_order = ['GC ', 'Pos', 'GC Pos']

out_dir = pathlib.Path(snakemake.output.dir)
out_dir.mkdir(exist_ok=True)

salmon = pandas.read_csv(snakemake.input.salmon, sep="\t")
beers_in = pandas.read_csv(snakemake.input.beers_in, sep="\t")
beers_out = pandas.read_csv(snakemake.input.beers_out, sep="\t")

# Merge the two alleles and pre/post mRNA together
# BEERS transcript ids are {sample}_{trancriptId}_{allele} and then maybe additional comments (particularly '_pre_mRNA')
#beers_in['allele'] = beers_in.TranscriptID.apply(lambda x: x.split('_')[2])
beers_in = beers_in.groupby(['TranscriptID', 'GeneID', 'sample', 'GC_bias', 'pos_3prime_bias']).sum().reset_index()
beers_in.rename(columns={"count": "BEERS_input_count"}, inplace=True)

beers_out = beers_out.groupby(['TranscriptID', 'GeneID', 'sample', 'GC_bias', 'pos_3prime_bias']).sum().reset_index()

# input molecule counts give the underlying true TPM, after normalizing
true_tpm = beers_in.set_index(["TranscriptID", "GeneID"]).groupby(['sample', 'GC_bias', 'pos_3prime_bias']).apply(lambda x: 1e6*x['BEERS_input_count'] / x['BEERS_input_count'].sum()).reset_index()
true_tpm.rename(columns={"BEERS_input_count": "true_tpm"}, inplace=True)

# output molecules from BEERS give the true fragment counts
true_counts = beers_out.copy()
true_counts.rename(columns={"count": "true_count"}, inplace=True)

# Join with the Salmon quantified values
data = pandas.merge(salmon, true_tpm, on=['sample', 'GC_bias', 'pos_3prime_bias', 'TranscriptID', 'GeneID'])
data = pandas.merge(data, true_counts, on=['sample', 'GC_bias', 'pos_3prime_bias', 'TranscriptID', 'GeneID'])
data = pandas.merge(data, beers_in, on=['sample', 'GC_bias', 'pos_3prime_bias', 'TranscriptID', 'GeneID'])

# Compute CPM values
data['true_cpm'] = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct'])['true_count'].apply( lambda x: x / x.sum() * 1_000_000)
data['salmon_cpm'] = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct'])['NumReads'].apply( lambda x: x / x.sum() * 1_000_000)

# Sum all transcripts in a gene to get gene-level data
gene_data = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct', 'GeneID']).apply(lambda x: x[['TPM', 'NumReads', 'true_count', 'true_tpm']].sum(axis=0)).reset_index()


# Summary stats
# Spearman correlations
count_corr = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: scipy.stats.spearmanr(x.NumReads, x.true_count, nan_policy='omit')[0])
tpm_corr = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: scipy.stats.spearmanr(x.TPM, x.true_tpm, nan_policy='omit')[0])
# Absolute median deviation - log2 scale
count_abs_med_dev = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: numpy.abs(numpy.log2((x['NumReads'] + 1) / (x['true_count']+1))).median())
tpm_abs_med_dev = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: numpy.abs(numpy.log2((x['TPM'] + 1) / (x['true_tpm']+1))).median())
cpm_abs_med_dev = data.groupby(['sample', 'GC_bias', 'pos_3prime_bias', 'GC_correct', 'Pos_correct']).apply(lambda x: numpy.abs(numpy.log2((x['salmon_cpm'] + 1) / (x['true_cpm']+1))).median())

stats = pandas.DataFrame({
    "count_corr": count_corr,
    "tpm_corr": tpm_corr,
    "count_abs_med_dev": count_abs_med_dev,
    "tpm_abs_med_dev": tpm_abs_med_dev,
    "cpm_abs_med_dev": cpm_abs_med_dev,
}).reset_index()
stats.to_csv(out_dir / "summary_stats.txt", sep="\t", index=None)

# Plot correlation summary stats
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
fig.savefig(out_dir / "count_corr.png", dpi=300)

# Plot TPM correlation summary stats
fig = sns.relplot(
    data=stats,
    x="GC_bias",
    y="tpm_corr",
    hue="GC_correct",
    style="Pos_correct",
    col="sample",
    row="pos_3prime_bias",
    kind="line",
)
fig.savefig(out_dir / "tpm_corr.png", dpi=300)

# Plot abs med dev summary stats
fig = sns.relplot(
    data=stats,
    x="GC_bias",
    y="count_abs_med_dev",
    hue="GC_correct",
    style="Pos_correct",
    col="sample",
    row="pos_3prime_bias",
    kind="line",
)
fig.savefig(out_dir / "count_abs_med_dev.png", dpi=300)

# Plot TPM abs med dev summary stats
fig = sns.relplot(
    data=stats,
    x="GC_bias",
    y="tpm_abs_med_dev",
    hue="GC_correct",
    style="Pos_correct",
    col="sample",
    row="pos_3prime_bias",
    kind="line",
)
fig.savefig(out_dir / "tpm_abs_med_dev.png", dpi=300)

# Plot CPM abs med dev summary stats
fig = sns.relplot(
    data=stats,
    x="GC_bias",
    y="cpm_abs_med_dev",
    hue="GC_correct",
    style="Pos_correct",
    col="sample",
    row="pos_3prime_bias",
    kind="line",
)
fig.savefig(out_dir / "cpm_abs_med_dev.png", dpi=300)


#  Plot in the no-bias no-correction case the effect of length on counts
d = data.query("GC_bias == 'none' and pos_3prime_bias == 'none' and GC_correct == False and Pos_correct == False")
fig = sns.relplot(
    data = d.melt(
        id_vars = ["sample", "Length"],
        value_vars = ["NumReads", "true_count", "BEERS_input_count"],
        var_name = "count_type",
        value_name = "count",
    ).reset_index().query("count > 0"),
    x = "Length",
    y = "count",
    hue = "count_type",
    col = "sample",
)
fig.set(yscale="log")
fig.savefig(out_dir /  "no_bias.length_effect.png", dpi=300)

# Plot the total number of reads by the different methods (i.e. the margins)
margins = data.melt(
        id_vars = ["sample", "GC_bias", "pos_3prime_bias", "GC_correct", "Pos_correct"],
        value_vars = ["NumReads", "true_count", "BEERS_input_count"],
        var_name = "count_type",
        value_name = "count"
).reset_index().groupby( ["sample", "GC_bias", "pos_3prime_bias", "GC_correct", "Pos_correct", "count_type"])['count'].sum().reset_index()
fig = sns.catplot(
    data = margins.query("GC_bias == 'none' and pos_3prime_bias == 'none'"),
    x = "GC_correct",
    y = "count",
    hue = "Pos_correct",
    row = "sample",
    col = "count_type",
    kind = "bar",
)
fig.savefig(out_dir / "margins.png", dpi=300)


# Compare the 4 salmon runs within each gene
# Baseline error will be the no-correction runs
# Compare each to that
d = data.copy()
d['abs_error'] = d['TPM'] - d['true_tpm']
baseline = d.query("GC_correct == False and Pos_correct == False")
beers_run = ['sample', 'GC_bias', 'pos_3prime_bias', "TranscriptID", "GeneID"]
d = pandas.merge(d, baseline[beers_run +["abs_error"]], left_on = beers_run, right_on = beers_run, suffixes = ('', '_baseline'))
d['rel_abs_error'] = d['abs_error'].abs() / d['abs_error_baseline'].abs()
d['log_rel_abs_error'] = numpy.log2(d.rel_abs_error)
d['correction type'] = d.apply(lambda x: ('GC ' if x['GC_correct'] else '') + ('Pos' if x['Pos_correct'] else ''), axis=1)
d = d.query("sample == 1 and (GC_correct != False  or Pos_correct != False)")

fig = sns.displot(
    x = "log_rel_abs_error",
    data = d,
    hue = "correction type",
    row = "GC_bias",
    col = "pos_3prime_bias",
    kind = "kde",
    clip = (-2.5, 2.5),
    row_order = bias_order,
    col_order = bias_order,
    hue_order = correction_type_order,
    common_norm = False,
)
fig.refline(x=0)
fig.set_axis_labels("log2( |TPM err| / |TPM baseline error|)", "Transcript Density")
fig.savefig(out_dir / "abs_error.compared_to_baseline.png", dpi=300)

# Same as above but gene-level
d = gene_data.copy()
d['abs_error'] = d['TPM'] - d['true_tpm']
baseline = d.query("GC_correct == False and Pos_correct == False")
beers_run = ['sample', 'GC_bias', 'pos_3prime_bias', "GeneID"]
d = pandas.merge(d, baseline[beers_run +["abs_error"]], left_on = beers_run, right_on = beers_run, suffixes = ('', '_baseline'))
d['rel_abs_error'] = d['abs_error'].abs() / d['abs_error_baseline'].abs()
d['log_rel_abs_error'] = numpy.log2(d.rel_abs_error)
d['correction type'] = d.apply(lambda x: ('GC ' if x['GC_correct'] else '') + ('Pos' if x['Pos_correct'] else ''), axis=1)
d = d.query("sample == 1 and (GC_correct != False  or Pos_correct != False)")

fig = sns.displot(
    x = "log_rel_abs_error",
    data = d,
    hue = "correction type",
    row = "GC_bias",
    col = "pos_3prime_bias",
    kind = "kde",
    clip = (-2.5, 2.5),
    row_order = bias_order,
    col_order = bias_order,
    hue_order = correction_type_order,
    common_norm = False,
)
fig.refline(x=0)
fig.set_axis_labels("log2( |TPM err| / |TPM baseline error|)", "Gene Density")
ax = fig.axes_dict[('none', 'none')]
print("Frac in range:", (d['log_rel_abs_error'].abs() < 2.5).mean())
fig.savefig(out_dir / "abs_error.compared_to_baseline.by_gene.png", dpi=300)


# Generate scatterplot figures
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
        for ax in axes[-1,:]:
            ax.set_xlabel("True")
        for ax in axes[:,0]:
            ax.set_ylabel("Est.")
        #fig.tight_layout()
        fig.savefig(f"{snakemake.output.dir}/{measure}.scatter.sample{sample}.png", dpi=150)
