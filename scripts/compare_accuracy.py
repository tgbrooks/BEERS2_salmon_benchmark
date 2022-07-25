import pathlib

import numpy
import scipy.stats
import pandas
import seaborn as sns
import pylab

from run_configs import run_configs, sample_ids, lanes_used

bias_order = ["none", "med", "high"]
run_order = run_configs.keys()
MIN_NUM_READS = 100

out_dir = pathlib.Path(snakemake.output.dir)
out_dir.mkdir(exist_ok=True)

salmon = pandas.read_csv(snakemake.input.salmon, sep="\t")
beers_in = pandas.read_csv(snakemake.input.beers_in, sep="\t")
beers_out = pandas.read_csv(snakemake.input.beers_out, sep="\t")

# Merge the two alleles and pre/post mRNA together
# BEERS transcript ids are {sample}_{trancriptId}_{allele} and then maybe additional comments (particularly '_pre_mRNA')
#beers_in['allele'] = beers_in.TranscriptID.apply(lambda x: x.split('_')[2])
beers_in = beers_in.groupby(['TranscriptID', 'GeneID', 'run', 'sample', 'GC_bias', 'pos_3prime_bias', "primer_bias"]).sum().reset_index()
beers_in.rename(columns={"count": "BEERS_input_count"}, inplace=True)

beers_out = beers_out.groupby(['TranscriptID', 'GeneID', 'run', 'sample', 'GC_bias', 'pos_3prime_bias', "primer_bias"]).sum().reset_index()

# input molecule counts give the underlying true TPM, after normalizing
true_tpm = beers_in.set_index(["TranscriptID", "GeneID"]).groupby(['run', 'sample', 'GC_bias', 'pos_3prime_bias', "primer_bias"]).apply(lambda x: 1e6*x['BEERS_input_count'] / x['BEERS_input_count'].sum()).reset_index()
true_tpm.rename(columns={"BEERS_input_count": "true_tpm"}, inplace=True)

# output molecules from BEERS give the true fragment counts
true_counts = beers_out.copy()
true_counts.rename(columns={"count": "true_count"}, inplace=True)

# Join with the Salmon quantified values
data = pandas.merge(salmon, true_tpm, on=['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'TranscriptID', 'GeneID'])
data = pandas.merge(data, true_counts, on=['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'TranscriptID', 'GeneID'])
data = pandas.merge(data, beers_in, on=['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'TranscriptID', 'GeneID'])
print(data.head())

# Compute CPM values
data['true_cpm'] = data.groupby(['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct'])['true_count'].apply( lambda x: x / x.sum() * 1_000_000)
data['salmon_cpm'] = data.groupby(['sample','run',  'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct'])['NumReads'].apply( lambda x: x / x.sum() * 1_000_000)

# Sum all transcripts in a gene to get gene-level data
gene_data = data.groupby(['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct', 'GeneID']).apply(lambda x: x[['TPM', 'salmon_cpm', 'NumReads', 'true_count', 'true_tpm', 'true_cpm']].sum(axis=0)).reset_index()

correction_type_order = ['baseline', 'GC', 'Pos', 'GC Pos', 'Seq', 'GC Seq', 'Pos Seq', 'GC Pos Seq']
def correction_type(row):
    ''' Summarize all Salmon correction options in one string '''
    corrections = [correction for correction in ['GC', 'Pos', 'Seq'] if row[f'{correction}_correct']]
    if corrections:
        return ' '.join(corrections)
    return 'baseline'

## Summary stats
## Spearman correlations
count_corr = data.groupby(['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct']).apply(lambda x: scipy.stats.spearmanr(x.NumReads, x.true_count, nan_policy='omit')[0])
tpm_corr = data.groupby(['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct']).apply(lambda x: scipy.stats.spearmanr(x.TPM, x.true_tpm, nan_policy='omit')[0])
# Absolute median deviation - log2 scale
count_abs_med_dev = data.groupby(['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct']).apply(lambda x: numpy.abs(numpy.log2((x['NumReads'] + 1) / (x['true_count']+1))).median())
tpm_abs_med_dev = data.groupby(['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct']).apply(lambda x: numpy.abs(numpy.log2((x['TPM'] + 1) / (x['true_tpm']+1))).median())
cpm_abs_med_dev = data.groupby(['sample', 'run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'GC_correct', 'Pos_correct', 'Seq_correct']).apply(lambda x: numpy.abs(numpy.log2((x['salmon_cpm'] + 1) / (x['true_cpm']+1))).median())

stats = pandas.DataFrame({
    "count_corr": count_corr,
    "tpm_corr": tpm_corr,
    "count_abs_med_dev": count_abs_med_dev,
    "tpm_abs_med_dev": tpm_abs_med_dev,
    "cpm_abs_med_dev": cpm_abs_med_dev,
}).reset_index()
stats.to_csv(out_dir / "summary_stats.txt", sep="\t", index=None)
stats['bias_correct'] = stats['GC_correct'].map({True: "GC ", False: ""}) + stats['Pos_correct'].map({True: "Pos ", False: ""}) + stats["Seq_correct"].map({True: "Seq", False: ""})

# Plot correlation summary stats
fig = sns.catplot(
    data=stats.query("sample == 1"),
    x="bias_correct",
    y="count_corr",
    col="run",
    col_wrap=3,
    col_order=run_order,
    kind="bar",
)
fig.savefig(out_dir / "count_corr.png", dpi=300)

# Plot TPM correlation summary stats
fig = sns.catplot(
    data=stats.query("sample == 1"),
    x="bias_correct",
    y="tpm_corr",
    col="run",
    col_wrap=3,
    col_order=run_order,
    kind="bar",
)
fig.savefig(out_dir / "tpm_corr.png", dpi=300)

# Plot abs med dev summary stats
fig = sns.catplot(
    data=stats.query("sample == 1"),
    x="bias_correct",
    y="count_abs_med_dev",
    col="run",
    col_wrap=3,
    col_order=run_order,
    kind="bar",
)
[ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right') for ax in fig.axes.flatten()]
fig.savefig(out_dir / "count_abs_med_dev.png", dpi=300)

# Plot TPM abs med dev summary stats
fig = sns.catplot(
    data=stats.query("sample == 1"),
    x="bias_correct",
    y="tpm_abs_med_dev",
    col="run",
    col_wrap=3,
    col_order=run_order,
    kind="bar",
)
[ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right') for ax in fig.axes.flatten()]
fig.savefig(out_dir / "tpm_abs_med_dev.png", dpi=300)

# Plot CPM abs med dev summary stats
fig = sns.catplot(
    data=stats.query("sample == 1"),
    x="bias_correct",
    y="cpm_abs_med_dev",
    col="run",
    col_wrap=3,
    col_order = run_order,
    kind="bar",
)
[ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right') for ax in fig.axes.flatten()]
fig.savefig(out_dir / "cpm_abs_med_dev.png", dpi=300)

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
fig.figure.suptitle("Overall counts (in no-bias simulations)")
fig.tight_layout()
fig.savefig(out_dir / "margins.png", dpi=300)

# Plot the total number of simulated reads by the different BEERS bias amounts
fig = sns.catplot(
    data = data.query("GC_correct == False and Pos_correct == False and Seq_correct == False").groupby(['run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'sample']).sum().reset_index(),
    x = "run",
    y = "true_count",
    col = "sample",
    kind = "bar",
    order=run_order,
)
fig.figure.suptitle("BEERS output counts")
[ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right') for ax in fig.axes.flatten()]
fig.tight_layout()
fig.savefig(out_dir / "margins.by_bias_amount.png", dpi=300)

# Plot the total number of quantified reads by the different BEERS bias amounts
fig = sns.catplot(
    data = data.query("GC_correct == False and Pos_correct == False and Seq_correct == False").groupby(['run', 'GC_bias', 'pos_3prime_bias', 'primer_bias', 'sample']).sum().reset_index(),
    x = "run",
    y = "NumReads",
    col = "sample",
    kind = "bar",
    order=run_order,
)
fig.figure.suptitle("Salmon total NumReads (without bias corrections)")
[ax.set_xticklabels(ax.get_xticklabels(), rotation=45, horizontalalignment='right') for ax in fig.axes.flatten()]
fig.tight_layout()
fig.savefig(out_dir / "margins.salmon.by_bias_amount.png", dpi=300)



# Compare the salmon runs within each gene
# Baseline error will be the no-correction runs
# Compare each to that
# TPM Transcript-level:
d = data.copy()
d['abs_error'] = d['TPM'] - d['true_tpm']
baseline = d.query("GC_correct == False and Pos_correct == False and Seq_correct == False")
beers_run = ['sample', 'run', 'TranscriptID', 'GeneID']
d = pandas.merge(d, baseline[beers_run +["abs_error"]], left_on = beers_run, right_on = beers_run, suffixes = ('', '_baseline'))
d['rel_abs_error'] = d['abs_error'].abs() / d['abs_error_baseline'].abs()
d['log_rel_abs_error'] = numpy.log2(d.rel_abs_error)
d['correction type'] = d.apply(correction_type, axis=1)
d = d.query("sample == 1 and (GC_correct != False  or Pos_correct != False or Seq_correct != False)")

fig = sns.displot(
    x = "log_rel_abs_error",
    data = d,
    hue = "correction type",
    col = "run",
    col_wrap = 3,
    kind = "kde",
    clip = (-2.5, 2.5),
    col_order = run_order,
    hue_order = correction_type_order,
    common_norm = False,
)
fig.refline(x=0)
fig.set_axis_labels("log2( |TPM err| / |TPM baseline error|)", "Transcript Density")
fig.savefig(out_dir / "abs_error.compared_to_baseline.png", dpi=300)

# Same as above but gene-level
d = gene_data.copy()
d['abs_error'] = d['TPM'] - d['true_tpm']
baseline = d.query("GC_correct == False and Pos_correct == False and Seq_correct == False")
beers_run = ['sample', 'run', 'GeneID']
d = pandas.merge(d, baseline[beers_run +["abs_error"]], left_on = beers_run, right_on = beers_run, suffixes = ('', '_baseline'))
d['rel_abs_error'] = d['abs_error'].abs() / d['abs_error_baseline'].abs()
d['log_rel_abs_error'] = numpy.log2(d.rel_abs_error)
d['correction type'] = d.apply(correction_type, axis=1)
d = d.query("sample == 1 and (GC_correct != False  or Pos_correct != False or Seq_correct != False)")

fig = sns.displot(
    x = "log_rel_abs_error",
    data = d,
    hue = "correction type",
    col = "run",
    col_wrap = 3,
    kind = "kde",
    clip = (-2.5, 2.5),
    col_order = run_order,
    hue_order = correction_type_order,
    common_norm = False,
)
fig.refline(x=0)
fig.set_axis_labels("log2( |TPM err| / |TPM baseline error|)", "Gene Density")
print("TPM Frac in range:", (d['log_rel_abs_error'].abs() < 2.5).mean())
fig.savefig(out_dir / "abs_error.compared_to_baseline.by_gene.png", dpi=300)

# CPM Transcript-level:
d = data.copy()
d['abs_error'] = d['salmon_cpm'] - d['true_cpm']
baseline = d.query("GC_correct == False and Pos_correct == False and Seq_correct == False")
beers_run = ['sample', 'run', 'TranscriptID', 'GeneID']
d = pandas.merge(d, baseline[beers_run +["abs_error"]], left_on = beers_run, right_on = beers_run, suffixes = ('', '_baseline'))
d['rel_abs_error'] = d['abs_error'].abs() / d['abs_error_baseline'].abs()
d['log_rel_abs_error'] = numpy.log2(d.rel_abs_error)
d['correction type'] = d.apply(correction_type, axis=1)
d = d.query("sample == 1 and (GC_correct != False  or Pos_correct != False or Seq_correct != False)")

fig = sns.displot(
    x = "log_rel_abs_error",
    data = d,
    hue = "correction type",
    col = "run",
    col_wrap = 3,
    kind = "kde",
    clip = (-2.5, 2.5),
    col_order = run_order,
    hue_order = correction_type_order,
    common_norm = False,
)
fig.refline(x=0)
fig.set_axis_labels("log2( |CPM err| / |CPM baseline error|)", "Transcript Density")
fig.savefig(out_dir / "abs_error.compared_to_baseline.CPM.png", dpi=300)

# Same as above but gene-level
d = gene_data.copy()
d['abs_error'] = d['salmon_cpm'] - d['true_tpm']
baseline = d.query("GC_correct == False and Pos_correct == False and Seq_correct == False")
beers_run = ['sample', 'run', 'GeneID']
d = pandas.merge(d, baseline[beers_run +["abs_error"]], left_on = beers_run, right_on = beers_run, suffixes = ('', '_baseline'))
d['rel_abs_error'] = d['abs_error'].abs() / d['abs_error_baseline'].abs()
d['log_rel_abs_error'] = numpy.log2(d.rel_abs_error)
d['correction type'] = d.apply(correction_type, axis=1)
d = d.query("sample == 1 and (GC_correct != False  or Pos_correct != False or Seq_correct != False)")

fig = sns.displot(
    x = "log_rel_abs_error",
    data = d,
    hue = "correction type",
    col = "run",
    col_wrap = 3,
    kind = "kde",
    clip = (-2.5, 2.5),
    col_order = run_order,
    hue_order = correction_type_order,
    common_norm = False,
)
fig.refline(x=0)
fig.set_axis_labels("log2( |CPM err| / |CPM baseline error|)", "Gene Density")
print("CPM Frac in range:", (d['log_rel_abs_error'].abs() < 2.5).mean())
fig.savefig(out_dir / "abs_error.compared_to_baseline.CPM.by_gene.png", dpi=300)


# Error broken down by bias-amount in the genes
# We can measure the effect of Salmon's bias-correcting procedures via comparing
# the EffectiveLength values across the different correction methods
# Genes with more effect from positional/GC bias correction will
# then have larger differences in EffectiveLength
# We break down the genes by the effects of the correction bias to see if
# genes with more bias-correction are quantified well
# Also just plot the bias factor distributions
# See https://github.com/COMBINE-lab/salmon/discussions/752

# For GC bias:
d = data.copy()
GC_baseline = d.query("pos_3prime_bias == 'none' and GC_correct == False and Pos_correct == False").set_index(['GC_bias', 'sample', 'TranscriptID'])
GC_corrected = d.query("pos_3prime_bias == 'none' and GC_correct == True and Pos_correct == False").set_index(['GC_bias', 'sample', 'TranscriptID'])
GC_correction = pandas.DataFrame({
    "GC_bias_factor":  GC_corrected.EffectiveLength / GC_baseline.EffectiveLength,
    "rel_abs_error": ((GC_corrected.TPM - GC_corrected.true_tpm) / (GC_baseline.TPM - GC_baseline.true_tpm)).abs(),
    "baseline_num_reads": GC_baseline.NumReads,
    "log10_baseline_num_reads": numpy.log10(GC_baseline.NumReads + 1),
}).reset_index()
fig = sns.relplot(
        x = "GC_bias_factor",
        y = "rel_abs_error",
        size = "log10_baseline_num_reads",
        data = GC_correction.query(f"baseline_num_reads > {MIN_NUM_READS}"),
        col = "sample",
        row = "GC_bias",
        row_order=bias_order,
        kind = "scatter",
        #kind = "kde",
        #common_norm = False,
)
fig.set(ylim=(0.0, 2.0))
fig.refline(y=1)
fig.savefig(out_dir / "rel_abs_err.by_GC_bias_factor.png", dpi=300)
fig = sns.displot(
        x = "GC_bias_factor",
        data = GC_correction.query(f"baseline_num_reads > {MIN_NUM_READS}"),
        col = "sample",
        hue = "GC_bias",
        hue_order=bias_order,
        kind = "kde",
        common_norm = False,
)
fig.savefig(out_dir / "GC_bias_factor.dist.png", dpi=300)

# Same for Positional bias
pos_baseline = d.query("GC_bias == 'none' and GC_correct == False and Pos_correct == False").set_index(['pos_3prime_bias', 'sample', 'TranscriptID'])
pos_corrected = d.query("GC_bias == 'none' and GC_correct == False and Pos_correct == True").set_index(['pos_3prime_bias', 'sample', 'TranscriptID'])
pos_correction = pandas.DataFrame({
    "pos_bias_factor":  pos_corrected.EffectiveLength / pos_baseline.EffectiveLength,
    "rel_abs_error": ((pos_corrected.TPM - pos_corrected.true_tpm) / (pos_baseline.TPM - pos_baseline.true_tpm)).abs(),
    "baseline_num_reads": pos_baseline.NumReads,
    "log10_baseline_num_reads": numpy.log10(pos_baseline.NumReads + 1),
}).reset_index()
fig = sns.relplot(
        x = "pos_bias_factor",
        y = "rel_abs_error",
        size = "log10_baseline_num_reads",
        data = pos_correction.query(f"baseline_num_reads > {MIN_NUM_READS}"),
        col = "sample",
        row = "pos_3prime_bias",
        row_order=bias_order,
        kind = "scatter",
        #kind = "kde",
        #common_norm = False,
)
fig.set(ylim=(0.0, 2.0))
fig.refline(y=1)
fig.savefig(out_dir / "rel_abs_err.by_pos_bias_factor.png", dpi=300)
fig = sns.displot(
        x = "pos_bias_factor",
        data = pos_correction.query(f"baseline_num_reads > {MIN_NUM_READS}"),
        col = "sample",
        hue = "pos_3prime_bias",
        hue_order=bias_order,
        kind = "kde",
        common_norm = False,
)
fig.savefig(out_dir / "pos_bias_factor.dist.png", dpi=300)




## Generate scatterplot figures
# NOTE: these aren't so useful for the full genome
#import pylab
#import math
#for measure in ['tpm', 'count']:
#    truth_var = {"count": "true_count", "tpm": "true_tpm"}[measure]
#    quant_var = {"count": "NumReads", "tpm": "TPM"}[measure]
#    for sample in data['sample'].unique():
#        by_study = data[data['sample'] == sample].set_index(["GC_bias", "pos_3prime_bias", "GC_correct", "Pos_correct"])
#        conditions = by_study.index.unique()
#        N_COLS = 4
#        N_ROWS = math.ceil(len(conditions) / N_COLS)
#        fig, axes = pylab.subplots(nrows=N_ROWS, ncols=N_COLS, figsize=(20,20), sharex=True, sharey=True)
#        for condition, ax in zip(conditions, axes.flatten()):
#            cond_data = by_study.loc[condition]
#            ax.loglog(cond_data[truth_var], cond_data[quant_var], color="k", linestyle='', marker='o')
#            #ax.set_aspect(1.0)
#            ax.set_title(condition)
#        for ax in axes[-1,:]:
#            ax.set_xlabel("True")
#        for ax in axes[:,0]:
#            ax.set_ylabel("Est.")
#        #fig.tight_layout()
#        fig.savefig(f"{snakemake.output.dir}/{measure}.scatter.sample{sample}.png", dpi=150)
