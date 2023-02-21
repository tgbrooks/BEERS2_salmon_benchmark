import pathlib

import numpy
import scipy.stats
import pandas


#lanes_used = [1]
#sample_ids = list(range(1,9))
lanes_used = snakemake.config['lanes_used']
sample_ids = snakemake.config['sample_ids']

bias_order = ["none", "med", "high"]
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
MIN_NUM_READS = 100

out_dir = pathlib.Path(snakemake.output.dir)
#out_dir = pathlib.Path("results/compare_accuracy")
out_dir.mkdir(exist_ok=True)

input_files = snakemake.input
#input_files = dict(
#        salmon = "results/salmon_quants.parquet",
#        beers_in = "results/beers_input_quants.parquet",
#        beers_out = "results/beers_output_quants.parquet",
#)
salmon = pandas.read_parquet(input_files['salmon'])
beers_in = pandas.read_parquet(input_files['beers_in'])
beers_out = pandas.read_parquet(input_files['beers_out'])

# Merge the two alleles and pre/post mRNA together
# BEERS transcript ids are {sample}_{trancriptId}_{allele} and then maybe additional comments (particularly '_pre_mRNA')
#beers_in['allele'] = beers_in.TranscriptID.apply(lambda x: x.split('_')[2])
beers_in = beers_in.groupby(['TranscriptID', 'GeneID', 'run', 'sample']).sum(numeric_only=True).reset_index()
beers_in.rename(columns={"count": "BEERS_input_count"}, inplace=True)

beers_out = beers_out.groupby(['TranscriptID', 'GeneID', 'run', 'sample']).sum(numeric_only=True).reset_index()

# input molecule counts give the underlying true TPM, after normalizing
true_tpm = beers_in.set_index(["TranscriptID", "GeneID"]).groupby(['run', 'sample']).apply(lambda x: 1e6*x['BEERS_input_count'] / x['BEERS_input_count'].sum()).reset_index()
true_tpm.rename(columns={"BEERS_input_count": "true_tpm"}, inplace=True)

# output molecules from BEERS give the true fragment counts
true_counts = beers_out.copy()
true_counts.rename(columns={"count": "true_count"}, inplace=True)

# Join with the Salmon quantified values
print(salmon.dtypes)
print(true_tpm.dtypes)
data = pandas.merge(salmon, true_tpm, on=['sample', 'run', 'TranscriptID', 'GeneID'])
data = pandas.merge(data, true_counts, on=['sample', 'run', 'TranscriptID', 'GeneID'])
data = pandas.merge(data, beers_in, on=['sample', 'run', 'TranscriptID', 'GeneID'])

# Compute CPM values
data['true_cpm'] = data.groupby(['sample', 'run'], group_keys=False)['true_count'].apply( lambda x: x / x.sum() * 1_000_000)
data['salmon_cpm'] = data.groupby(['sample','run'], group_keys=False)['NumReads'].apply( lambda x: x / x.sum() * 1_000_000)

print(data.head())
print(data.shape)

def correction_type(row):
    ''' Summarize all Salmon correction options in one string '''
    corrections = [correction for correction in ['GC', 'Pos', 'Seq'] if row[f'{correction}_correct']]
    if corrections:
        return ' '.join(corrections)
    return 'baseline'
data['correction type'] = data.apply(correction_type, axis=1)

# Sum all transcripts in a gene to get gene-level data
#gene_data = data.groupby(['sample', 'run', 'GeneID']).apply(lambda x: x[['TPM', 'salmon_cpm', 'NumReads', 'true_count', 'true_tpm', 'true_cpm']].sum(axis=0)).reset_index()

#### Summary stats ####

## Spearman correlations
count_corr = data.groupby(['sample', 'run', 'correction type']).apply(lambda x: scipy.stats.spearmanr(x.NumReads, x.true_count, nan_policy='omit')[0])
tpm_corr = data.groupby(['sample', 'run', 'correction type']).apply(lambda x: scipy.stats.spearmanr(x.TPM, x.true_tpm, nan_policy='omit')[0])

## Absolute median deviation - log2 scale
count_abs_med_dev = data.groupby(['sample', 'run', 'correction type']).apply(lambda x: numpy.abs(numpy.log2((x['NumReads'] + 1) / (x['true_count']+1))).median())
tpm_abs_med_dev = data.groupby(['sample', 'run', 'correction type']).apply(lambda x: numpy.abs(numpy.log2((x['TPM'] + 1) / (x['true_tpm']+1))).median())
cpm_abs_med_dev = data.groupby(['sample', 'run', 'correction type']).apply(lambda x: numpy.abs(numpy.log2((x['salmon_cpm'] + 1) / (x['true_cpm']+1))).median())

## Median absolute relative difference (MARD)
# Transcript-level MADR in both counts and TPMs compared to truth
data['ARD_count'] = (data['NumReads'] - data['true_count']).abs() / (data['NumReads'] + data['true_count'] + 1e-10) # Small epislon means that it goes to zero properly
data['ARD_cpm'] = (data['salmon_cpm'] - data['true_cpm']).abs() / (data['salmon_cpm'] + data['true_cpm'] + 1e-10)
data['ARD_tpm'] = (data['TPM'] - data['true_tpm']).abs() / (data['TPM'] + data['true_tpm'] + 1e-10)
MARD_count = data.groupby(['sample', 'run', 'correction type']).ARD_count.mean()
MARD_count.name = "MARD count"
MARD_tpm = data.groupby(['sample', 'run', 'correction type']).ARD_tpm.mean()
MARD_tpm.name = "MARD TPM"
MARD_cpm = data.groupby(['sample', 'run', 'correction type']).ARD_cpm.mean()
MARD_cpm.name = "MARD CPM"

# Gather all the summary stats together and save to disk
stats = pandas.DataFrame({
    "count_corr": count_corr,
    "tpm_corr": tpm_corr,
    "count_abs_med_dev": count_abs_med_dev,
    "tpm_abs_med_dev": tpm_abs_med_dev,
    "cpm_abs_med_dev": cpm_abs_med_dev,
    "MARD_count":  MARD_count,
    "MARD_cpm": MARD_cpm,
    "MARD_tpm": MARD_tpm,
})
stats_keys = stats.columns
stats = stats.reset_index()
stats.to_csv(snakemake.input.stats_file, sep="\t", index=None)
