import pandas
import seaborn as sns

pcr_dupe_counts = []
for run, path in zip(snakemake.params.runs, snakemake.input.pcr_dupes):
    df = pandas.read_csv(path, sep="\t")
    df['run'] = run
    pcr_dupe_counts.append(df)
pcr_dupe_counts = pandas.concat(pcr_dupe_counts)

# Transform so that uniques are 0 duplicates rather than 1 copy
pcr_dupe_counts['num_PCR_dupes'] = pcr_dupe_counts['num_PCR_dupes'] - 1

fig = sns.catplot(
        x = "num_PCR_dupes",
        y = "count",
        col = "run",
        col_wrap = 4,
        data = pcr_dupe_counts[pcr_dupe_counts['num_PCR_dupes'] < 5],
        kind = "bar",
        height = 2.5,
        sharey=False,

)
fig.set(xlabel="PCR duplicates")
fig.savefig(snakemake.output.pcr_dupes, dpi=300)
