import json

import pandas
import seaborn as sns

READ_LENGTH = 100
BASES = ['A', 'C', 'G', 'T']
biases = ['none', 'med', 'high']

json_data = json.load(open(snakemake.input.seq_frequencies))
frequencies = pandas.DataFrame([{
    "GC_bias": entry['GC_bias'],
    "pos_3prime_bias": entry['pos_3prime_bias'],
    "base": base,
    "position": read_pos,
    "read_direction": fwd_rev,
    "frequency": entry['fwd_frequencies' if fwd_rev == 'fwd' else 'rev_frequencies'][base][read_pos],
}  
    for entry in json_data
    for base in BASES
    for read_pos in range(READ_LENGTH)
    for fwd_rev in ['fwd', 'rev']]
)
print(frequencies.head())

fig = sns.relplot(
        x = 'position',
        y = 'frequency',
        hue = "base",
        col = "pos_3prime_bias",
        row = "GC_bias",
        data = frequencies.query("read_direction == 'fwd'"),
        kind = 'line',
        row_order = biases,
        col_order = biases,
        hue_order = BASES,
)
fig.savefig(snakemake.output.fwd_frequencies)

fig = sns.relplot(
        x = 'position',
        y = 'frequency',
        hue = "base",
        col = "pos_3prime_bias",
        row = "GC_bias",
        data = frequencies.query("read_direction == 'rev'"),
        kind = 'line',
        row_order = biases,
        col_order = biases,
        hue_order = BASES,
)
fig.savefig(snakemake.output.rev_frequencies)
