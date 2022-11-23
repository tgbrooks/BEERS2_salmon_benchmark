import pathlib
import json
import numpy

run_configs = snakemake.config['run_configs']
lanes_used = snakemake.config['lanes_used']
sample_ids = snakemake.config['sample_ids']

pathlib.Path('results/seq_bias/').mkdir(exist_ok=True)
biases = ['none', 'med', 'high']

#TODO: just the first lane? If we use multiple lanes, they won't all be represented here
sam_paths = {(run): f"data/{run}/beers/results/S1_L1.sam"
                    for run in run_configs.keys()}

BASES = ['A', 'C', 'G', 'T']
COMPLEMENT_BASES = ['T', 'G', 'C', 'A']
READ_LENGTH = 100
results = []
for run, sam_path in sam_paths.items():
    fwd_count_matrix = numpy.zeros((4, READ_LENGTH))
    rev_count_matrix = numpy.zeros((4, READ_LENGTH))
    num_molecules = 0
    with open(sam_path) as sam:
        for line in sam:
            if line.startswith("@"):
                continue # Skip metadata
            id, flag, chrom, pos, _, cigar, _, _, _, seq, qual = line.split("\t")
            # SAM files store reads in reverse complement if they align to the negative strand
            negative_strand = int(flag) & 0x10 != 0
            rev_read = int(flag) & 0x80 != 0
            count_matrix = rev_count_matrix if rev_read else fwd_count_matrix
            if negative_strand:
                # Negative strand - reverse complement
                count_matrix += [[1 if x == char else 0 for x in seq[::-1]] for char in COMPLEMENT_BASES]
            else:
                # positive strand - no reverse complement
                count_matrix += [[1 if x == char else 0 for x in seq] for char in BASES]
            if not rev_read:
                num_molecules += 1
    results.append({ 
            "run": run,
            "GC_bias": run_configs[run]['GC_bias'],
            "pos_3prime_bias": run_configs[run]['pos_3prime_bias'],
            "primer_bias": run_configs[run]['primer_bias'],
            "sample": 1, # TODO: multiple samples?
            # Frequencies of bases by position for the forward read
            "fwd_frequencies": {base: list(fwd_count_matrix[i] / num_molecules) for i, base in enumerate(BASES)},
            # Frequencies of bases by position for the reverse read
            "rev_frequencies": {base: list(rev_count_matrix[i] / num_molecules) for i, base in enumerate(BASES)},
    })
json.dump(results, open(snakemake.output.seq_frequencies, 'w'))
