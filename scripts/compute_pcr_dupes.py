import collections
import pathlib
import re
import numpy
import pandas

run_configs = snakemake.config['run_configs']

pathlib.Path('results/polya_tail/').mkdir(exist_ok=True)

beersdir = pathlib.Path(snakemake.input.flag_file).parent
sam_path = beersdir / "results/S1_L1.sam"

# Read ids look like:
# BEERS:1:1:1:1101:11211:5592:INPUT_ID.72.3243.cdna.cdna.10	
# and we can extract the fragment ID from that as everything up to the 'cdna.cdna' section
fragment_id_re = re.compile(r"([a-zA-Z0-9_-][\.\d+]+).cdna.cdna")

BASES = ['A', 'C', 'G', 'T']
COMPLEMENT_BASES = ['T', 'G', 'C', 'A']
READ_LENGTH = 100
id_counts = {}
num_molecules = 0
with open(sam_path) as sam:
    for line in sam:
        if line.startswith("@"):
            continue # Skip metadata
        id, flag, chrom, pos, _, cigar, _, _, _, seq, qual = line.split("\t")
        # SAM files store reads in reverse complement if they align to the negative strand
        negative_strand = int(flag) & 0x10 != 0
        rev_read = int(flag) & 0x80 != 0
        if not rev_read:
            # Check just one read - they're separate lines in the SAM file
            # since that is at the 3' end of the molecule
            match = fragment_id_re.search(id)
            if match:
                frag_id = match.groups()[0]
                id_counts[frag_id] = id_counts.get(frag_id, 0) + 1
            else:
                print(f"Failed to identify fragment on line:\n{line}\nbreaking!")
                break
            num_molecules += 1
dupe_counts = collections.Counter(id_counts.values())
results = pandas.DataFrame({ 
        "num_PCR_dupes": list(dupe_counts.keys()),
        "count": list(dupe_counts.values()),
}).sort_values(by=["num_PCR_dupes"])
results.to_csv(snakemake.output.pcr_dupes, sep="\t", index=False)
