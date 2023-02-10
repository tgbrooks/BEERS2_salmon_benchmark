import pathlib
import re
import itertools
import gzip

import pandas
import numpy

import beers_utils.cigar

run_configs = snakemake.config['run_configs']
lanes_used = snakemake.config['lanes_used']
sample_ids = snakemake.config['sample_ids']

data_dir = pathlib.Path(snakemake.params.beers_output_data_folder)
print(data_dir)

# Load annotations
import re
exons = {}
strands = {}
with gzip.open('/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/Mus_musculus.GRCm38.93.gtf.gz') as annotations:
    transcript_id_re = re.compile(r'transcript_id "([a-zA-Z0-9]*)";')
    biotype_re = re.compile(r'gene_biotype "([a-zA-Z0-9_]*)";')
    for line in annotations:
        line = line.decode()
        if line.startswith("#"): continue
        chrom, source, type_, start, end, _, strand, *tail = line.split("\t")
        transcript_id = transcript_id_re.search(line)
        if type_ == 'transcript':
            transcript_id = transcript_id.groups()[0]
            exons[transcript_id] = []
            strands[transcript_id] = strand
        elif type_ == 'exon':
            transcript_id = transcript_id.groups()[0]
            exons[transcript_id].append((int(start), int(end)))
transcript_lens = {ID:sum(end - start for start,end in my_exons) for ID, my_exons in exons.items()}

transcript_id_re = re.compile(r"_(ENSMUST[0-9]+)_")

# Prepare position bins: we don't store at the base-level coverage
bins = numpy.linspace(0,1, 41)
bin_starts = bins[:-1]
bin_ends = bins[1:]

# Compute the coverage by position
results = []
for sample in snakemake.params.samples:
    sam_files = list(data_dir.glob(f"S{sample}_*.sam"))
    transcript_coverage = {}
    for sam_file in sam_files:
        print("Processing", sam_file)
        #counts_by_bin = numpy.zeros_like(bins, dtype=int)
        with open(sam_file) as sam:
            id = None
            i = 0
            for line in sam:
                i += 1
                if line.startswith("@"): continue # Comment at start
                transcript_id = transcript_id_re.search(line).groups()[0]
                read_name, _, chrome, start, _, cigar, *tail = line.split("\t")
                if transcript_id not in transcript_coverage:
                    transcript_coverage[transcript_id] = numpy.zeros_like(bin_starts)
                cov = transcript_coverage[transcript_id]

                my_exons = exons[transcript_id]
                transcript_len = transcript_lens[transcript_id]
                starts = [start for start,end in my_exons]
                ends = [end for start,end in my_exons]
                introns = [end - start for end, start in zip(ends[1:], starts[:-1])]
                start_exon = numpy.searchsorted(starts, start)

                # Compute start and end position relative to the transcript, normalized from 0 to 1
                # Note that it is assumed that the read is entirely contiguous within the transcript
                # after ignoring small deletes
                align_start = int(start) - starts[0] - sum(introns[:start_exon])
                cigar_split = beers_utils.cigar.split_cigar(cigar)
                align_end = align_start
                for op, length in cigar_split:
                    if op in "MD":
                        # Don't count 'N' in cigar as introns aren't part of the transcript
                        align_end += length
                align_start /= transcript_len
                align_end /= transcript_len

                # Mark the appropriate part of the transcript
                # Each bin is incremented proportional to the overlap of the read with the bin
                overlap =  numpy.clip(align_end, bin_starts, bin_ends) - numpy.clip(align_start, bin_starts, bin_ends)
                if strands[transcript_id] == '-':
                    overlap = overlap[::-1]
                cov += overlap

    cov_df = pandas.DataFrame.from_dict(
            transcript_coverage,
            orient = "index",
            columns = bin_starts,
    )
    cov_df.index.name = "transcript_id"
    cov_df['run'] = snakemake.wildcards.run
    cov_df["GC_bias"] = run_configs[snakemake.wildcards.run]['GC_bias']
    cov_df["pos_3prime_bias"] = run_configs[snakemake.wildcards.run]['pos_3prime_bias']
    cov_df["primer_bias"] = run_configs[snakemake.wildcards.run]['primer_bias']
    cov_df["sample"] = sample

    results.append(cov_df)
cov_results = pandas.concat(results, axis=0)
cov_results.to_csv(snakemake.output.pos_results, sep="\t")
