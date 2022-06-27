import pathlib
import re
import itertools

import pandas
import numpy

from run_configs import run_configs

data_dir = pathlib.Path(snakemake.params.beers_output_data_folder)
print(data_dir)

transcript_id_re = re.compile(r"_(ENSMUST[0-9]+)_")

bins = numpy.linspace(0,1, 101)
bins[-1] = 1.001 # make it right-inclusive on the last bin

results = []
for sample in snakemake.params.samples:
    fastq_files = list(data_dir.glob(f"S{sample}_*.fastq"))
    print(f"Processing {fastq_files}")
    for fastq_file in fastq_files:
        print("Processing", fastq_file)
        counts_by_bin = numpy.zeros_like(bins, dtype=int)
        with open(fastq_file) as fastq:
            id = None
            i = 0
            while True:
                line = fastq.readline()
                if line == '':
                    break
                i += 1
                assert line.startswith("@"), f"Error in {fastq_file} on line {i}: expected @ not found instead received {repr(line)}"
                transcript_id = transcript_id_re.search(line).groups()[0]
                seq = fastq.readline().strip()
                GC_content = (seq.count('C') + seq.count('G') + seq.count('N')/2) / len(seq)
                fastq.readline() # + -- skip
                fastq.readline() # quality -- skip
                bin_ = numpy.searchsorted(bins, GC_content)
                counts_by_bin[bin_] += 1
    results.append(pandas.DataFrame({
        "run": snakemake.wildcards.run,
        "GC_bias": run_configs[snakemake.wildcards.run]['GC_bias'],
        "pos_3prime_bias": run_configs[snakemake.wildcards.run]['pos_3prime_bias'],
        "primer_bias": run_configs[snakemake.wildcards.run]['primer_bias'],
        "sample": sample,
        "GC_content": bins,
        "count": counts_by_bin,
    }))
GC_content = pandas.concat(results, axis=0)
GC_content.to_csv(snakemake.output.gc_content, sep="\t", index=False)
