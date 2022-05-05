import pathlib
import re
import itertools

import pandas
import numpy

data_dir = pathlib.Path(snakemake.params.beers_output_data_folder)
print(data_dir)

transcript_id_re = re.compile(r"_(ENSMUST[0-9]+)_")

results = []
for sample in snakemake.params.samples:
    fastq_files = list(data_dir.glob(f"S{sample}_*.fastq"))
    print(f"Processing {fastq_files}")
    for fastq_file in fastq_files:
        print("Processing", fastq_file)
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
                results.append({
                    "GC_bias": snakemake.wildcards.GC_bias,
                    "pos_3prime_bias": snakemake.wildcards.pos_3prime_bias,
                    "sample": sample,
                    "TranscriptID": transcript_id, 
                    "GC_content": GC_content,
                })
GC_content = pandas.DataFrame(results)
GC_by_transcript = GC_content.groupby(["GC_bias", "pos_3prime_bias", "sample", "TranscriptID"]).GC_content.mean()
GC_by_transcript.reset_index().to_csv(snakemake.output.gc_content, sep="\t")
