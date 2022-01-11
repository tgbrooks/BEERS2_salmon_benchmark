# BEERS2_salmon_benchmark
Benchmark of Salmon bias correction using BEERS2

BEERS2 simulates RNA-seq at the molecule level throughout library preparation and sequencing.
This allows the inclusion of realistic biases in the generated data, such as a GC-content bias.
Salmon is an RNA-seq pseudo-aligner that contains options for correcting for certain biases, in particular GC content bias, positional bias (like 3' bias), and sequence bias of start positions.
We therefore use BEERS2 to benchmark the performance of these bias-correction methods to assess how them may affect transcript quantification.

Pipeline written with Snakemake, run `run_pipeline.sh` to execute the pipeline on an LSF cluster.
