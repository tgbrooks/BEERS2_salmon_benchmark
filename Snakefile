import pathlib

sample_ids = [1,2]
GC_corrects = ['yes', 'no']
Pos_corrects = ['yes', 'no']

GC_bias = dict(
    none = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": 0.0,
    },
    med = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": -30,
    },
    high = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": -100,
    },
)

pos_3prime_bias = dict(
    none = {
        "breakpoint_prob_per_base": 0.0,
    },
    med = {
        "breakpoint_prob_per_base": 0.001,
    },
    high = {
        "breakpoint_prob_per_base": 0.005,
    },
)

wildcard_constraints:
    GC_bias = "[A-Za-z]+",
    pos_3prime_bias = "[A-Za-z]+",
    sample = "[0-9]+",
    GC_correct = "[A-Za-z]+",
    Pos_correct = "[A-Za-z]+",

rule all:
    input:
        salmon_quant = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}/quant.genes.sf",
                    GC_bias = GC_bias.keys(),
                    pos_3prime_bias=pos_3prime_bias.keys(),
                    sample = sample_ids,
                    GC_correct = GC_corrects,
                    Pos_correct = Pos_corrects),
        beers_input_quants = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/input_quant.txt",
                    GC_bias = GC_bias.keys(),
                    pos_3prime_bias=pos_3prime_bias.keys(),
                    sample = sample_ids),
        beers_output_quants = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/output_quant.txt",
                    GC_bias = GC_bias.keys(),
                    pos_3prime_bias=pos_3prime_bias.keys(),
                    sample = sample_ids),

rule run_beers:
    input:
        config_template = "config/config.template.json",
    output:
        config = "config/generated/config.GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}.json",
        beers = directory("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/"),
    run:
        #Generate config template
        import string
        config_template = string.Template(open(input.config_template, "r").read())
        config_fill_values = dict()
        config_fill_values.update(GC_bias[wildcards.GC_bias])
        config_fill_values.update(pos_3prime_bias[wildcards.pos_3prime_bias])
        config_fill_values['output_dir'] = output.beers +"/run"
        config = config_template.substitute(config_fill_values)
        with open(output.config, "w") as f:
            f.write(config)
        # Run beers
        beers_cmd =  f"run_beers.py -c {output.config} -r 1 -d -m serial prep_and_sequence_pipeline"
        print(beers_cmd)
        shell(beers_cmd)

rule get_beers_quant_files:
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/",
    output:
        input_quants = expand("data/GC_bias={{GC_bias}}.pos_3prime_bias={{pos_3prime_bias}}/sample{sample}/input_quant.txt", sample=sample_ids),
        output_quants = expand("data/GC_bias={{GC_bias}}.pos_3prime_bias={{pos_3prime_bias}}/sample{sample}/output_quant.txt", sample=sample_ids),
    run:
        import pandas
        # Sum the quantification files provided by beers per packet and then demultiplex by sample
        input_quantfiles = pathlib.Path(input[0]).glob("run_run1/library_prep_pipeline/logs/Pipeline/*/library_prep_pipeline_molecule_pkt*.quant_file")
        output_quantfiles = pathlib.Path(input[0]).glob("run_run1/library_prep_pipeline/data/*/library_prep_pipeline_result_molecule_pkt*.quant_file")

        input_quants_list = []
        for input_file in input_quantfiles:
            # Load all the packet quants
            input_quant = pandas.read_csv(input_file, sep="\t", header=None, index_col=0)
            input_quants_list.append(input_quant)
        input_quants = pandas.concat(input_quants_list, axis=0)
        for sample, filepath in zip(sample_ids, output.input_quants):
            # De-multiplex the samples
            # Transcripts are identified with sample# at the start
            prefix = f"{sample}_"
            sample_input_quants = input_quants[input_quants.index.str.startswith(prefix)]
            sample_input_quants.to_csv(filepath, sep="\t", header=None)

        output_quants_list = []
        for output_file in output_quantfiles:
            # Load all the packet quants
            output_quant = pandas.read_csv(output_file, sep="\t", header=None, index_col=0)
            output_quants_list.append(output_quant)
        output_quants = pandas.concat(output_quants_list, axis=0)
        for sample, filepath in zip(sample_ids, output.output_quants):
            # De-multiplex the samples
            # Transcripts are per-sample and identified with sample# at the start
            prefix = f"{sample}_"
            sample_output_quants = output_quants[output_quants.index.str.startswith(prefix)]
            sample_output_quants.to_csv(filepath, sep="\t", header=None)

rule generate_salmon_index:
    output:
        directory("index/baby_mouse")
    threads: 16
    resources:
        mem_mb=70000,
    shell:
        "salmon index -p 16 -i {output} \
             -t resources/baby.genome.cDNA.fasta "
             # TODO: decoys? Benchmark decoys?

rule run_salmon:
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/run_run1/controller/data/",
        "index/baby_mouse"
    output:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}/quant.genes.sf",
    params:
        gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf",# TODO: okay to use the full genome gtf?,
        out_folder = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}/",
    threads: 6
    resources:
        mem_mb=25000,
    run:
        args = "-l A --softclip --softclipOverhangs --biasSpeedSamp 10 -p 6"
        if wildcards.GC_correct == 'yes':
            args += ' --gcBias'
        if wildcards.Pos_correct == 'yes':
            args += ' --posBias'
        #TODO: handle --seqBias option
        fastqdir=pathlib.Path(input[0])
        if len(list(fastqdir.glob("*_2.fastq")))>0:
            shell(f"salmon quant -i {input[1]} -g {params.gtf_file} {args} -1 {input[0]}/S{wildcards.sample}_*_R1.fastq -2 {input[0]}/S{wildcards.sample}_*_R2.fastq -o {params.out_folder}")
        else:
            shell(f"salmon quant -i {input[1]} -g {params.gtf_file} {args} -r {input[0]}/S{wildcards.sample}_*_R1.fastq -o {params.out_folder}")

