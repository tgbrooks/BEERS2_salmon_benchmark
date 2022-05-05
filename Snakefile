import pathlib
import itertools
import util

sample_ids = [1,2]
GC_corrects = ['yes', 'no']
Pos_corrects = ['yes', 'no']
lanes_used = [1]

GC_bias = dict(
    none = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": 0.0,
    },
    med = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": -10,
    },
    high = {
        "gc_bias_constant": 1.0,
        "gc_bias_linear": 0.0,
        "gc_bias_quadratic": -50,
    },
)

pos_3prime_bias = dict(
    none = {
        "breakpoint_prob_per_base": 0.0,
    },
    med = {
        "breakpoint_prob_per_base": 0.0002,
    },
    high = {
        "breakpoint_prob_per_base": 0.001,
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
        "results/accuracy/",
        "results/salmon_quants.txt",
        "results/gc_summary.txt",
        "results/gc_content/",
        "results/coverage/",
        "/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/baby_genome/run_1/CAMPAREE/data/",
        # For comparison, align just a few with STAR
        "data/GC_bias=none.pos_3prime_bias=none/sample1/STAR/",
        "data/GC_bias=high.pos_3prime_bias=high/sample1/STAR/",
        salmon_quant = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}/quant.sf",
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
        beers_output_bam = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.bam",
                    GC_bias = GC_bias.keys(),
                    pos_3prime_bias=pos_3prime_bias.keys(),
                    sample = sample_ids),
        #renamed_gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.renamed_chromes.gtf",

rule run_CAMPAREE:
    input:
        #'config/camparee_config.yaml',
        'config/baby_genome.camparee_config.yaml',
    output:
        directory('/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/baby_genome/run_1/CAMPAREE/data/'),
    params:
        '/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/baby_genome/run_1/'
    resources:
        mem_mb = 10_000,
    shell:
        'rm -r {params[0]} && run_camparee.py -r 1 -d -c {input[0]}'

rule run_beers:
    input:
        config_template = "config/config.template.json",
        camparee_output = "/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/baby_genome/run_1/CAMPAREE/data/",
        #reference_genome = "/project/itmatlab/for_tom/BEERS2_benchmark/resources/mm10/Mus_musculus.GRCm38.dna.oneline_seqs.fa", # Full genome
        reference_genome = "/home/thobr/BEERS2/CAMPAREE/resources/baby_genome.mm10/baby_genome.mm10.oneline_seqs.fa", # Baby genome
    output:
        flag_file = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/finished_flag"
    params:
        beers_dir = directory("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/"),
        config = "config/generated/config.GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}.json",
    resources:
        mem_mb = 18_000
    run:
        #Generate config template
        import string
        config_template = string.Template(open(input.config_template, "r").read())
        config_fill_values = dict()
        config_fill_values.update(GC_bias[wildcards.GC_bias])
        config_fill_values.update(pos_3prime_bias[wildcards.pos_3prime_bias])
        config_fill_values['camparee_output'] = input.camparee_output
        config_fill_values['output_dir'] = params.beers_dir +"/run"
        config_fill_values['reference_genome'] = input.reference_genome
        config = config_template.substitute(config_fill_values)
        with open(params.config, "w") as f:
            f.write(config)
        # Run beers
        beers_cmd =  f"run_beers.py -c {params.config} -r 1 -d -m serial prep_and_sequence_pipeline"
        print(beers_cmd)
        shell(beers_cmd)
        shell("touch {output.flag_file}")

rule get_beers_quant_files:
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/finished_flag",
    params:
        beers_dir = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/",
    output:
        input_quants = expand("data/GC_bias={{GC_bias}}.pos_3prime_bias={{pos_3prime_bias}}/sample{sample}/input_quant.txt", sample=sample_ids),
        output_quants = expand("data/GC_bias={{GC_bias}}.pos_3prime_bias={{pos_3prime_bias}}/sample{sample}/output_quant.txt", sample=sample_ids),
    run:
        import pandas
        # Sum the quantification files provided by beers per packet and then demultiplex by sample
        input_quantfiles = pathlib.Path(params.beers_dir).glob("run_run1/library_prep_pipeline/logs/Pipeline/*/library_prep_pipeline_molecule_pkt*.quant_file")
        output_quantfiles = pathlib.Path(params.beers_dir).glob("run_run1/library_prep_pipeline/data/*/library_prep_pipeline_result_molecule_pkt*.quant_file")

        input_quants_list = []
        for input_file in input_quantfiles:
            # Load all the packet quants
            input_quant = pandas.read_csv(input_file, sep="\t", header=None, index_col=0)
            input_quant.index.name = "BeersTranscriptID"
            input_quant.columns = ['count']
            input_quants_list.append(input_quant)
        input_quants = pandas.concat(input_quants_list, axis=0)
        for sample, filepath in zip(sample_ids, output.input_quants):
            # De-multiplex the samples
            # Transcripts are identified with sample# at the start
            prefix = f"{sample}_"
            sample_input_quants = input_quants[input_quants.index.str.startswith(prefix)].reset_index()
            sample_input_quants.groupby('BeersTranscriptID').sum().to_csv(filepath, sep="\t")

        output_quants_list = []
        for output_file in output_quantfiles:
            # Load all the packet quants
            output_quant = pandas.read_csv(output_file, sep="\t", header=None, index_col=0)
            output_quant.index.name = "BeersTranscriptID"
            output_quant.columns = ['count']
            output_quants_list.append(output_quant)
        output_quants = pandas.concat(output_quants_list, axis=0)
        for sample, filepath in zip(sample_ids, output.output_quants):
            # De-multiplex the samples
            # Transcripts are per-sample and identified with sample# at the start
            prefix = f"{sample}_"
            sample_output_quants = output_quants[output_quants.index.str.startswith(prefix)].reset_index()
            sample_output_quants.groupby('BeersTranscriptID').sum().to_csv(filepath, sep="\t")

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

rule strip_polya_tails:
    # Since we believe that excess simulated PolyA tail could be throwing off Salmon,
    # this step can be used to remove it. However, it did not appear to improve Salmon
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/finished_flag",
    output:
        directory("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers_no_polya/"),
    params:
        data_folder = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/run_run1/controller/data",
    resources:
        mem_mb = 6_000
    run:
        import re
        TAIL_LENGTH = 15# Check last 15 bases for a tail
        polya = re.compile(r"(T+)$")
        outdir = pathlib.Path(output[0])
        outdir.mkdir(exist_ok=True)
        outidr = outdir / "beers"
        outdir.mkdir(exist_ok=True)
        outdir = outdir / "run_run1"
        outdir.mkdir(exist_ok=True)
        outdir = outdir / "controller"
        outdir.mkdir(exist_ok=True)
        outdir = outdir / "data"
        outdir.mkdir(exist_ok=True)
        for sample in sample_ids:
            for lane in lanes_used:
                R1_path = pathlib.Path(params.data_folder) / f"S{sample}_L{lane}_R1.fastq"
                R2_path = pathlib.Path(params.data_folder) / f"S{sample}_L{lane}_R2.fastq"
                R1_out_path = outdir / f"S{sample}_L{lane}_R1.fastq"
                R2_out_path = outdir / f"S{sample}_L{lane}_R2.fastq"
                print(f"Running {sample} {lane}")
                shell(f"python scripts/strip_polya_tails.py {R1_path} {R2_path} {R1_out_path} {R2_out_path}")



rule run_salmon:
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/finished_flag",
        #"index/baby_mouse"
    output:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}/quant.sf",
    params:
        beers_dir = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/",
        data_folder = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/run_run1/controller/data",
        gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf",
        out_folder = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}/",
        index = "/project/itmatlab/index/SALMON-1.4.0_indexes/Mus_musculus.GRCm38.75/Mus_musculus.GRCm38.75/",
    threads: 6
    resources:
        mem_mb=25000,
    run:
        args = "-l A --softclip --softclipOverhangs -p 6"
        if wildcards.GC_correct == 'yes':
            args += ' --gcBias'
        if wildcards.Pos_correct == 'yes':
            args += ' --posBias'
        #TODO: handle --seqBias option
        shell(f"salmon quant -i {params.index} -g {params.gtf_file} {args} -1 {params.data_folder}/S{wildcards.sample}_*_R1.fastq -2 {params.data_folder}/S{wildcards.sample}_*_R2.fastq -o {params.out_folder}")


rule gather_quants:
    input:
        salmon_quants = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}/quant.sf",
                        GC_bias = GC_bias.keys(),
                        pos_3prime_bias = pos_3prime_bias.keys(),
                        GC_correct = GC_corrects,
                        Pos_correct = Pos_corrects,
                        sample = sample_ids),
        beers_input_quants = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/input_quant.txt",
                        GC_bias = GC_bias.keys(),
                        pos_3prime_bias = pos_3prime_bias.keys(),
                        sample=sample_ids),
        beers_output_quants = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/output_quant.txt",
                        GC_bias = GC_bias.keys(),
                        pos_3prime_bias = pos_3prime_bias.keys(),
                        sample=sample_ids),
        gene_annotations = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf",
    output:
        salmon = "results/salmon_quants.txt",
        beers_input = "results/beers_input_quants.txt",
        beers_output = "results/beers_output_quants.txt",
    resources:
        mem_mb = 12_000
    run:
        # Load transcript ID to gene ID mappings
        import re
        transcript_to_gene = {}
        with open(input.gene_annotations) as annotations:
            gene_id_re = re.compile(r'gene_id "([a-zA-Z0-9]*)";')
            transcript_id_re = re.compile(r'transcript_id "([a-zA-Z0-9]*)";')
            for line in annotations:
                gene_id = gene_id_re.search(line)
                transcript_id = transcript_id_re.search(line)
                if gene_id and transcript_id:
                    transcript_to_gene[transcript_id.groups()[0]] = gene_id.groups()[0]

        import pandas
        salmon_quants = []
        for quant_file, (GC_bias_, pos_3prime_bias_, GC_correct, Pos_correct, sample) in zip(input.salmon_quants, itertools.product(GC_bias.keys(), pos_3prime_bias.keys(), GC_corrects, Pos_corrects, sample_ids)):
            quants = pandas.read_csv(quant_file, sep="\t", index_col=0)
            #Name	Length	EffectiveLength	TPM	NumReads
            quants.index.name = "TranscriptID"
            quants['GC_bias'] = GC_bias_
            quants['sample'] = sample
            quants['pos_3prime_bias'] = pos_3prime_bias_
            quants['GC_correct'] = GC_correct == 'yes'
            quants['Pos_correct'] = Pos_correct == 'yes'
            salmon_quants.append(quants)
        salmon_quants = pandas.concat(salmon_quants, axis=0)
        salmon_quants['GeneID'] = salmon_quants.index.map(transcript_to_gene)
        salmon_quants.to_csv(output.salmon, sep="\t")

        beers_input_quants = []
        for quant_file, (GC_bias_, pos_3prime_bias_, sample) in zip(input.beers_input_quants, itertools.product(GC_bias.keys(), pos_3prime_bias.keys(), sample_ids)):
            quants = pandas.read_csv(quant_file, sep="\t", index_col=0)
            quants['sample'] = sample
            quants['GC_bias'] = GC_bias_
            quants['pos_3prime_bias'] = pos_3prime_bias_
            beers_input_quants.append(quants)
        beers_input_quants = pandas.concat(beers_input_quants, axis=0)
        beers_input_quants['TranscriptID'] = beers_input_quants.index.map(util.strip_beers_transcript_id)
        beers_input_quants['GeneID'] = beers_input_quants.TranscriptID.map(transcript_to_gene)
        beers_input_quants.to_csv(output.beers_input, sep="\t")

        beers_output_quants = []
        for quant_file, (GC_bias_, pos_3prime_bias_, sample) in zip(input.beers_output_quants, itertools.product(GC_bias.keys(), pos_3prime_bias.keys(), sample_ids)):
            quants = pandas.read_csv(quant_file, sep="\t", index_col=0)
            quants['sample'] = sample
            quants['GC_bias'] = GC_bias_
            quants['pos_3prime_bias'] = pos_3prime_bias_
            beers_output_quants.append(quants)
        beers_output_quants = pandas.concat(beers_output_quants, axis=0)
        beers_output_quants['TranscriptID'] = beers_output_quants.index.map(util.strip_beers_transcript_id)
        beers_output_quants['GeneID'] = beers_output_quants.TranscriptID.map(transcript_to_gene)
        beers_output_quants.to_csv(output.beers_output, sep="\t")

rule compare_accuracy:
    input:
        salmon = "results/salmon_quants.txt",
        beers_in = "results/beers_input_quants.txt",
        beers_out = "results/beers_output_quants.txt",
    output:
        dir = directory("results/accuracy/")
    resources:
        mem_mb = 24_000
    script:
        "scripts/compare_accuracy.py"

rule compute_gc_content:
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/finished_flag"
    output:
        gc_content = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/gc_content.txt",
    params:
        beers_output_data_folder = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/run_run1/controller/data",
        GC_biases = list(GC_bias.keys()),
        pos_3prime_biases = list(pos_3prime_bias.keys()),
        samples = sample_ids,
    resources:
        mem_mb = 18_000,
    script:
        "scripts/compute_gc_content.py"

rule gather_gc_content:
    input:
        gc_content = expand("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/gc_content.txt", GC_bias=GC_bias.keys(), pos_3prime_bias=pos_3prime_bias.keys())
    output:
        gc_content = "results/gc_content.txt",
        gc_summary = "results/gc_summary.txt"
    params:
        GC_biases = GC_bias.keys(),
        pos_3prime_bias = pos_3prime_bias.keys()
    run:
        import pandas
        import numpy
        res = []
        for gc_content in input.gc_content:
            d = pandas.read_csv(gc_content, sep="\t")
            res.append(d)
        data = pandas.concat(res)
        data.to_csv(output.gc_content, index=False)

        bins = numpy.linspace(0,1, 101)
        bins[-1] = 1.001 # make it right-inclusive on the last bin
        def get_GC_content(data):
            distribution,_ = numpy.histogram(data.GC_content, bins)
            return pandas.Series(distribution / numpy.sum(distribution), index=bins[:-1])
        GC_summary = data.groupby(["GC_bias", "pos_3prime_bias", "sample"]).apply(get_GC_content)
        GC_summary.to_csv(output.gc_summary, sep="\t")

rule plot_gc_content:
    input:
        gc_summary = "results/gc_summary.txt",
    output:
        directory("results/gc_content/")
    script:
        "scripts/plot_gc_content.py"

rule run_STAR:
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/finished_flag",
    output:
        out_dir = directory("data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/STAR/"),
    params:
        beers_dir = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/",
        data_folder = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/run_run1/controller/data",
        gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf",
        index = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/",
    threads: 6
    resources:
        mem_mb=40000,
    run:
        args = "--runThreadN 6"
        #R1s = ','.join(str(x) for x in pathlib.Path(params.data_folder).glob(f"S{wildcards.sample}_*_R1.fastq"))
        #R2s = ','.join(str(x) for x in pathlib.Path(params.data_folder).glob(f"S{wildcards.sample}_*_R2.fastq"))
        R1s = f"{params.data_folder}/S{wildcards.sample}_L1_R1.fastq" # Only use the first lane - star doesn't seem to like combining these files
        R2s = f"{params.data_folder}/S{wildcards.sample}_L1_R2.fastq"
        cmd = f"STAR --genomeDir {params.index} {args} --readFilesIn {R1s} {R2s} --outFileNamePrefix {output.out_dir}/STAR"
        print(cmd)
        shell(cmd)


rule generate_BAM:
    input:
        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/finished_flag",
    output:
        bam = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.bam",
    params:
        data_folder = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/beers/run_run1/controller/data",
    run:
        import uuid
        tmpdir = pathlib.Path(resources.tmpdir) / str(uuid.uuid4())
        tmpdir.mkdir(exist_ok=True)
        files = pathlib.Path(params.data_folder).glob(f"S{wildcards.sample}_L*.sam")
        print(f"Temp dir: {tmpdir}")
        for sam in files:
            print(f"Sorting on {sam}")
            shell(f"samtools sort -o {tmpdir}/{sam.name}.sorted.bam -O bam {sam}")
        shell(f"samtools merge {output.bam} {tmpdir}/*.sorted.bam")

rule generate_BAI:
    input:
        bam = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.bam",
    output:
        bai = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.bam.bai",
    shell:
        "samtools index -b {input} {output}"

rule plot_coverage:
    input:
        bams = expand('data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.bam',
                        GC_bias = GC_bias.keys(),
                        pos_3prime_bias = pos_3prime_bias.keys(),
                        sample = sample_ids),
    output:
        directory("results/coverage/")
    script:
        "scripts/plot_coverage.py"

# These rules were included for coverage plot generation but now we go straight from BAM files
#rule generate_bigwig:
#    input:
#        bam = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.bam",
#        bai = "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.bam.bai",
#    output:
#        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.forward.bigwig",
#        "data/GC_bias={GC_bias}.pos_3prime_bias={pos_3prime_bias}/sample{sample}/BEERS_output.reverse.bigwig",
#    shell:
#        "bamCoverage --filterRNAstrand forward --binSize 1 -b {input.bam} -o {output[0]} && bamCoverage --filterRNAstrand reverse -b {input.bam} -o {output[1]}"
#
#rule rename_chromes_in_gtf:
#    # Rename chromosomes in a gtf file as 'chr1',...'chrM' instead of '1', ..., 'MT'
#    # Used for coolbox plotting, but also just consistent with the generated camparee data
#    input:
#        gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf",
#    output:
#        renamed_gtf_file = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.renamed_chromes.gtf",
#    run:
#        with open(input.gtf_file) as input_gtf:
#            with open(output.renamed_gtf_file, "w") as output_gtf:
#                for line in input_gtf:
#                    if line.startswith("#"):
#                        output_gtf.write(line)
#                    elif line.startswith("MT"):
#                        output_gtf.write("chrM" + line[2:])
#                    else:
#                        output_gtf.write("chr" + line)
