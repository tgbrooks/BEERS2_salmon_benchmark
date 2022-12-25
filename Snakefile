import pathlib
import itertools
import util

configfile:
    "config.yaml"

cfg = config
run_configs = config['run_configs']
run_ids = run_configs.keys()
lanes_used = config['lanes_used']
sample_ids = config['sample_ids']

RCLONE_REMOTE = "aws_igv_data" # Must first run `rclone config` and set up a remote with this name for uploading trackhubs to
BUCKET_NAME = "itmat.igv.data" # Bucket to upload to with rclone for trackhubs
BUCKET_DIR = f"BEERS2_SALMON_BENCHMARK/DATA"


BABY_GENOME = False
if BABY_GENOME:
    CAMPAREE_CONFIG = 'config/baby_genome.camparee_config.yaml'
    REFERENCE_GENOME = '/home/thobr/BEERS2/CAMPAREE/resources/baby_genome.mm10/baby_genome.mm10.oneline_seqs.fa'
    CAMPAREE_OUTPUT = '/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/baby_genome/run_1/CAMPAREE/data/'
    CAMPAREE_DIR = '/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/baby_genome/run_1/'
    SALMON_INDEX = 'index/baby_mouse'
    # NOTE: we use the full annotations since baby genome is a subset. This is wasteful but works
    GENE_ANNOTATIONS = '/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/Mus_musculus.GRCm38.93.gtf.gz'
    STAR_INDEX = '/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/'
    GENOME = "baby_mouse"
else:
    CAMPAREE_CONFIG = 'config/camparee_config.yaml'
    REFERENCE_GENOME = '/project/itmatlab/for_tom/BEERS2_benchmark/resources/mm10/Mus_musculus.GRCm38.dna.oneline_seqs.fa'
    CAMPAREE_OUTPUT = '/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/run_1/CAMPAREE/data/'
    CAMPAREE_DIR = '/project/itmatlab/for_tom/BEERS2_benchmark/CAMPAREE_out/run_1/'
    SALMON_INDEX = '/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/index/'
    GENE_ANNOTATIONS = '/project/itmatlab/index/SALMON-1.9.0_indexes/GRCm38.ensemblv93/Mus_musculus.GRCm38.93.gtf.gz'
    STAR_INDEX = '/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/' #TODO: switch this to v93 too? Do we still need this command?
    GENOME = "mm10"

LATENCY_WAIT = 60 # Seconds


GC_corrects = ['yes', 'no']
Pos_corrects = ['yes', 'no']
Seq_corrects = ['yes', 'no']

wildcard_constraints:
    run = "[A-Za-z_0-9]+",
    sample = "[0-9]+",
    GC_correct = "[A-Za-z]+",
    Pos_correct = "[A-Za-z]+",
    Seq_correct = "[A-Za-z]+",

rule all:
    input:
        "results/seq_bias/fwd_seq_frequencies.png",
        "results/accuracy/",
        "results/salmon_quants.txt",
        "results/gc_content/",
        #"results/coverage/",
        "results/igv/url.txt",
        beers_input_quants = expand("data/{run}/sample{sample}/input_quant.txt",
                    run = run_ids,
                    sample = sample_ids),
        beers_output_quants = expand("data/{run}/sample{sample}/output_quant.txt",
                    run = run_ids,
                    sample = sample_ids),
        beers_output_bam = expand("data/{run}/sample{sample}/BEERS_output.bam",
                    run = run_ids,
                    sample = sample_ids),

rule run_CAMPAREE:
    input:
        CAMPAREE_CONFIG
    output:
        directory(CAMPAREE_DIR)
    params:
        CAMPAREE_DIR
    resources:
        mem_mb = 10_000,
    shell:
        'rm -r {params[0]} && run_camparee.py -r 1 -d -c {input[0]}'

rule run_beers:
    input:
        config_template = "config/config.template.yaml",
        camparee_output = CAMPAREE_OUTPUT,
        reference_genome = REFERENCE_GENOME,
    output:
        flag_file = "data/{run}/beers/finished_flag"
    params:
        beers_dir = directory("data/{run}/beers/"),
        config = "data/{run}/beers.config.yaml",
    resources:
        mem_mb = 6_000
    run:
        #Generate config template
        import string
        config_template = string.Template(open(input.config_template, "r").read())
        config_fill_values = dict()
        config_fill_values.update(run_configs[wildcards.run])
        config_fill_values['camparee_output'] = input.camparee_output
        config_fill_values['reference_genome'] = input.reference_genome
        config_fill_values.update(cfg['pos_3prime_bias'][config_fill_values['pos_3prime_bias']])
        config_fill_values.update(cfg['GC_bias'][config_fill_values['GC_bias']])
        config_fill_values.update(cfg['primer_bias'][config_fill_values['primer_bias']])
        config = config_template.substitute(config_fill_values)
        with open(params.config, "w") as f:
            f.write(config)
        # Run beers
        beers_cmd =  f"run_beers --configfile {params.config} --directory {params.beers_dir} --profile lsf -j 30 -c 30 --latency-wait {LATENCY_WAIT}"
        print(beers_cmd)
        shell(beers_cmd)
        shell("touch {output.flag_file}")

rule get_beers_quant_files:
    input:
        "data/{run}/beers/finished_flag",
    params:
        beers_dir = "data/{run}/beers/library_prep_pipeline/sample{sample}/",
    output:
        input_quants = "data/{run}/sample{sample}/input_quant.txt",
        output_quants = "data/{run}/sample{sample}/output_quant.txt",
    resources:
        mem_mb = 6_000,
    run:
        sample = wildcards.sample
        import pandas
        # Sum the quantification files provided by beers per packet and then demultiplex by sample
        input_quantfiles = list(pathlib.Path(params.beers_dir).glob("from_*/input_quants/library_prep_pipeline_molecule_pkt*.quant_file"))
        output_quantfiles = list(pathlib.Path(params.beers_dir).glob("from_*/library_prep_pipeline_result_molecule_pkt*.quant_file"))
        print(f"Processing {len(input_quantfiles)} input quant files and {len(output_quantfiles)} output quant files")

        input_quants_list = []
        for input_file in input_quantfiles:
            # Load all the packet quants
            input_quant = pandas.read_csv(input_file, sep="\t", index_col=0)
            input_quant.index.name = "BeersTranscriptID"
            input_quant.columns = ['count']
            input_quants_list.append(input_quant)
        input_quants = pandas.concat(input_quants_list, axis=0).reset_index()
        input_quants.groupby('BeersTranscriptID')['count'].sum().to_csv(output.input_quants, sep="\t")
        #input_quants['mature'] = ~input_quants.index.str.endswith("_pre_mRNA")
        #pandas.DataFrame({
        #    "count": input_quants.groupby('BeersTranscriptID')['count'].sum(),
        #    "mature_mRNA_count": input_quants.query("mature").groupby('BeersTranscriptID')['count'].sum(),
        #}).fillna(0).to_csv(output.input_quants, sep="\t")

        output_quants_list = []
        for output_file in output_quantfiles:
            # Load all the packet quants
            output_quant = pandas.read_csv(output_file, sep="\t", index_col=0)
            output_quant.index.name = "BeersTranscriptID"
            output_quant.columns = ['count']
            output_quants_list.append(output_quant)
        output_quants = pandas.concat(output_quants_list, axis=0).reset_index()
        output_quants.groupby('BeersTranscriptID')['count'].sum().to_csv(output.output_quants, sep="\t")
        #output_quants['mature'] = ~output_quants.index.str.endswith("_pre_mRNA")
        #pandas.DataFrame({
        #    "count": output_quants.groupby('BeersTranscriptID')['count'].sum(),
        #    "mature_mRNA_count": output_quants.query("mature").groupby('BeersTranscriptID')['count'].sum(),
        #}).fillna(0).to_csv(output.output_quants, sep="\t")

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
        "data/{run}/beers/finished_flag",
    output:
        directory("data/{run}/beers_no_polya/"),
    params:
        data_folder = "data/{run}/beers/results",
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
        outdir = outdir / "results"
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
        "data/{run}/beers/finished_flag",
        "index/baby_mouse"
    output:
        "data/{run}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}.Seq_correct={Seq_correct}/quant.sf",
    params:
        beers_dir = "data/{run}/beers/",
        data_folder = "data/{run}/beers/results",
        gtf_file = GENE_ANNOTATIONS,
        out_folder = "data/{run}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}.Seq_correct={Seq_correct}/",
        index = SALMON_INDEX,
    threads: 6
    resources:
        mem_mb=25000,
    run:
        args = "-l A --softclip --softclipOverhangs -p 6"
        if wildcards.GC_correct == 'yes':
            args += ' --gcBias'
        if wildcards.Pos_correct == 'yes':
            args += ' --posBias'
        if wildcards.Seq_correct == 'yes':
            args += ' --seqBias'
        #TODO: handle --seqBias option
        shell(f"salmon quant -i {params.index} -g {params.gtf_file} {args} -1 {params.data_folder}/S{wildcards.sample}_*_R1.fastq -2 {params.data_folder}/S{wildcards.sample}_*_R2.fastq -o {params.out_folder} --numBootstraps 50")


rule gather_quants:
    input:
        salmon_quants = expand("data/{run}/sample{sample}/salmon/GC_correct={GC_correct}.Pos_correct={Pos_correct}.Seq_correct={Seq_correct}/quant.sf",
                        run = run_ids,
                        GC_correct = GC_corrects,
                        Pos_correct = Pos_corrects,
                        Seq_correct = Seq_corrects,
                        sample = sample_ids),
        beers_input_quants = expand("data/{run}/sample{sample}/input_quant.txt",
                        run = run_ids,
                        sample=sample_ids),
        beers_output_quants = expand("data/{run}/sample{sample}/output_quant.txt",
                        run = run_ids,
                        sample=sample_ids),
        gene_annotations = GENE_ANNOTATIONS,
    output:
        salmon = "results/salmon_quants.txt",
        beers_input = "results/beers_input_quants.txt",
        beers_output = "results/beers_output_quants.txt",
    resources:
        mem_mb = 12_000
    run:
        # Load transcript ID to gene ID mappings
        import gzip
        import re
        transcript_to_gene = {}
        open_func = gzip.open if input.gene_annotations.endswith(".gz") else open
        with open_func(input.gene_annotations, mode="rt") as annotations:
            gene_id_re = re.compile(r'gene_id "([a-zA-Z0-9]*)";')
            transcript_id_re = re.compile(r'transcript_id "([a-zA-Z0-9]*)";')
            for line in annotations:
                gene_id = gene_id_re.search(line)
                transcript_id = transcript_id_re.search(line)
                if gene_id and transcript_id:
                    transcript_to_gene[transcript_id.groups()[0]] = gene_id.groups()[0]

        import pandas
        salmon_quants = []
        for quant_file, (run, GC_correct, Pos_correct, Seq_correct, sample) in zip(input.salmon_quants, itertools.product(run_ids, GC_corrects, Pos_corrects, Seq_corrects, sample_ids), strict=True):
            quants = pandas.read_csv(quant_file, sep="\t", index_col=0)
            #Name	Length	EffectiveLength	TPM	NumReads
            quants.index = quants.index.map(lambda x: x.split(".")[0]) # Some annotations have ".###" at the end to indicate the version number - we remove those
            quants.index.name = "TranscriptID"
            quants['run'] = run
            quants['GC_bias'] = run_configs[run]['GC_bias']
            quants['pos_3prime_bias'] = run_configs[run]['pos_3prime_bias']
            quants['primer_bias'] = run_configs[run]['primer_bias']
            quants['sample'] = sample
            quants['GC_correct'] = GC_correct == 'yes'
            quants['Pos_correct'] = Pos_correct == 'yes'
            quants['Seq_correct'] = Seq_correct == 'yes'
            salmon_quants.append(quants)
        salmon_quants = pandas.concat(salmon_quants, axis=0)
        salmon_quants['GeneID'] = salmon_quants.index.map(transcript_to_gene)
        salmon_quants.to_csv(output.salmon, sep="\t")

        beers_input_quants = []
        for quant_file, (run, sample) in zip(input.beers_input_quants, itertools.product(run_ids, sample_ids)):
            quants = pandas.read_csv(quant_file, sep="\t", index_col=0)
            quants['run'] = run
            quants['sample'] = sample
            quants['GC_bias'] = run_configs[run]['GC_bias']
            quants['pos_3prime_bias'] = run_configs[run]['pos_3prime_bias']
            quants['primer_bias'] = run_configs[run]['primer_bias']
            beers_input_quants.append(quants)
        beers_input_quants = pandas.concat(beers_input_quants, axis=0)
        beers_input_quants['TranscriptID'] = beers_input_quants.index.map(util.strip_beers_transcript_id)
        beers_input_quants['GeneID'] = beers_input_quants.TranscriptID.map(transcript_to_gene)
        beers_input_quants.to_csv(output.beers_input, sep="\t")

        beers_output_quants = []
        for quant_file, (run, sample) in zip(input.beers_output_quants, itertools.product(run_ids, sample_ids)):
            quants = pandas.read_csv(quant_file, sep="\t", index_col=0)
            quants['run'] = run
            quants['sample'] = sample
            quants['GC_bias'] = run_configs[run]['GC_bias']
            quants['pos_3prime_bias'] = run_configs[run]['pos_3prime_bias']
            quants['primer_bias'] = run_configs[run]['primer_bias']
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
        "data/{run}/beers/finished_flag"
    output:
        gc_content = "data/{run}/gc_content.txt",
    params:
        beers_output_data_folder = "data/{run}/beers/results",
        samples = sample_ids,
    resources:
        mem_mb = 18_000,
    script:
        "scripts/compute_gc_content.py"

rule gather_gc_content:
    input:
        gc_content = expand("data/{run}/gc_content.txt", run = run_ids),
    output:
        gc_content = "results/gc_content.txt",
    run:
        import pandas
        import numpy
        res = []
        for gc_content in input.gc_content:
            d = pandas.read_csv(gc_content, sep="\t")
            res.append(d)
        data = pandas.concat(res)
        data.to_csv(output.gc_content, index=False, sep="\t")

rule plot_gc_content:
    input:
        gc_summary = "results/gc_content.txt",
    output:
        directory("results/gc_content/")
    script:
        "scripts/plot_gc_content.py"

rule run_STAR:
    input:
        "data/{run}/beers/finished_flag",
    output:
        out_dir = directory("data/{run}/sample{sample}/STAR/"),
    params:
        beers_dir = "data/{run}/beers/",
        data_folder = "data/{run}/beers/results",
        gtf_file = GENE_ANNOTATIONS,
        index = STAR_INDEX,
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
        "data/{run}/beers/finished_flag",
    output:
        bam = "data/{run}/sample{sample}/BEERS_output.bam",
    params:
        data_folder = "data/{run}/beers/results",
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
        bam = "data/{run}/sample{sample}/BEERS_output.bam",
    output:
        bai = "data/{run}/sample{sample}/BEERS_output.bam.bai",
    shell:
        "samtools index -b {input} {output}"

#rule plot_coverage:
#    input:
#        bams = expand('data/{run}/sample{sample}/BEERS_output.bam',
#                        run = run_ids,
#                        sample = sample_ids),
#    output:
#        directory("results/coverage/")
#    script:
#        "scripts/plot_coverage.py"

rule compute_seq_bias:
    input:
        flag_file = expand("data/{run}/beers/finished_flag",
                        run = run_ids,
                    ),
    output:
        seq_frequencies = "results/seq_bias/seq_frequencies.json",
    script:
        "scripts/compute_seq_bias.py"

rule plot_seq_bias:
    input:
        seq_frequencies = "results/seq_bias/seq_frequencies.json",
    output:
        fwd_frequencies = "results/seq_bias/fwd_seq_frequencies.png",
        rev_frequencies = "results/seq_bias/rev_seq_frequencies.png",
    script:
        "scripts/plot_seq_bias.py"


rule make_igv_view:
    input:
        bams = expand("data/{run}/sample{sample}/BEERS_output.bam", run=run_ids, sample=sample_ids),
        bais = expand("data/{run}/sample{sample}/BEERS_output.bam.bai", run=run_ids, sample=sample_ids),
        camparee_config = CAMPAREE_CONFIG
    output:
        "results/igv/url.txt",
    params:
        out_dir = "results/igv/",
    run:
        import gzip
        out_dir = pathlib.Path(params.out_dir)
        out_dir.mkdir(exist_ok=True)

        import yaml
        camparee_config = yaml.safe_load(open(input.camparee_config))
        genome_dir = pathlib.Path(camparee_config['resources']['directory_path']) / camparee_config['resources']['species_model']
        genome_file = (genome_dir /camparee_config['resources']['reference_genome_filename']).resolve()

        genome_index_file = pathlib.Path("results/igv/genome.fai")
        if not genome_index_file.exists():
            shell(f"samtools faidx -o {genome_index_file} {genome_file}")


        # Indexed annotations
        genome_gtf = GENE_ANNOTATIONS
        bgz_gtf_file = pathlib.Path('results/igv/annotation.gtf.bgz')
        if not bgz_gtf_file.exists():
            shell(f'''(grep "^#" {genome_gtf}; grep -v "^#" {genome_gtf} | sort -t"`printf '\t'`" -k1,1 -k4,4n) | bgzip -c  > {bgz_gtf_file}''')
        index_gtf_file = pathlib.Path('results/igv/annotation.gtf.bgz.tbi')
        if not index_gtf_file.exists():
            shell(f"tabix -p gff {bgz_gtf_file}")

        print("Uploading to bucket via rclone")
        rclone_cmd = f"rclone copyto {genome_file} {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.fa"
        print(rclone_cmd)
        shell(rclone_cmd)
        rclone_cmd = f"rclone copyto {genome_index_file} {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.fai"
        print(rclone_cmd)
        shell(rclone_cmd)
        rclone_cmd = f"rclone copyto {bgz_gtf_file} {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.gtf.bgz"
        print(rclone_cmd)
        shell(rclone_cmd)
        rclone_cmd = f"rclone copyto {index_gtf_file} {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.gtf.bgz.tbi"
        print(rclone_cmd)
        shell(rclone_cmd)

        for run in run_ids:
            for sample_id in sample_ids:
                rclone_cmd = f"rclone copy data/{run}/sample{sample_id}/BEERS_output.bam {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/{run}/sample{sample_id}/"
                print(rclone_cmd)
                shell(rclone_cmd)
                rclone_cmd = f"rclone copy data/{run}/sample{sample_id}/BEERS_output.bam.bai {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/{run}/sample{sample_id}/"
                print(rclone_cmd)
                shell(rclone_cmd)

        # Generate the session object - it tells IGV.js how to load everything we generated
        genome = {
             "fastaURL": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.fa",
             "indexURL": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.fai",
             "tracks": [
                {
                    "name": "Genes",
                    "type": "annotation",
                    "format": "gtf",
                    "url": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.gtf.bgz",
                    "indexURL": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/GENOME/genome.gtf.bgz.tbi",
                    "nameField": "gene_id",
                    "altColor": "rgb(0,100,100)",
                },
            ],
        }

        session = {
            "genome": genome,
            "tracks": [
                {
                    "type": "alignment",
                    "name": f"{run_id}.sample{sample_id}",
                    "format": "bed",
                    "url": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/{run_id}/sample{sample_id}/BEERS_output.bam",
                    "indexURL": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/{run_id}/sample{sample_id}/BEERS_output.bam.bai",
                }
                for run_id in run_ids
                for sample_id in sample_ids
            ],
        }
        with open("results/igv/session.json", "w") as session_file:
            json.dump(session, session_file)
        rclone_cmd = f"rclone copyto results/igv/session.json {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/session.json"
        print(rclone_cmd)
        shell(rclone_cmd)

        # Genereate a simplified session showing only some data
        simple_session = {
            "genome": genome,
            "tracks": [
                {
                    "type": "alignment",
                    "name": f"{run_id}.sample{sample_id}",
                    "format": "bed",
                    "url": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/{run_id}/sample{sample_id}/BEERS_output.bam",
                    "indexURL": f"https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/{run_id}/sample{sample_id}/BEERS_output.bam.bai",
                }
                for run_id in ['unbiased', 'all_bias']
                for sample_id in [1]
            ],
        }
        with open("results/igv/simple_session.json", "w") as session_file:
            json.dump(simple_session, session_file)
        rclone_cmd = f"rclone copyto results/igv/simple_session.json {RCLONE_REMOTE}:{BUCKET_NAME}/{BUCKET_DIR}/simple_session.json"
        print(rclone_cmd)
        shell(rclone_cmd)


        # Write out the URL of the trackhub for easy use in the genome browser
        URLS = f'''
        https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/session.json\n
        https://s3.amazonaws.com/{BUCKET_NAME}/{BUCKET_DIR}/simple_session.json\n
        '''
        (out_dir / "url.txt").write_text(URLS)
        print("Uploaded sessions to:")
        print(URLS)
