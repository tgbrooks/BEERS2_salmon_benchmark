import pathlib
import io
import subprocess as subp
import coolbox
import coolbox.api as cb
import numpy
import pandas
import matplotlib.patches as patches
from dna_features_viewer import GraphicFeature, GraphicRecord

outdir = pathlib.Path(snakemake.output[0])
outdir.mkdir(exist_ok=True)

#genome_range = "chr6:125161854-125166467" # Gapdh
#genome_range = cb.GenomeRange("chr6:125160000-125170000") # Gapdh
#genome_range = cb.GenomeRange("6:125160000-125170000") # Gapdh
#genome_range = cb.GenomeRange("chr5:90460897-90476602") # Alb
#genome_range = cb.GenomeRange("chr9:110281287-110281409") # Mir6236
#genome_range = cb.GenomeRange("9:110281287-110281409") # Mir6236
#genome_range = cb.GenomeRange("9:110280000-110282000") # Mir6236
#genome_range = cb.GenomeRange("12:100199435-100209814")
genome_range = cb.GenomeRange("chr1:768000-816000") # Cflar
gtf_path = "/home/thobr/BEERS2/CAMPAREE/resources/baby_genome.mm10/baby_genome.mm10.annotation.gtf"
#gtf_path = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.renamed_chromes.gtf"
#gtf_path = "/project/itmatlab/index/STAR-2.7.6a_indexes/GRCm38.ensemblv102/Mus_musculus.GRCm38.102.gtf"
bam_paths = {
    (gc_bias, pos_bias): f"data/GC_bias={gc_bias}.pos_3prime_bias={pos_bias}/sample1/BEERS_output.bam"
        for gc_bias in ['none', 'med', 'high']
        for pos_bias in ['none', 'med', 'high']
}

class CustomGTF(cb.GTF):
    CDS_height = 1
    noncoding_height = 0.3
    intron_height = 0.03
    transcript_padding = 0.2
    min_gap = 50 # in base-pairs
    exon_color = "black"
    def plot(self, ax, gr: cb.GenomeRange, **kwargs):
        self.ax = ax
        df = self.fetch_plot_data(gr)
        if self.has_prop("row_filter"):
            filters = self.properties["row_filter"]
            for filter_ in filters.split(";"):
                if filter_ == '':
                    continue
                df = df.query(filter_)
        region_length = gr.end - gr.start
        len_ratio_th = self.properties["length_ratio_thresh"]
        df = df[(df["end"] - df["start"]) > region_length * len_ratio_th]

        # Group by transcripts
        transcript_ids = df.attribute.str.extract(r'transcript_id "(.*?)";').iloc[:,0].dropna().unique()
        #features = []
        max_y = 0
        furthest_right_per_line = {}
        for transcript in transcript_ids:
            transcript_df = df[df.attribute.str.contains(transcript)]
            transcript_left = min(transcript_df['start'])
            transcript_right = max(transcript_df['end'])
            color = numpy.random.choice(self.colors)
            # Use the first line we fit in
            # or start a new line if we don't fit
            i = min((i for i, right_end in furthest_right_per_line.items() if right_end + self.min_gap < transcript_left), default=len(furthest_right_per_line) + 1)
            y = (self.CDS_height + self.transcript_padding)*i
            furthest_right_per_line[i] = transcript_right
            for _, row in transcript_df.query("feature == 'CDS'").iterrows():
                # Draw the CDSs
                ax.add_patch(patches.Rectangle(
                    (row['start'], y),row['end'] - row['start'], self.CDS_height,
                    facecolor = self.exon_color,
                ))
            for _, row in transcript_df.query("feature.isin(['five_prime_utr', 'three_prime_utr'])").iterrows():
                # Draw the CDSs
                ax.add_patch(patches.Rectangle(
                    (row['start'], y),row['end'] - row['start'], self.noncoding_height,
                    facecolor = self.exon_color,
                ))
            max_y = (self.CDS_height + self.transcript_padding) * len(furthest_right_per_line)
            # Draw the introns
            ax.add_patch(patches.Rectangle(
                (transcript_left, y),
                transcript_right - transcript_left,
                self.intron_height,
                facecolor = self.exon_color,
            ))
        ax.set_ylim(0, max_y)
        ax.set_xlim(gr.start, gr.end)
        self.plot_label()

def depth_by_samtools(bam_path, region):
    cmd = ["samtools", "depth", bam_path, "-r", region, '-s', '-d', '0', '-a']
    p = subp.Popen(cmd, stdout=subp.PIPE)
    results = p.stdout.read().decode('utf-8')
    if results == '':
        return pandas.DataFrame([[genome_range.chrom, genome_range.start, 0]], columns=["chr", "pos", "score"])
    covs = pandas.read_csv(io.StringIO(results), sep="\\s+", header=None)
    covs.columns = ['chr', 'pos', 'score']
    return covs

class BAMDepth(cb.HistBase):
    """
    Alignment read depth track. Per-base read levels.
    Parameters
    ----------
    file: str
        File path of bam file.
    """

    DEFAULT_PROPERTIES = {
        "height": 3,
        "style": cb.HistBase.STYLE_FILL,
        "color": "#6688ff",
    }

    def __init__(self, file, **kwargs):
        properties = BAMDepth.DEFAULT_PROPERTIES.copy()
        properties.update({
            'file': file,
            **kwargs,
        })
        super().__init__(**properties)
        self.indexed_bam = coolbox.utilities.bam.process_bam(file)

    def fetch_data(self, gr: cb.GenomeRange, **kwargs) -> pandas.DataFrame:
        return depth_by_samtools(
            self.indexed_bam,
            str(genome_range),
        )

# Load GTF data and process so that transcripts are on separate lines
gtf_gene = cb.GTF(gtf_path)

gtf_transcripts = CustomGTF(gtf_path, row_filter='')

bams = {index: BAMDepth(bam_path) for index, bam_path in bam_paths.items()}
#bws = {(gc_bias, pos_bias): (cb.BigWig(f"data/GC_bias={gc_bias}.pos_3prime_bias={pos_bias}/sample1/BEERS_output.forward.bigwig"),
#                                cb.BigWig(f"data/GC_bias={gc_bias}.pos_3prime_bias={pos_bias}/sample1/BEERS_output.reverse.bigwig"))
#                    for gc_bias in ['none', 'med', 'high'] for pos_bias in ['none', 'med', 'high']}

frame =  cb.XAxis() + gtf_gene + cb.TrackHeight(2) + gtf_transcripts + cb.TrackHeight(5) + cb.XAxis()
#for ((gc_bias, pos_bias), (fwd_bw, rev_bw)) in bws.items():
for ((gc_bias, pos_bias), bam) in bams.items():
    frame += bam
    #frame += fwd_bw
    frame += cb.Title(f"GC_bias: {gc_bias}\nPos bias: {pos_bias}")
    #frame += cb.Color('#000000')
    frame += cb.MinValue(0)
    frame += cb.MaxValue(bam.fetch_data(genome_range).score.max())
    #frame += rev_bw
    ##frame += cb.Color('#FF00000')
    #frame += cb.MinValue(0)
    #frame += cb.MaxValue(rev_bw.fetch_data(genome_range).score.max())

plot = frame.plot(genome_range)
plot.savefig(outdir / "coverage.png", dpi=300)
