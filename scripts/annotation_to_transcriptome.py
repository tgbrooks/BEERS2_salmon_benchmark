import argparse

import pandas

import camparee.camparee_utils
import beers_utils.general_utils

parser = argparse.ArgumentParser(description="Convert a genome and annotation file to a transcriptome fasta file")
parser.add_argument("--genome", help="genome fasta file")
parser.add_argument("--annotation", help="tab-separated one-line-per-transcript annotation file such as used by CAMPAREE")
parser.add_argument("--output", help="path to output the transcriptome fasta file")


args = parser.parse_args()

genome = camparee.camparee_utils.CampareeUtils.create_genome(args.genome)

annotation = pandas.read_csv(args.annotation, index_col="transcriptID", sep="\t")
# Contains the following columns:
#chrom	strand	txStart	txEnd	exonCount	exonStarts	exonEnds	transcriptID	geneID	geneSymbol	biotype

with open(args.output, 'w') as output:
    for transcript_id, transcript in annotation.iterrows():
        exon_starts = [int(x) for x in transcript['exonStarts'].split(",")]
        exon_ends = [int(x) for x in transcript['exonEnds'].split(",")]
        chrom = genome[transcript['#chrom']]
        strand = transcript.strand
        sequence = ''.join(chrom[start-1:end] for start, end in zip(exon_starts, exon_ends))

        if strand == '-':
            sequence = beers_utils.general_utils.GeneralUtils.create_complement_strand(sequence)
            
        output.write(f">{transcript_id} {transcript.geneID}\n")
        output.write(sequence)
        output.write("\n")
