def strip_beers_transcript_id(beers_id):
    # Beers gives transcript IDs with extra information
    # the format is:
    # {sample}_{transcript_id}_{allele}
    return beers_id.split("_")[1]
