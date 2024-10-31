import argparse
import pysam
import re
import itertools
from Bio import SeqIO, Seq
from Bio.Data.IUPACData import ambiguous_dna_values
from Bio.Align import PairwiseAligner

parser = argparse.ArgumentParser()

parser.add_argument('-i',
            '--infile',
            help='Input bamfile',
            type=str,
            required=True)

parser.add_argument('-o',
            '--outfile',
            help='Output bamfile',
            type=str,
            required=True)

parser.add_argument('-fq',
            '--outfastq',
            help='Output fastq file for trimmed reads',
            type=str,
            required=True)

parser.add_argument('-l', 
                '--minlength',
                help='Minimum length of fastq read after trimming allowed',
                default=None,
                type=int,
               required = False)

parser.add_argument('-p',
            '--primerfile',
            help='Input primer file in fasta format',
            type=str,
            required=True)

parser.add_argument('-r',
            '--referencefile',
            help='Input reference file in fasta format',
            type=str,
            required=True)

parser.add_argument('-fwd',
            '--fwdkey',
            help='Keyword that indicates forward primer',
            type=str,
            required=True)

parser.add_argument('-rev',
            '--revkey',
            help='Keyword that indicates reverse primer',
            type=str,
            required=True)

parser.add_argument('-x', 
                '--padding',
                help='Number of nucleotides allowed to be upstream of the primer while clipping',
                default=10,
                type=int,
               required = False)

parser.add_argument('-m', 
                '--mismatch',
                help='Number of mismatches allowed while matching and determining the position of the primers',
                default=2,
                type=int,
               required = False)

class Region(object):
    """Primer alignment region class"""
    def __init__(self, start, end, strand):
        self.start = int(start)
        self.end = int(end)
        self.strand = strand

def calculate_overlap(read, primer, padding):
    """Determines overlap between aligned read and primer object"""
    n_clip = 0
    if primer.strand == 'left':
        if read.get_overlap(primer.start, primer.end + 1) > 0:
            if (read.reference_start > (primer.start - padding)):
                n_clip = primer.end - read.reference_start
    elif primer.strand == 'right':
        if read.get_overlap(primer.start, primer.end) > 0:
            if (read.reference_end < (primer.end + padding)):
                n_clip = read.reference_end - primer.start + 1
    if n_clip < 0:
        raise ValueError("Something went wrong: Trying to clip a negative number of nucleotides in read "+read.query_name)
    return(n_clip)

def clip_read(read, n_clip, side):
    '''
    Softclips number of nucleotides of either left or right side of an aligned read.

    (For reference) Pysam cigar codes:
    M   BAM_CMATCH  0
    I   BAM_CINS    1
    D   BAM_CDEL    2
    N   BAM_CREF_SKIP   3 (not handled)
    S   BAM_CSOFT_CLIP  4
    H   BAM_CHARD_CLIP  5
    P   BAM_CPAD    6 (not handled)
    =   BAM_CEQUAL  7 (not handled)
    X   BAM_CDIFF   8 (not handled)
    B   BAM_CBACK   9 (not handled)
    '''
    
    #If all aligned bases are clipped completely clip cigarstring
    if read.qlen <= n_clip:
        read.cigartuples = [(4,read.query_length)]
        return

    current_cigar = read.cigartuples
    #Reverse cigar if we have to clip the right side
    if side == "right":
        current_cigar.reverse()

    #Expand cigarstring
    cigar_expanded = ''.join([j*str(i) for i,j in current_cigar])

    #Clip cigar until no more "clip" left:
    clip_leftover = n_clip
    n = 0
    cigar_replacement = ''
    while clip_leftover > 0:
        cig = int(cigar_expanded[n])
        if cig == 0:
            cigar_replacement += "4" #Replace matches with softclipped base
            clip_leftover -= 1
        elif cig == 1:
            cigar_replacement += "4" #Replace insertions with softclipped base
        elif cig == 2:
            n_clip += 1 #Increase n_clip to increase reference_start shift at the end in case of deletions
        elif cig == 4:
            cigar_replacement += "4" #Do not replace softclipped
        elif cig == 5:
            cigar_replacement += "5" #Do not replace hardclipped
        else:
            raise ValueError("Something went wrong: do not know how to clip " + str(cig) + " in cigarstring")
        n += 1

    cigar_expanded = cigar_replacement + cigar_expanded[n:]

    #Recreate tuples:
    clipped_cigar = []
    c = 0
    nprev=cigar_expanded[0]
    for n in cigar_expanded:
        if n==nprev:
            c+=1
        else:
            clipped_cigar.append((int(nprev),c))
            c=1
            nprev=n
    clipped_cigar.append((int(n),c))

    #Un-reverse cigar if clipped on right side
    if side == 'right':
        clipped_cigar.reverse()

    #Replace cigar tuples of the read
    read.cigartuples = clipped_cigar

    #Shift alignment start if clipped to the right
    if side == 'left':
        read.reference_start = read.reference_start + n_clip

def trim_read(args, read):
    '''
    Trims softclipped nucleotides of the reads from the left and right side and returns read in fastq format.
    '''

    read_name = read.qname
    read_sequence = read.query_sequence
    read_qualities = ''.join(map(lambda x: chr( x+33 ), read.query_qualities))

    #Trim left region
    if read.cigartuples[0][0] == 4:
        n_trim = read.cigartuples[0][1]
        read_sequence = read_sequence[n_trim:]
        read_qualities = read_qualities[n_trim:]
    
    #Trim right region
    if read.cigartuples[-1][0] == 4:
        n_trim = read.cigartuples[-1][1]
        read_sequence = read_sequence[:-n_trim]
        read_qualities = read_qualities[:-n_trim]

    #Check if trimmed read is above minlength threshold
    if args.minlength and len(read_sequence) < args.minlength:
        return(None)

    return(f"@{read_name}\n{read_sequence}\n+\n{read_qualities}\n")

def find_primer_position(args, primer, reference):
    
    def generate_unambiguous_variants(query):
        #Create list of variable and non-variable positions 
        tuples_list = [[i for i in ambiguous_dna_values[j]] for j in query]
        #Create all combinations of all variants
        variants = list(itertools.product(*tuples_list))
        #Combine into strings and create Seq object
        variants = [Seq.Seq(''.join(var)) for var in variants]
        return(variants)

    def create_region_from_alignment(alignment):
        aln_ref, aln_str, aln_primer, _ = str(alignment).split("\n")
        #Find the aligned region
        start, end = re.search("[ACTGRYWSMKHBDVXN-]+", aln_primer).span()
        return(Region(start, end, strand))

    def find_alignment_regions(variants):
        regions = []
        for var in variants:
            for aln in aligner.align(target, var):
                if (aln.score + allowed_mismatch) >= len(query):
                    regions.append(create_region_from_alignment(aln))
            #Also check reverse complement
            for aln in aligner.align(target, var.reverse_complement()):
                if (aln.score + allowed_mismatch) >= len(query):
                    regions.append(create_region_from_alignment(aln))
        return(regions)
    
    query = primer.upper().seq
    target = reference.upper().seq

    if args.fwdkey in primer.id:
        strand = 'left'
    else:
        strand = 'right'

    AMBIGUOUS_DNA_LETTERS = list(ambiguous_dna_values.keys())[4:]

    allowed_mismatch = args.mismatch

    #Setup aligner object
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.match_score = 1
    aligner.mismatch_score = 0
    aligner.open_gap_score = -100 #I Don't allow gaps
    aligner.extend_gap_score = -100 #I Don't allow gaps

    #Check for ambiguous positions in the primer
    if any([nt in AMBIGUOUS_DNA_LETTERS for nt in query]):
        variants = generate_unambiguous_variants(query)
    else:
        variants = [query]

    regions = find_alignment_regions(variants)

    if len(regions) == 0:
        print("Could not find a good match between "+primer.id+" and "+reference.id)
    
    return(regions)

def main():
    args = parser.parse_args()

    primer_records = list(SeqIO.parse(args.primerfile, "fasta"))
    reference = SeqIO.read(args.referencefile, "fasta")

    trim_regions = []

    for primer in primer_records:
        if not ((args.fwdkey in primer.id) | (args.revkey in primer.id)):
            raise ValueError("Neither "+args.fwdkey+" nor "+args.revkey+" could be found in "+primer.id+" which is necessary to determine its orientation")

        for region in find_primer_position(args, primer, reference):
            trim_regions.append(region)

    with pysam.AlignmentFile(args.infile, "rb") as infile, pysam.AlignmentFile(args.outfile, "wb", header=infile.header) as outfile, open(args.outfastq, "w") as outfastq:
        for read in infile.fetch():
            if read.is_unmapped | read.is_secondary | read.is_supplementary:
                continue
            for region in trim_regions:
                
                overlap = calculate_overlap(read, region, padding=args.padding)
                if overlap > 0:
                    try:
                        clip_read(read, overlap, region.strand)
                    except Exception as e:
                        print(read)
                        print(region)
                        raise(e)
            _ = outfile.write(read)
            trimmed_read = trim_read(args, read)
            if trimmed_read:
                _ = outfastq.write(trimmed_read)

if __name__ == "__main__":
    main()