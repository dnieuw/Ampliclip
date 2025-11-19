import argparse
import pysam
import re
import itertools
import logging
from ampliclip import __version__
from datetime import datetime
from collections import defaultdict
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

parser.add_argument('-log',
                '--logfile',
                help="Name of logfile",
                type=str,
                required = False)

parser.add_argument('-stats',
            '--trimstats',
            help='File with stats about primers and trimming',
            type=str,
            required=False)

parser.add_argument('--quiet',
                    action='store_true',
                    default = False,
                    required = False)

parser.add_argument(
        "--version",
        action="version",
        version=f"%(prog)s {__version__}"
    )

class Region(object):
    """Primer alignment region class"""
    def __init__(self, start, end, strand, primer_id):
        self.start = int(start)
        self.end = int(end)
        self.strand = strand
        self.primer_id = primer_id

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

    #Expand cigarstring - use list for efficiency
    cigar_expanded = []
    for i, j in current_cigar:
        cigar_expanded.extend([i] * j)

    #Clip cigar until no more "clip" left:
    clip_leftover = n_clip
    n = 0
    while clip_leftover > 0:
        cig = cigar_expanded[n]
        if cig == 0:
            cigar_expanded[n] = 4  #Replace matches with softclipped base
            clip_leftover -= 1
        elif cig == 1:
            cigar_expanded[n] = 4  #Replace insertions with softclipped base
        elif cig == 2:
            n_clip += 1  #Increase n_clip to increase reference_start shift at the end in case of deletions
        elif cig == 4:
            pass  #Do not replace softclipped
        elif cig == 5:
            pass  #Do not replace hardclipped
        else:
            raise ValueError("Something went wrong: do not know how to clip " + str(cig) + " in cigarstring")
        n += 1

    #Recreate tuples:
    clipped_cigar = []
    c = 0
    nprev = cigar_expanded[0]
    for n in cigar_expanded:
        if n == nprev:
            c += 1
        else:
            clipped_cigar.append((nprev, c))
            c = 1
            nprev = n
    clipped_cigar.append((n, c))

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
    # Optimize quality score conversion using list comprehension
    read_qualities = ''.join([chr(x + 33) for x in read.query_qualities])

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

# Create a cached aligner object to reuse across all primer alignments
_cached_aligner = None

def get_aligner():
    """Return a cached PairwiseAligner object for reuse"""
    global _cached_aligner
    if _cached_aligner is None:
        _cached_aligner = PairwiseAligner()
        _cached_aligner.mode = 'local'
        _cached_aligner.match_score = 1
        _cached_aligner.mismatch_score = 0
        _cached_aligner.open_gap_score = -100  # Don't allow gaps
        _cached_aligner.extend_gap_score = -100  # Don't allow gaps
    return _cached_aligner

def find_primer_position(args, primer, reference):
    
    def generate_unambiguous_variants(query):
        #Create list of variable and non-variable positions 
        tuples_list = [[i for i in ambiguous_dna_values[j]] for j in query]
        
        # Calculate total number of variants to avoid exponential explosion
        total_variants = 1
        for options in tuples_list:
            total_variants *= len(options)
            # Limit to prevent exponential explosion (e.g., 1000 variants max)
            if total_variants > 1000:
                logging.warning(f"Primer {primer.id} has too many ambiguous bases ({total_variants} variants). "
                               "Using only first 1000 variants for alignment.")
                break
        
        #Create all combinations of all variants (limit to prevent memory issues)
        variants = []
        for i, var in enumerate(itertools.product(*tuples_list)):
            if i >= 1000:  # Stop after 1000 variants
                break
            variants.append(Seq.Seq(''.join(var)))
        return(variants)

    def find_alignment_regions(variants):
        regions = []
        for var in variants:
            for aln in aligner.align(target, var):
                if (aln.score + allowed_mismatch) >= len(query):
                    start, end = aln.coordinates[0]
                    if strand == "right": #Shift by 1 for the right size primer because of 0 indexing
                        regions.append(Region(start+1, end+1, strand, primer.id))
                    else:
                        regions.append(Region(start, end, strand, primer.id))
            #Also check reverse complement
            for aln in aligner.align(target, var.reverse_complement()):
                if (aln.score + allowed_mismatch) >= len(query):
                    start, end = aln.coordinates[0]
                    if strand == "right": #Shift by 1 for the right size primer because of 0 indexing
                        regions.append(Region(start+1, end+1, strand, primer.id))
                    else:
                        regions.append(Region(start, end, strand, primer.id))
        return(regions)
    
    query = primer.upper().seq
    target = reference.upper().seq

    if args.fwdkey in primer.id:
        strand = 'left'
    else:
        strand = 'right'

    AMBIGUOUS_DNA_LETTERS = list(ambiguous_dna_values.keys())[4:]

    allowed_mismatch = args.mismatch

    # Use cached aligner object instead of creating new one
    aligner = get_aligner()

    #Check for ambiguous positions in the primer
    if any([nt in AMBIGUOUS_DNA_LETTERS for nt in query]):
        variants = generate_unambiguous_variants(query)
    else:
        variants = [query]

    regions = find_alignment_regions(variants)
    if len(regions) == 0:
        print("Could not find a good match between "+primer.id+" and "+reference.id)
    
    return(regions)

def log_summary(total_reads_processed, reads_clipped_count, primer_clip_counts, primer_order):

    if total_reads_processed > 0:
        percent_clipped = (reads_clipped_count / total_reads_processed) * 100
    else:
        percent_clipped = 0

    logging.info(f"--- Ampliclip Run Summary ---\nTotal reads processed: {total_reads_processed}\nReads clipped: {reads_clipped_count} ({percent_clipped:.2f}%)")

    if not primer_clip_counts:
        logging.info("""--- Clipping events by primer ---
        No primers caused any clipping events...""")
    else:
        log_text = []
        for primer_id in primer_order:
            count = primer_clip_counts.get(primer_id, 0)
            log_text.append(f"{primer_id:<30} {count}")
        log_text = '\n'.join(log_text)
        logging.info(f"--- Clipping events by primer ---\n{log_text}")

    logging.info("--- End of Summary ---")

def main():
    args = parser.parse_args()

    handlers = []

    # Console handler
    console_handler = logging.StreamHandler()
    if args.quiet:
        console_handler.setLevel(logging.CRITICAL)  # suppress INFO
    else:
        console_handler.setLevel(logging.INFO)
    handlers.append(console_handler)

    # File handler (always INFO, full logs)
    if args.logfile:
        file_handler = logging.FileHandler(args.logfile, mode="w")
        file_handler.setLevel(logging.INFO)
        handlers.append(file_handler)

    logging.basicConfig(level=logging.DEBUG,
                            handlers=handlers,
                            format='%(levelname)s %(asctime)s %(message)s',
                            datefmt='%Y-%m-%d %H:%M:%S')

    primer_records = list(SeqIO.parse(args.primerfile, "fasta"))
    reference = SeqIO.read(args.referencefile, "fasta")

    primer_order = []
    for primer in primer_records:
        primer_order.append(primer.id)

    trim_regions = []

    logging.info("Aligning primers to reference...")
    
    for primer in primer_records:
        if not ((args.fwdkey in primer.id) | (args.revkey in primer.id)):
            raise ValueError("Neither "+args.fwdkey+" nor "+args.revkey+" could be found in "+primer.id+" which is necessary to determine its orientation")

        for region in find_primer_position(args, primer, reference):
            trim_regions.append(region)

    # Initialize counters
    total_reads_processed = 0
    reads_clipped_count = 0
    primer_clip_counts = defaultdict(int)

    logging.info("Trimming primers...")

    with pysam.AlignmentFile(args.infile, "rb") as infile, pysam.AlignmentFile(args.outfile, "wb", header=infile.header) as outfile, open(args.outfastq, "w") as outfastq:
        for read in infile.fetch():
            if read.is_unmapped | read.is_secondary | read.is_supplementary:
                continue
            
            # Increment total read counter and set up a flag for this read
            total_reads_processed += 1
            read_was_clipped = False

            # Check if read overlaps with any primer region before detailed processing
            # Quick boundary check to skip reads that don't overlap any region
            if trim_regions:
                read_start = read.reference_start
                read_end = read.reference_end
                has_potential_overlap = False
                
                for region in trim_regions:
                    # Quick range check with padding
                    if not (read_end < region.start - args.padding or read_start > region.end + args.padding):
                        has_potential_overlap = True
                        break
                
                if has_potential_overlap:
                    for region in trim_regions:
                        overlap = calculate_overlap(read, region, padding=args.padding)
                        if overlap > 0:
                            try:
                                clip_read(read, overlap, region.strand)
                                # If a clip happened, set flag and count the primer
                                read_was_clipped = True
                                primer_clip_counts[region.primer_id] += 1
                            except Exception as e:
                                print(read)
                                print(region)
                                raise(e)
            # Increment the clipped read counter if the flag was set
            if read_was_clipped:
                reads_clipped_count += 1

            _ = outfile.write(read)
            trimmed_read = trim_read(args, read)
            if trimmed_read:
                _ = outfastq.write(trimmed_read)

    log_summary(total_reads_processed, reads_clipped_count, primer_clip_counts, primer_order)

if __name__ == "__main__":
    main()