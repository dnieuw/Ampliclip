import argparse
import pysam
import re
import itertools
import logging
import os
import sqlite3
import tempfile
from contextlib import ExitStack, closing
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
            required=False)

parser.add_argument('-fq1', '--outfastq1', type=str,
                    help='Output trimmed R1 FASTQ (requires -fq2; excludes -fq)')
parser.add_argument('-fq2', '--outfastq2', type=str,
                    help='Output trimmed R2 FASTQ (requires -fq1; excludes -fq)')

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

parser.add_argument('-t', '--threads', help='Number of threads for BAM I/O', default=1, type=int)

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
    # If all aligned bases are clipped completely clip cigarstring
    if read.qlen <= n_clip:
        read.cigartuples = [(4, read.query_length)]
        return

    current_cigar = read.cigartuples
    if side == "right":
        current_cigar.reverse()

    clipped_cigar = []
    clip_leftover = n_clip
    softclips_to_add = 0
    ref_shift = 0  # To track deletions for reference_start shifting

    for op, length in current_cigar:
        if clip_leftover > 0:
            if op == 4:  # Soft clip (Does not consume reference)
                softclips_to_add += length
            elif op == 5:  # Hard clip
                if softclips_to_add > 0:
                    clipped_cigar.append((4, softclips_to_add))
                    softclips_to_add = 0
                clipped_cigar.append((op, length))
            elif op == 1:  # Insertion (Does not consume reference)
                softclips_to_add += length
            elif op == 2:  # Deletion (Consumes reference, but not query)
                ref_shift += length
            elif op == 0:  # Match
                if length <= clip_leftover:
                    softclips_to_add += length
                    clip_leftover -= length
                else:
                    softclips_to_add += clip_leftover
                    clipped_cigar.append((4, softclips_to_add))
                    softclips_to_add = 0
                    
                    clipped_cigar.append((0, length - clip_leftover))
                    clip_leftover = 0
            else:
                raise ValueError("Something went wrong: do not know how to clip " + str(op) + " in cigarstring")
        else:
            # clip_leftover is exhausted, flush any pending softclips immediately
            if softclips_to_add > 0:
                clipped_cigar.append((4, softclips_to_add))
                softclips_to_add = 0
            clipped_cigar.append((op, length))

    # Flush any lingering softclips at the very end
    if softclips_to_add > 0:
        clipped_cigar.append((4, softclips_to_add))

    # Merge consecutive identical operations to keep CIGAR valid (e.g., 5S followed by 10S becomes 15S)
    merged_cigar = []
    for op, length in clipped_cigar:
        if length == 0:
            continue
        if merged_cigar and merged_cigar[-1][0] == op:
            merged_cigar[-1] = (op, merged_cigar[-1][1] + length)
        else:
            merged_cigar.append((op, length))

    if side == 'right':
        merged_cigar.reverse()

    read.cigartuples = merged_cigar

    if side == 'left':
        # Shift reference start by the number of clipped bases + any deletions we passed
        read.reference_start = read.reference_start + n_clip + ref_shift

def trim_read(args, read, paired=False):
    '''
    Trims softclipped nucleotides of the reads from the left and right side and returns read in fastq format.
    '''

    if paired:
        # BAM stores reverse-strand reads in reference orientation. Slice in
        # that orientation first, then restore both sequence and quality order.
        if read.query_sequence is None or read.query_qualities is None:
            return None
        start, end = read.query_alignment_start, read.query_alignment_end
        sequence = read.query_sequence[start:end]
        qualities = read.query_qualities[start:end]
        if not sequence or (args.minlength and len(sequence) < args.minlength):
            return None
        if read.is_reverse:
            sequence = str(Seq.Seq(sequence).reverse_complement())
            qualities = qualities[::-1]
        quality_string = ''.join(chr(q + 33) for q in qualities)
        return f"@{read.query_name}\n{sequence}\n+\n{quality_string}\n"

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

def log_summary(total_reads_processed, reads_clipped_count, primer_clip_counts, primer_order, primer_coordinates):

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
        # Add header for the columns
        log_text.append(f"{'Primer ID':<30} {'Count':<10} {'Position':<15}")
        
        for primer_id in primer_order:
            count = primer_clip_counts.get(primer_id, 0)
            # Retrieve coordinates, using set to remove duplicates if variants mapped to same location
            coords = list(set(primer_coordinates.get(primer_id, ["N/A"])))
            # Sort for consistency
            coords.sort() 
            coords_str = ",".join(coords)
            
            log_text.append(f"{primer_id:<30} {count:<10} {coords_str}")
            
        log_text = '\n'.join(log_text)
        logging.info(f"--- Clipping events by primer ---\n{log_text}")

    logging.info("--- End of Summary ---")

def parse_args(argv=None):
    args = parser.parse_args(argv)
    paired = args.outfastq1 is not None or args.outfastq2 is not None
    if paired:
        if args.outfastq:
            parser.error('-fq cannot be combined with -fq1 or -fq2')
        if not args.outfastq1 or not args.outfastq2:
            parser.error('paired-end output requires both -fq1 and -fq2')
    elif not args.outfastq:
        parser.error('provide -fq, or both -fq1 and -fq2')
    # Reject aliases before opening anything for writing (including symlinks
    # and existing hard links).
    if paired:
        inputs = [args.infile, args.primerfile, args.referencefile]
        outputs = [args.outfile, args.outfastq1, args.outfastq2]
        if args.logfile:
            outputs.append(args.logfile)
        def same_file(a, b):
            return (os.path.realpath(a) == os.path.realpath(b) or
                    (os.path.exists(a) and os.path.exists(b) and os.path.samefile(a, b)))
        for index, output in enumerate(outputs):
            if any(same_file(output, other) for other in inputs + outputs[:index]):
                parser.error('output paths must be distinct and must not overwrite inputs')
    return args


class PairedFastqWriter:
    """Disk-backed mate matching, independent of BAM sort order.

    Keep even rejected records until finalization so duplicate primary mate
    records are detected rather than paired arbitrarily.
    """
    def __init__(self, connection, r1, r2):
        self.db = connection
        self.r1, self.r2 = r1, r2
        self.skipped = 0
        self.db.execute('CREATE TABLE mates (name TEXT, rg TEXT, mate INTEGER, '
                        'fastq TEXT, PRIMARY KEY (name, rg, mate))')

    def add(self, read, fastq):
        if not read.is_paired or read.is_read1 == read.is_read2 or not read.query_name:
            self.skipped += 1
            return
        mate = 1 if read.is_read1 else 2
        rg = read.get_tag('RG') if read.has_tag('RG') else ''
        try:
            self.db.execute('INSERT INTO mates VALUES (?, ?, ?, ?)',
                            (read.query_name, rg, mate, fastq))
        except sqlite3.IntegrityError as exc:
            raise ValueError('Duplicate primary R{} for read {!r} (RG {!r})'
                             .format(mate, read.query_name, rg)) from exc

    def finish(self):
        self.db.commit()
        pairs = 0
        for r1, r2 in self.db.execute(
                'SELECT a.fastq, b.fastq FROM mates a JOIN mates b '
                'ON a.name = b.name AND a.rg = b.rg '
                'WHERE a.mate = 1 AND b.mate = 2 '
                'AND a.fastq IS NOT NULL AND b.fastq IS NOT NULL ORDER BY a.rowid'):
            self.r1.write(r1)
            self.r2.write(r2)
            pairs += 1
        templates = self.db.execute(
            'SELECT COUNT(*) FROM (SELECT name, rg FROM mates GROUP BY name, rg)'
        ).fetchone()[0]
        logging.info('Paired FASTQ: %s pairs written; %s incomplete or filtered '
                     'templates omitted; %s records with invalid/missing mate flags skipped',
                     pairs, templates - pairs, self.skipped)


def main():
    args = parse_args()

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
    # Dictionary to store coordinates for the log summary
    primer_coordinates = defaultdict(list)

    logging.info("Aligning primers to reference...")
    
    for primer in primer_records:
        if not ((args.fwdkey in primer.id) | (args.revkey in primer.id)):
            raise ValueError("Neither "+args.fwdkey+" nor "+args.revkey+" could be found in "+primer.id+" which is necessary to determine its orientation")

        for region in find_primer_position(args, primer, reference):
            trim_regions.append(region)
            # Store the coordinate range found for this primer
            primer_coordinates[region.primer_id].append(f"{region.start}-{region.end}")

    # Initialize counters
    total_reads_processed = 0
    reads_clipped_count = 0
    primer_clip_counts = defaultdict(int)

    logging.info("Trimming primers...")

    paired = args.outfastq1 is not None
    with ExitStack() as stack:
        infile = stack.enter_context(pysam.AlignmentFile(args.infile, "rb", threads=args.threads))
        outfile = stack.enter_context(pysam.AlignmentFile(
            args.outfile, "wb", header=infile.header, threads=args.threads))
        if paired:
            r1 = stack.enter_context(open(args.outfastq1, 'w'))
            r2 = stack.enter_context(open(args.outfastq2, 'w'))
            tempdir = stack.enter_context(tempfile.TemporaryDirectory(prefix='ampliclip-'))
            db = stack.enter_context(closing(sqlite3.connect(os.path.join(tempdir, 'mates.db'))))
            pair_writer = PairedFastqWriter(db, r1, r2)
        else:
            outfastq = stack.enter_context(open(args.outfastq, 'w'))
        for read in infile:
            if read.is_secondary or read.is_supplementary:
                continue
            if read.is_unmapped:
                if paired:
                    pair_writer.add(read, None)
                continue
            
            # Increment total read counter and set up a flag for this read
            total_reads_processed += 1
            read_was_clipped = False

            for region in trim_regions:
                if read.reference_end < region.start - args.padding or read.reference_start > region.end + args.padding:
                    continue
                
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
            trimmed_read = trim_read(args, read, paired=paired)
            if paired:
                pair_writer.add(read, trimmed_read)
            elif trimmed_read:
                _ = outfastq.write(trimmed_read)

        if paired:
            pair_writer.finish()

    log_summary(total_reads_processed, reads_clipped_count, primer_clip_counts, primer_order, primer_coordinates)

if __name__ == "__main__":
    main()
