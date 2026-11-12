# Ampliclip
Tool to softclip single or paired reads in bam files based on amplicon primers

## Installation

1. Clone this repository and ``cd Ampliclip``
2. ``conda env create -f environment.yml --name ampliclip``
3. ``conda activate ampliclip``
4. ``pip install .``

## Full usage

```
usage: ampliclip [-h] -i INFILE -o OUTFILE [-fq OUTFASTQ]
                    [-fq1 OUTFASTQ1] [-fq2 OUTFASTQ2] [-l MINLENGTH] -p
                    PRIMERFILE -r REFERENCEFILE -fwd FWDKEY -rev REVKEY
                    [-x PADDING] [-m MISMATCH]

optional arguments:
  -h, --help            show this help message and exit
  -i INFILE, --infile INFILE
                        Input bamfile
  -o OUTFILE, --outfile OUTFILE
                        Output bamfile
  -fq OUTFASTQ, --outfastq OUTFASTQ
                        Output fastq file for trimmed reads
  -fq1 OUTFASTQ1, --outfastq1 OUTFASTQ1
                        Output trimmed R1 FASTQ (requires -fq2; excludes -fq)
  -fq2 OUTFASTQ2, --outfastq2 OUTFASTQ2
                        Output trimmed R2 FASTQ (requires -fq1; excludes -fq)
  -l MINLENGTH, --minlength MINLENGTH
                        Minimum length of fastq read after trimming allowed
  -p PRIMERFILE, --primerfile PRIMERFILE
                        Input primer file in fasta format
  -r REFERENCEFILE, --referencefile REFERENCEFILE
                        Input reference file in fasta format
  -fwd FWDKEY, --fwdkey FWDKEY
                        Keyword that indicates forward primer
  -rev REVKEY, --revkey REVKEY
                        Keyword that indicates reverse primer
  -x PADDING, --padding PADDING
                        Number of nucleotides allowed to be upstream of the
                        primer while clipping
  -m MISMATCH, --mismatch MISMATCH
                        Number of mismatches allowed while matching and
                        determining the position of the primers
  -t THREADS, --threads
                        Number of threads for BAM I/O
```

Use either `-fq OUTFASTQ` (single-read mode), or **both**
`-fq1 OUTFASTQ1` and `-fq2 OUTFASTQ2` (paired-end mode). Mixing modes or supplying
only one mate output is an error.
