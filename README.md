# Ampliclip
Tool to softclip reads in bam files based on amplicon primers

## Installation

## Installation

1. Clone this repository and ``cd Ampliclip``
2. ``conda env create -f environment.yml --name ampliclip``
3. ``conda activate ampliclip``
4. ``pip install .``

## Full usage

```
usage: ampliclip [-h] -i INFILE -o OUTFILE -fq OUTFASTQ [-l MINLENGTH] -p
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
```