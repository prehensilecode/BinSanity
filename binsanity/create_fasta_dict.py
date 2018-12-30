#!/usr/bin/env python
from Bio import SeqIO

def create_fasta_dict(fastafile):
    """makes dictionary using fasta files to be binned"""
    fasta_dict = {}
    for record in SeqIO.parse(fastafile, "fasta"):
        fasta_dict[record.id] = record.seq
    return fasta_dict
