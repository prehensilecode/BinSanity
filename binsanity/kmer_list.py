#!/usr/bin/env python
import os
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

### Binsanity-refine
def kmer_list(dna, k):
    """Makes list of k-mers based on input of k and the dna sequence"""
    result = []
    dna = dna.upper()
    dna_edit = Seq(dna, IUPAC.unambiguous_dna)
    reverse_complement = (dna_edit).reverse_complement()
    dna_edit = str(dna_edit)
    reverse_complement = str(reverse_complement)    

    for x in range(len(dna_edit) + 1 - k):
        result.append(dna_edit[x:x+k])

    for x in range(len(reverse_complement)+1-k):
        result.append(reverse_complement[x:x+k])

    result = [s for s in result if not s.strip('AGTC')]
    return result

### Binsanity-wf
def kmer_list(dna, k):
    """Makes list of k-mers based on input of k and the dna sequence"""
    result = []
    dna = dna.upper()   
    dna_edit = Seq(str(dna), IUPAC.unambiguous_dna)
    reverse_complement = (dna_edit).reverse_complement()
    dna_edit = str(dna_edit)
    reverse_complement = str(reverse_complement)

    for x in range(len(dna)+1-k):
        result.append(dna[x:x+k])

    for x in range(len(reverse_complement)+1-k):
        result.append(reverse_complement[x:x+k])

    result = [s for s in result if not s.strip('AGTC')]
    return result


### Binsanity-lc
def kmer_list(dna, k):
    """Makes list of k-mers based on input of k and the dna sequence"""
    result = []
    dna = dna.upper()
    dna_edit = Seq(str(dna))
    reverse_complement = (dna_edit).reverse_complement()

    for x in range(len(dna)+1-k):
        result.append(dna[x:x+k])

    for x in range(len(reverse_complement)+1-k):
        result.append(reverse_complement[x:x+k])

    result = [s for s in result if not s.strip('AGTC')]
    return result

