#!/usr/bin/env python
import collections
from Bio import SeqIO


### Binsanity-refine
def kmer_counts(kmer, b, name):
    kmer_dict = {}
    for record in SeqIO.parse(os.path.join(b, name), "fasta"):
        id_ = str(record.id)
        seq = str(record.seq)
        length = len(seq)
        tetra = kmer_list(seq, kmer)
        c = collections.Counter(tetra)
        c = dict(c)
        val = list(c.values())
        val_edit = []
        for freq in val:
            freq= (float(freq)/float(length))
            val_edit.append(freq)

        keys = []
        for key in c:
            keys.append(str(key))
        c_edit = dict(zip(keys, val_edit))
        kmer_dict.setdefault(str(id_), [])
        kmer_dict[str(id_)].append(c_edit)
    return kmer_dict

### Binsanity-wf
def kmer_counts(kmer,b,name):
    kmer_dict = {}
    for record in SeqIO.parse(os.path.join(b, name), "fasta"):
        id_ = record.id
        seq = record.seq
        length = len(seq)
        tetra = kmer_list(seq, kmer)
        c = collections.Counter(tetra)
        c = dict(c)
        val = list(c.values())
        val_edit = []
        for freq in val:
            freq= (float(freq)/float(length))
            val_edit.append(freq)

        keys = []
        for key in c:
            keys.append(str(key))
        c_edit = dict(zip(keys,val_edit))
        kmer_dict.setdefault(str(id_), [])
        kmer_dict[str(id_)].append(c_edit)
    return kmer_dict

### Binsanity-lc
def kmer_counts(kmer, b, name):
    kmer_dict = {}
    for record in SeqIO.parse(os.path.join(b,name), "fasta"):
        id_ = record.id
        seq = record.seq
        length = len(seq)
        tetra = kmer_list(seq, kmer)
        c = collections.Counter(tetra)
        c = dict(c)
        val = list(c.values())
        val_edit = []
        for freq in val:
            freq = (float(freq)/float(length))
            val_edit.append(freq)

        keys = []
        for key in c:
            keys.append(str(key))
        c_edit = dict(zip(keys,val_edit))
        kmer_dict.setdefault(str(id_),[])
        kmer_dict[str(id_)].append(c_edit)
    return kmer_dict

