#!/usr/bin/env python
import os
import numpy as np
from Bio import SeqIO

def cov_array(a,path,file_name,size):
    """Computes a coverage array based on the contigs in files and the contig names associated with the coverage file"""
    count = 0
    names = []
    c = 0
    cov_array = []
    for record in SeqIO.parse(os.path.join(path, file_name), "fasta"):
        count = count + 1
        if len(record.seq) >= int(size):
            if record.id in a.keys():
                data = a[record.id]
                names.append(data[0])
                data.remove(data[0])
                line = " ".join(data)
                if c == 1:
                    temparray = np.fromstring(line, dtype=float, sep=' ')
                    cov_array = np.vstack((cov_array, temparray))
                if c == 0:
                    cov_array = np.fromstring(line, dtype=float, sep=' ')
                    c += 1
    print "          %s" % (cov_array.shape,)
    return cov_array, names

