#!/usr/bin/env python
import os
import numpy as np
import pandas as pd

### Binsanity-refine
def get_contigs(c, outdir, kmer):
    gc = open(os.path.join(outdir, "GC_count.txt"), "r")

    dataframe = pd.read_csv(gc, sep="\t", index_col=0, header=None)
    dataframe = dataframe.drop("contig")
    dataframe = dataframe.apply(pd.to_numeric)

    print dataframe.head()

    dataframe *= 100
    dataframe += 1  
    dataframe = np.log10(dataframe)

    tetra = open(os.path.join(outdir, "%smer-frequencies.txt" % (str(kmer))), "r")

    dataframe2 = pd.read_csv(tetra,sep="\t",index_col=0,header=None)
    dataframe2 = dataframe2.drop('contig')
    dataframe2 = dataframe2.apply(pd.to_numeric)
    dataframe2 *= 100
    dataframe2 += 1
    dataframe2 = np.log10(dataframe2)

    cov = open(str(c), "r")

    reader =list(csv.reader(cov,delimiter='\t'))
    titles = ["contig"]
    reader.insert(0, titles)

    dataframe3 = pd.DataFrame(reader)
    dataframe3.columns = dataframe3.iloc[0]
    dataframe3 = dataframe3[1:]
    dataframe3 = dataframe3.set_index("contig")

    result = pd.concat([dataframe, dataframe2, dataframe3], axis=1, join='inner')
    result.to_csv(path_or_buf=str(os.path.join(outdir, "kmer-GC-out.txt")), header=None, sep ="\t")


### Binsanity-wf
def get_contigs(c, prefix, kmer, location):
    gc = open(os.path.join(location, str(prefix) + "_GC_count.txt"), "r")

    dataframe = pd.read_csv(gc,sep="\t")
    dataframe = dataframe.set_index("contig")
    dataframe *= 100
    dataframe += 1
    dataframe = np.log10(dataframe)

    tetra = open(os.path.join(location, str(prefix) + '_%smer_frequencies.txt' % (str(kmer))), "r")

    dataframe2 = pd.read_csv(tetra, sep="\t")
    dataframe2 = dataframe2.set_index("contig")
    dataframe2 *= 100
    dataframe2 += 1
    dataframe2 = np.log10(dataframe2)

    cov = open(str(c),"r")
    reader =list(csv.reader(cov,delimiter='\t'))
    titles = ["contig"]
    reader.insert(0,titles)

    dataframe3 = pd.DataFrame(reader)
    dataframe3.columns = dataframe3.iloc[0]
    dataframe3 = dataframe3[1:]
    dataframe3 = dataframe3.set_index("contig")

    result = pd.concat([dataframe, dataframe2, dataframe3], axis=1, join='inner')
    result.to_csv(path_or_buf=os.path.join(location, str(prefix) + "_kmerGC.txt"), header=None, sep ="\t")


### Binsanity-lc
def get_contigs(c, prefix, kmer, location):
    gc = open(os.path.join(location, str(prefix) + "-GC_count.txt"), "r")
    dataframe = pd.read_csv(gc, sep="\t")
    dataframe = dataframe.set_index("contig")
    dataframe *= 100
    dataframe += 1
    dataframe = np.log10(dataframe)

    tetra = open(os.path.join(location, str(prefix) + "-%smer-frequencies.txt" % (str(kmer))), "r")

    dataframe2 = pandas.read_csv(tetra, sep="\t")
    dataframe2 = dataframe2.set_index("contig")
    dataframe2 *= 100
    dataframe2 += 1
    dataframe2 = np.log10(dataframe2)

    cov = open(str(c),"r")
    reader = list(csv.reader(cov,delimiter='\t'))
    titles = ["contig"]
    reader.insert(0,titles)

    dataframe3 = pd.DataFrame(reader)
    dataframe3.columns = dataframe3.iloc[0]
    dataframe3 = dataframe3[1:]
    dataframe3 = dataframe3.set_index("contig")

    result = pd.concat([dataframe, dataframe2, dataframe3], axis=1, join='inner')
    result.to_csv(path_or_buf=str(os.path.join(location, str(prefix) + "-kmerGC.txt")), header=None, sep ="\t")

