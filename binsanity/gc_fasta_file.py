#!/usr/bin/env python
import os
from Bio import SeqIO
from Bio import GC

def gc_fasta_file(b,name,prefix,location):
    """Makes a tab-delimeted output indicating the GC count associated with each contig"""
    gc_dict = {}
    input_file = open(os.path.join(b,name), 'r') 
    output_file = open(os.path.join(location, str(prefix) + '-GC_count.txt'), 'w')
    output_file.write("%s\t%s\n" % ('contig','GC'))

    for cur_record in SeqIO.parse(input_file, "fasta") :
        gene_name = cur_record.id 
        gc_percent = float(GC(cur_record.seq))
        gc_dict.setdefault(gene_name,[])
        gc_dict[gene_name].append(gc_percent) 
        output_line = '%s\t%i\n' % \
        (gene_name, gc_percent) 
        output_file.write(output_line)

    output_file.close() 
    input_file.close()

