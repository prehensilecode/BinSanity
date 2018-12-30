#!/usr/bin/env python
import sys
import os
import subprocess


### Binsanity-wf
def run_checkM(bins,threads,prefix,directory):
    output = open(os.path.join(directory,str(prefix)+"_checkm_lineagewf-results.txt"),"w")   
    subprocess.call(["checkm","lineage_wf","-x","fna","-t",str(threads),str(bins),os.path.join(directory,str(prefix)+"_binsanity_checkm")],stdout=output)

### Binsanity-lc
def run_checkM(bins,threads,prefix):
    output = open(str(prefix)+"-checkm_lineagewf-binsanity.out","w")   
    subprocess.call(["checkm","lineage_wf","-x","fna","-t",str(threads),'--pplacer_threads',str(threads),str(bins),str(prefix)+"-binsanity-checkm"],stdout=output)
