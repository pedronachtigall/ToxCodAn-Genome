#!/usr/bin/env python3
'''
Script designed to convert a GTF format file output by ToxCodAn-Genome (which contains CDSs of toxin genes) into a tab-delimited file to be used to plot toxin loci using R and gggenes (or other package and softwares).
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

import sys

def _ParseGTF_(gtf, out):
    ch = {}
    start = {}
    end = {}
    strand = {}
    a = open(gtf,"r")
    for line in a:
        if not line.startswith("#") and not line.startswith("\n"):
            line1 = line.strip().split("\t")
            id = line1[-1].split(";")[0].replace("\"","").replace("ID=","").replace("transcript_id","").replace("gene_id","").replace(";","").strip()
            ch[id] = line1[0]
            start.setdefault(id,[])
            start[id].append(line1[3])
            end.setdefault(id,[])
            end[id].append(line1[4])
            strand[id] = line1[6]
    a.close()
    OUT = open(out, "w")
    for id in ch.keys():
        fam = id.split("-")[0]
        OUT.write(id+"\t"+fam+"\t"+ch[id]+"\t"+min(start[id])+"\t"+max(end[id])+"\t"+strand[id]+"\n")
    OUT.close()

def _main_():
    if len (sys.argv) != 3:
        print("Basic usage: fromCDStoGENE.py toxin_annotation.gtf toxin_annotation_GENE.tsv")
        quit()

    gtf = sys.argv[1]
    out = sys.argv[2]
    _ParseGTF_(gtf, out)

_main_()

#END