#!/usr/bin/env python3
'''
Script designed to adjust fgenesh+ output from a specific region of the genome to fit into the main genome.
 - It is a built-in script from ToxCodAn-Genome
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

import sys

def _Run_(fgene, geneid, start):
    chr = start.split("-")[0]
    start = int(start.split("-")[1])
    a = open(fgene,"r")
    for line in a:
        if "CDSi" in line or "CDSf" in line or "CDSl" in line:
            line1 = line.strip().split()
            #print(line1)
            strand = line1[1]
            st = int(line1[4])
            end = int(line1[6])
            print(chr+"\tfgenesh\tCDS\t"+str(start+st-1)+"\t"+str(start+end-1)+"\t.\t"+strand+"\t.\tgene_id \""+geneid+"\"; transcript_id \""+geneid+"\";")
    a.close()

def _main_():
    if len (sys.argv) != 4:
        print("Basic usage: AdjustingFgenesh.py fgenesh.txt geneid ch-start")
        print("\t> fgenesh.txt: output of fgenesh+ in txt format")
        print("\t> geneid: id of the gene desired in the output")
        print("\t> ch-start: chromosome/contig and starting nucleotide of the region used to annotate the gene within fgenesh+ (e.g., contig_1-1988)")
        quit()

    fgene = sys.argv[1]
    geneid = sys.argv[2]
    start = sys.argv[3]
    _Run_(fgene, geneid, start)

_main_()

#END
