#!/usr/bin/env python3
'''
Script designed to retrieve CDS seqeunces from a GenBank file format to be used as a database by ToxCodAn-Genome.
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

import sys
from Bio import SeqIO

def _RetrieveCDS_(gb, gene, out):
    final = {}
    for record in SeqIO.parse(gb,"genbank"):
        ID = str(record.id)
        organism = str(record.annotations['organism']).replace(" ","_")
        for feature in record.features:
            if feature.type == "CDS":
                try:
                    proteinID = str(feature.qualifiers["protein_id"][0])
                except:
                    proteinID = gene
                CDS = str(feature.location.extract(record).seq)
                #get only full CDSs
                if CDS[:3] == "ATG" and CDS[-3:] in ["TAA", "TAG", "TGA"] and "N" not in CDS:
                    final[ID] = [organism, proteinID, CDS]
    output = open(out,"w")
    for k in final.keys():
        alfa = final[k]
        output.write(">"+k+"||"+alfa[1]+"||"+alfa[0]+"_"+gene+"\n"+alfa[2]+"\n")
    output.close()

def _main_():

    if len (sys.argv) != 4:
        print("Basic usage: ParseGenBank.py inut.gb output_cds.fasta gene_string")
        print("\t> input.gb: input with CDSs annotation in GENBANK format")
        print("\t> output_cds.fasta: output CDSs in FASTA format")
        print("\t> gene_string: string to determine the gene family at the end of header [\"_gene\"]")
        quit()

    gb = sys.argv[1]
    out = sys.argv[2]
    gene = sys.argv[3]
    _RetrieveCDS_(gb, gene, out)

_main_()

#END
