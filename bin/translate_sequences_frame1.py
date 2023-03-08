#!/usr/bin/env python3
'''
script to translate coding sequences to peptide - using frame 1 only
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

import sys
from Bio import SeqIO
from Bio.Seq import Seq

def _Translate_(IN,OUT):
    translated = open(OUT,"w")
    for record in SeqIO.parse(IN, "fasta"):
        ID = str(record.id)
        protein = str(record.seq.translate())
        translated.write(">"+ID+"\n"+protein+"\n")
    translated.close()

def _main_():
    if len (sys.argv) != 3:
        print("Basic usage: translate_sequences_frame1.py input_fasta translated_output_fasta")
        print("\t> input_fasta: curated toxins in fasta format")
        print("\t> translated_output_fasta: translated toxins in fasta format")
        quit()

    fasta = sys.argv[1]
    out = sys.argv[2]
    _Translate_(fasta,out)

_main_()

#END
