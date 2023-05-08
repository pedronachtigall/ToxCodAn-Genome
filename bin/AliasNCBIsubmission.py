#!/usr/bin/env python3
'''
Script designed to convert gene symbols into gene names accepted in NCBI. It was designed to help prepare the annotation files for submission to NCBI.
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

import os
import sys
from optparse import OptionParser

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

def _AliasDefault_():
    final = {"SVMP":"snake venom metalloproteinase",
    "SVSP":"snake venom serine protease",
    "PLA2":"phospolipase A2",
    "CTL":"C-type lectin",
    "VEGF":"vascular endothelial growth factor",
    "VEGFF":"vascular endothelial growth factor",
    "HYAL":"hyaluronidase",
    "NGF":"nerve growth factor",
    "NUC":"nucleotidase",
    "SVLIPA":"snake venom lipase",
    "Vespryn":"vespryn",
    "Ohanin":"Vespryn",
    "Waprin":"snake venom waprin",
    "KUN":"snake venom kunitz",
    "3FTx":"three-finger toxin",
    "3FTX":"three-finger toxin",
    "BPP":"bradykinin-potentiating peptide",
    "CNP":"C-type natriuretic peptide",
    "CRISP":"cysteine-rich secretory protein",
    "CVF":"venom factor",
    "VF":"venom factor"}
    return final

def _ParseAlias_(alias):
    final = {}
    a = open(alias,"r")
    for line in a:
        line1 = line.strip().split("\t")
        final[line1[0]] = line1[1]
    a.close()
    return final

## add options
 #design a default database alias, but also accepts a file designed by the user

def _Alias_(gtf, alias, output):
    if alias != None:
        A = _ParseAlias_(alias)
    if alias == None:
        A = _AliasDefault_()

    NotAlias = set([])

    a = open(gtf,"r")
    OUT = open(output,"w")
    for line in a:
        line1 = line.strip().split("\t")
        if line1[2] == "CDS":
            id = line1[-1].split(";")[0].replace("\"","").replace("ID=","").replace("transcript_id","").replace("gene_id","").replace(";","").strip()
            if "-" in id:
                gs = id.split("-")[0]
            if "-" not in id:
                gs = id
            if gs in A.keys():
                OUT.write(line.strip()+" product \""+A[gs]+"\"\n")
            if not gs in A.keys():
                NotAlias.add(id)
                OUT.write(line.strip()+" product \"hypothetical protein\"\n")
    a.close()
    OUT.close()

    if len(NotAlias) >= 1:
        print("The following entries have no gene descrition in the alias used:")
        for i in NotAlias:
            print(i)
        print("\nConsider complementing your alias file and re-run this script.")
        print("See https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide for further details.")
        print("Otherwise, all of the entries above will have the \"product\" feature marked as \"hypothetical protein\".")

##>>>>Options
def __main__():
    parser = OptionParser()
    parser.add_option("-g", "--gtf", dest="gtf", help="Mandatory - GTF file to be used as template for the output. [e.g., input.gtf]", metavar="gtf", default=None)
    parser.add_option("-a", "--alias", dest="alias", help="Optional - user-supplied file in TXT format indicating an alias to add a \"product=\" feature in the description column of GTF for submission to NCBI. The alias must be a two-column TXT file with gene symbol in the first column and gene description in the second column. Of note, the built-in alias was designed for snake toxins, consider designing your own alias file to complement the final annotation file. [default=None]", metavar="file", default=None)
    parser.add_option("-o", "--output", dest="output", help="Optional - output file in GTF format. If not specified it will be output in the same directory as the input GTF file with a suffix \"_NCBI\". [default=input_NCBI.gtf]", metavar="file", default=None)

    (options, args) = parser.parse_args()

    if options.gtf == None:
        print("""
BASIC USAGE:
AliasNCBIsubmission.py -g input.gtf -a alias.txt -o output.gtf

Use -h for help!

        """)
        quit()

    if options.output == None:
        options.output = options.gtf.split(".")[0]+"_NCBI.gtf"

    if options.gtf != None:
        _Alias_(options.gtf, options.alias, options.output)
        print("\n\n>>> Finished <<<")

#__main__()

_Alias_("/home/nachtigall/Desktop/posdoc_butantan/toxcodan-genome/blast_uniprot/ToxCodAnGenome_output/toxin_annotation.gtf", None, "/home/nachtigall/Desktop/posdoc_butantan/toxcodan-genome/blast_uniprot/ToxCodAnGenome_output/toxin_annotation_NCBI.gtf")

#END
