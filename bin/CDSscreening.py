#!/usr/bin/env python3
'''
CDSscreening - Survey for CDSs in Transcriptome assembly using a specific database
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

##Modules
import os
import subprocess as sp
import warnings
from optparse import OptionParser
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio import BiopythonWarning
    warnings.simplefilter('ignore', BiopythonWarning)
except:
    print('''ToxCodAn was not able to run due to the ERROR below:
    The biopython package is not properly installed or it is not active in your actual environment.
    Please, install biopython properly (check https://biopython.org/ for more details), or active the environment where biopython is installed!''')
    quit()

##Functions
#parse fasta file
def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta,"fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

#translate cds
def _Translate_(seq):
    coding_dna = Seq(seq)
    translated = str(coding_dna.translate())
    return translated

#parse blast search
def _ParseBlast_(blast, fasta):
    final = {}
    processed = []
    F = _ParseFasta_(fasta)
    a = open(blast,"r")
    for line in a:
        line1 = line.strip().split("\t")
        qID = str(line1[0])
        qLEN = int(line1[1])
        percid = float(line1[4])
        seqID = str(line1[2])
        st = int(line1[9])
        end = int(line1[10])
        SEQ = F[seqID]
        CDS = ""
        if seqID not in processed:
            #hit at plus strand
            if st < end:
                if SEQ[st-1:st+2] == "ATG":
                    CDS += SEQ[st-1:st+2]
                    for n in range(st+3, end+66, 3):
                        codon = SEQ[n-1:n+2]
                        CDS += codon
                        if SEQ[n-1:n+2] in ["TAA", "TAG", "TGA"]:
                            break
                    if CDS[-3:] in ["TAA", "TAG", "TGA"]:
                        PEP = _Translate_(CDS)
                        final.setdefault(seqID,[qID, str(percid), CDS, PEP])
                        processed.append(seqID)
            #hit at minus strand
            if st > end:
                #adjust positions and reverse complement sequence
                newst = len(SEQ)-st
                newend = newst+(st-end)
                newSEQ = str(Seq(SEQ).reverse_complement())
                if newSEQ[newst:newst+3] == "ATG":
                    CDS += newSEQ[newst:newst+3]
                    for n in range(newst+3, newend+66, 3):
                        codon = newSEQ[n:n+3]
                        CDS += codon
                        if newSEQ[n:n+3] in ["TAA", "TAG", "TGA"]:
                            break
                    if CDS[-3:] in ["TAA", "TAG", "TGA"]:
                        PEP = _Translate_(CDS)
                        final.setdefault(seqID,[qID, str(percid), CDS, PEP])
                        processed.append(seqID)
    a.close()
    return final

def _RunScreening_(transcripts, database, outF, percid, length, cpu):
    #blast search
    sp.check_output("makeblastdb -dbtype nucl -in "+transcripts+" -out "+outF+"blastDB/dna", shell=True)
    sp.call("blastn -num_threads "+str(cpu)+" -query "+database+" -db "+outF+"blastDB/dna -out "+outF+"temp_blast.out -perc_identity "+str(percid)+" -qcov_hsp_perc 0.95 -outfmt \'6 qseqid qlen sseqid slen pident length mismatch qstart qend sstart send evalue bitscore\'", shell=True)
    os.system("sort -k5 -nr "+outF+"temp_blast.out > "+outF+"blast.out")
    os.system("rm "+outF+"temp_blast.out")
    blast = outF+"blast.out"
    B = _ParseBlast_(blast, transcripts)
    CDS = open(outF+"cds_screening.fasta","w")
    PEP = open(outF+"pep_screening.fasta","w")
    for k in B.keys():
        if len(B[k][2]) >= length:
            CDS.write(">"+k+"||"+B[k][1]+"||"+str(B[k][0])+"\n"+B[k][2]+"\n")
            PEP.write(">"+k+"||"+B[k][1]+"||"+str(B[k][0])+"\n"+B[k][3]+"\n")
    PEP.close()
    CDS.close()
    os.system("rm -r "+outF+"blastDB/")

##Options
def __main__():
    parser = OptionParser()
    parser.add_option("-t", "--transcripts", dest="transcripts", help="Mandatory - transcripts recovered from RNA-seq data for the individual/species using de-novo/genome-guided assembly methods in FASTA format, /path/to/transcripts.fasta", metavar="fasta", default=None)
    parser.add_option("-d", "--database", dest="database", help="Mandatory - database with coding sequences (CDSs) in FASTA format, /path/to/cds.fasta", metavar="fasta", default=None)
    parser.add_option("-o", "--output", dest="output", help="Optional - output folder, /path/to/output_folder; if not defined, the output folder will be set in the current directory with the following name [screening]", metavar="folder", default=None)
    parser.add_option("-p", "--percid", dest="percid", help="Optional - threshold value used as the minimum percent identity between CDS in the database and transcript [default=80]", metavar="int", default="80")
    parser.add_option("-l", "--length", dest="length", help="Optional - minimum size of a CDS; it will remove any identified CDS shorter than the specified threshold [default=200]", metavar="int", default="200")
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used in each step [default=1]", metavar="int", default="1")

    (options, args) = parser.parse_args()

    if options.transcripts == None or options.database == None:
        print(
        """
>>>> Screening transcripts using a CDS database <<<<
      ****Use -h for help!****

USAGE:
CDSscreening.py -t transcripts.fasta -d CDS_database.fasta
        """)
        quit()

    if options.output != None:
        if not options.output.endswith("/"):
            options.output += "/"
    if options.output == None:
        CWD = os.getcwd()
        options.output = CWD+"/screening/"

    if os.path.isdir(options.output) == False:
        os.mkdir(options.output)

    if options.transcripts != None and options.database != None:

        _RunScreening_(options.transcripts,
                    options.database,
                    options.output,
                    str(options.percid),
                    int(options.length),
                    str(options.cpu))

if __name__ == '__main__':
    __main__()

#END