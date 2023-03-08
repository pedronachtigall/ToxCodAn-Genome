#!/usr/bin/env python3
'''
ToxCodAn-Genome - Toxin Coding Sequence Annotation in Genomic Data
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

##Modules
import os
import datetime as dt
import subprocess as sp
import csv
import warnings
from optparse import OptionParser
try:
    import pandas as pd
    pd.set_option('mode.chained_assignment', None)
except:
    print('''ToxCodAn-Genome was not able to run due to the ERROR below:
    The pandas package is not properly installed or it is not active in your actual environment.
    Please, install pandas properly (check https://pandas.pydata.org/ for more details), or active the environment where pandas is installed!''')
try:
    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio import BiopythonWarning
    warnings.simplefilter('ignore', BiopythonWarning)
except:
    print('''ToxCodAn-Genome was not able to run due to the ERROR below:
    The biopython package is not properly installed or it is not active in your actual environment.
    Please, install biopython properly (check https://biopython.org/ for more details), or active the environment where biopython is installed!''')
    quit()


##Functions
#parse fasta file
def _ParseFasta_(fasta):
    final = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ID = str(record.id)
        SEQ = str(record.seq)
        final[ID] = SEQ
    return final

#parse GTF file
def _ParseGTF_(gtf):
    G = {}
    a = open(gtf, "r")
    for line in a:
        line1 = line.strip().split("\t")
        id = line1[-1].split("; ")[0].replace("gene_id ","").replace("\"","")
        G.setdefault(id, [])
        G[id].append(line1)
    a.close()
    return G

#translate CDS considering frame 1
#write warning report about sequences that need to be checked with caution
def _TranslateCDS_(outF):
    WARNING = open(outF+"annotation_warning.txt","w")
    translated = open(outF+"toxin_annotation_pep.fasta","w")
    for record in SeqIO.parse(outF+"toxin_annotation_cds.fasta", "fasta"):
        ID = str(record.id)
        protein = str(record.seq.translate())
        translated.write(">"+ID+"\n"+protein+"\n")
        if protein.count("*") > 1:
            WARNING.write(ID+"\tThis annotation needs to be manually checked. It may represents a truncated paralog, a pseudogene or an erroneous annotation.\n")
        if protein.count("*") == 1 and protein[-1] != "*":
            WARNING.write(ID+"\tThis annotation needs to be manually checked. It may represents a truncated paralog, a pseudogene or an erroneous annotation.\n")
    translated.close()
    WARNING.close()

#remove partial annotations
def _CheckPartials_(gtf, fasta, Fout, CDSlength):
    F = _ParseFasta_(fasta)
    G = _ParseGTF_(gtf)
    ToxLociFinal = {}
    ToxLociN = {}
    ToxLociPartial = {}
    for k in F.keys():
        #correct codon frame and CDS longer than threshold
        if len(F[k]) % 3 == 0 and len(F[k]) >= CDSlength:
            toxfam = k.split("-")[0]
            ToxLociN.setdefault(toxfam, 0)
            ToxLociN[toxfam] += 1
            geneid = toxfam+"-"+str(ToxLociN[toxfam])
            ToxLociFinal.setdefault(geneid, [])
            for line in G[k]:
                line[-1] = "gene_id \""+geneid+"\"; transcript_id \""+geneid+"\";"
                ToxLociFinal[geneid].append("\t".join(line))
        #correct codon frame but CDS shorter than threshold
        if len(F[k]) % 3 == 0 and len(F[k]) < CDSlength:
            for line in G[k]:
                line[-1] = line[-1].replace("-","-TooShort-")
                ToxLociPartial.setdefault(k, [])
                ToxLociPartial[k].append("\t".join(line))
        #wrong codon frame but CDS longer than threshold
        if len(F[k]) % 3 != 0 and len(F[k]) >= CDSlength:
            for line in G[k]:
                line[-1] = line[-1].replace("-","-Partial-")
                ToxLociPartial.setdefault(k, [])
                ToxLociPartial[k].append("\t".join(line))
        #wrong codon frame and CDS longer than threshold
        if len(F[k]) % 3 != 0 and len(F[k]) < CDSlength:
            for line in G[k]:
                line[-1] = line[-1].replace("-","-PartialTooShort-")
                ToxLociPartial.setdefault(k, [])
                ToxLociPartial[k].append("\t".join(line))
    GTFout = open(Fout+"toxin_annotation.gtf", "w")
    for k in ToxLociFinal.keys():
        GTFout.write("\n".join(ToxLociFinal[k])+"\n")
    GTFout.close()
    if len(ToxLociPartial.keys()) != 0:
        GTFout = open(Fout+"annotation_removed.txt", "w")
        for k in ToxLociPartial.keys():
            GTFout.write("#"+k+" -> It may represent an erroneous annotation removed by ToxCodAn-Genome. Incorrect codon frame and/or CDS shorter than the specified -length- threshold\n")
            GTFout.write("\n".join(ToxLociPartial[k])+"\n")
            GTFout.write("\n")
        GTFout.close()
    return ToxLociN

#generate final gtf annotation
def _OutputFinalGTF_(Fout, ToxLoci):
    GTF = {}
    TOX = {}
    for i in sorted(os.listdir(Fout+"matched_GTFs/")):
        a = open(Fout+"matched_GTFs/"+i,"r")
        loci = i.replace(".gtf","")
        toxfam = ToxLoci[loci]
        TOX.setdefault(toxfam, 0)
        TOX[toxfam] += 1
        geneid = toxfam+"-"+str(TOX[toxfam])
        GTF.setdefault(i, [])
        for line in a:
            if "\tCDS\t" in line:
                line1 = line.strip().split("\t")
                line1[1] = "ToxCodAn"
                line1[-1] = "gene_id \""+geneid+"\"; transcript_id \""+geneid+"\";"
                GTF[i].append("\t".join(line1))
        a.close()
    GTFout = open(Fout+"pre_toxin_annotation.gtf", "w")
    for k in GTF.keys():
        GTFout.write("\n".join(GTF[k])+"\n")
    GTFout.close()

#output GTF from each loci
def _OutputMatchedGTF_(region, gtf, Fout):
    if os.path.isdir(Fout) == False:
        os.mkdir(Fout)
    contig = region.split("--")[0]
    startR = int(region.split("--")[1].split("-")[0])-100
    df = pd.read_csv(gtf,sep='\t',names=['seqid','mode','feature','st','end','sc1','strand','frame','desc'])
    df["st"] = df["st"] + startR
    df["end"] = df["end"] + startR
    df['seqid'] = contig
    df['mode'] = "TG:best_hit"
    df.to_csv(Fout+region+".gtf", sep="\t", index=False, header=False, quoting=csv.QUOTE_NONE)

#parse exonerate output and convert to gtf
def _ParseExonerate_(folder):
    gff = {}
    score = {}
    for i in os.listdir(folder):
        if i.endswith("_exonerate.out"):
            pos = i.replace("_exonerate.out","")
            if os.path.getsize(folder+i) > 600:
                a = open(folder+i,"r")
                for line in a:
                    if "\tgene\t" in line:
                        line1 = line.strip().split("\t")
                        score[pos] = int(line1[5])
                        strand = line1[6]
                        line2 =  line1[8].split(" ; ")
                        for feature in line2:
                            if "sequence" in feature:
                                protein = feature.replace("sequence ","")
                            if "identity" in feature:
                                identity = feature.replace("identity ","")
                            if "similarity" in feature:
                                similarity = feature.replace("similarity ","")
                        if strand == "+":
                            st = line1[3]
                            en = line1[4]
                            startline = line1[0]+"\t"+line1[1]+"\tstart_codon\t"+line1[3]+"\t"+str(int(st)+2)+"\t.\t+\t.\tgene_id \""+protein+"\"; transcript_id \""+protein+".t\";"
                            stopline = line1[0]+"\t"+line1[1]+"\tstop_codon\t"+str(int(en)-2)+"\t"+line1[4]+"\t.\t+\t.\tgene_id \""+protein+"\"; transcript_id \""+protein+".t\";"
                            gff.setdefault(pos,[])
                            gff[pos].append(startline)
                            gff[pos].append(stopline)
                        if strand == "-":
                            st = line1[3]
                            en = line1[4]
                            stopline = line1[0]+"\t"+line1[1]+"\tstop_codon\t"+line1[3]+"\t"+str(int(st)+2)+"\t.\t-\t.\tgene_id \""+protein+"\"; transcript_id \""+protein+".t\";"
                            startline = line1[0]+"\t"+line1[1]+"\tstart_codon\t"+str(int(en)-2)+"\t"+line1[4]+"\t.\t-\t.\tgene_id \""+protein+"\"; transcript_id \""+protein+".t\";"
                            gff.setdefault(pos,[])
                            gff[pos].append(startline)
                            gff[pos].append(stopline)
                    if "\texon\t" in line:
                        line1 = "\t".join(line.strip().split("\t")[:8]).replace("\texon\t","\tCDS\t")
                        gff.setdefault(pos,[])
                        gff[pos].append(line1+"\tgene_id \""+protein+"\"; transcript_id \""+protein+".t\";")
                a.close()
        if i.endswith("_dna.fasta"):
            seqID = i.replace("_dna.fasta","")
            region = SeqIO.to_dict(SeqIO.parse(folder+i,"fasta"))[seqID].seq

    #remove hits with no canonical start and stop codons
    toremove = set([])
    if len(gff.keys()) >= 1:
        for k in gff.keys():
            for l in gff[k]:
                if "start_codon" in l:
                    l1 = l.split("\t")
                    if l1[6] == "+":
                        if str(region)[int(l1[3])-1:int(l1[4])] != "ATG":
                            toremove.add(k)
                    if l1[6] == "-":
                        if str(region)[int(l1[3])-1:int(l1[4])] != "CAT":
                            toremove.add(k)
                if "stop_codon" in l:
                    l1 = l.split("\t")
                    if l1[6] == "+":
                        if str(region)[int(l1[3])-1:int(l1[4])] not in ["TAA", "TAG", "TGA"]:
                            toremove.add(k)
                    if l1[6] == "-":
                        if str(region)[int(l1[3])-1:int(l1[4])] not in ["TTA", "CTA", "TCA"]:
                            toremove.add(k)

    for k in list(toremove):
        gff.pop(k)
        score.pop(k)

    #check the best annotation among all (alignment score comparison)
    if len(score.keys()) >= 1:
        highscore = (0, "")
        for k in score.keys():
            if score[k] > highscore[0]:
                highscore = (score[k], k)
        selected = highscore[1]
    if len(score.keys()) == 0:
        gff["missing"] = "missing"
        selected = "missing"
    return selected, gff[selected]

#ToxCodAn-genome
def _ToxCodAn_(genome, cds, outF, cpu, percid, mingenesize, maxgenesize, length, keeptemp):
    #blast search
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> searching for toxin regions in genome using the database file...")
    sp.check_output("makeblastdb -dbtype nucl -in "+genome+" -out "+outF+"blastDB/dna", shell=True)
    sp.call("blastn -num_threads "+str(cpu)+" -evalue 0.01 -query "+cds+" -db "+outF+"blastDB/dna -out "+outF+"blast.out -outfmt \'6 qseqid sseqid pident evalue qstart qend qlen sstart send qseq sseq\'", shell=True)

    dna = SeqIO.to_dict(SeqIO.parse(genome, 'fasta'))
    cds = SeqIO.to_dict(SeqIO.parse(cds, 'fasta'))

    #filtering match regions
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> retrieving and filtering matched regions...")
    if os.path.getsize(outF+"blast.out") == 0:
        print('''

ToxCodAn-Genome stopped to run, due to the following issue:
    No Blast hits was found between the genome and the toxin database used as input.
    ''')
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> ToxCodAn-Genome suddenly finished, please check the issue above!")
        quit()
    blast = pd.read_csv(outF+"blast.out",sep='\t',names=['qseqid','sseqid','pident','evalue','qstart','qend','qlen','sstart','send','qseq','sseq'])
    blast = blast[blast['pident'] >= int(percid)]
    blast.sort_values(by=['sseqid','qseqid','sstart'], ascending=True, inplace=True)
    blast['strand'] = 'strand'
    removerows = []
    for index, row in blast.iterrows():
        if row['sstart'] < row['send']:
            blast.at[index,'strand'] = "plus"
        if row['sstart'] > row['send']:
            blast.at[index,'strand'] = "minus"
        if row['qstart'] == 1 and row['sseq'][:3] != "ATG":
            removerows.append(index)
        if row['qend'] == row['qlen'] and row['sseq'][-3:] not in ["TAA", "TAG", "TGA"]:
            removerows.append(index)

    blast.drop(removerows, axis=0, inplace=True)
    plus = blast.loc[blast['strand'] == "plus"]
    plus.sort_values(by=['sseqid','qseqid','sstart'], ascending=True, inplace=True)
    minus = blast.loc[blast['strand'] == "minus"]
    minus.sort_values(by=['sseqid','qseqid','send'], ascending=False, inplace=True)

    #get regions matching with full-length CDSs
    matchFL = {}
    temp = []
    #plus strand
    for index, row in plus.iterrows():
        if row['qstart'] == 1:
            temp.append((index,row['qstart'],row['qend'],row['sstart'],row['send'],row['strand']))
        if row['qstart'] > 1:
            if row['qend'] >= row['qlen']:
                if temp != []:
                    matchstart = temp[-1][3]
                    matchend = row['send']
                    if abs(matchstart-matchend) >= int(mingenesize) and abs(matchstart-matchend) <= int(maxgenesize):
                        matchFL.setdefault(row['qseqid'],[])
                        matchFL[row['qseqid']].append((row['sseqid'],matchstart,matchend,row['strand']))
                        temp = []
        #entire CDS inside an unique exon
        if row['qstart'] == 1 and row['qend'] >= row['qlen']:
            matchstart = row['sstart']
            matchend = row['send']
            matchFL.setdefault(row['qseqid'],[])
            matchFL[row['qseqid']].append((row['sseqid'],matchstart,matchend,row['strand']))
            temp = []

    #minus strand
    for index, row in minus.iterrows():
        if row['qstart'] == 1:
            temp.append((index,row['qstart'],row['qend'],row['sstart'],row['send'],row['strand']))
        if row['qstart'] > 1:
            if row['qend'] >= row['qlen']:
                if temp != []:
                    matchstart = row['send']
                    matchend = temp[-1][3]
                    if abs(matchstart-matchend) >= int(mingenesize) and abs(matchstart-matchend) <= int(maxgenesize):
                        matchFL.setdefault(row['qseqid'],[])
                        matchFL[row['qseqid']].append((row['sseqid'],matchstart,matchend,row['strand']))
                        temp = []
        #entire CDS inside an unique exon
        if row['qstart'] == 1 and row['qend'] >= row['qlen']:
            matchstart = row['send']
            matchend = row['sstart']
            matchFL.setdefault(row['qseqid'],[])
            matchFL[row['qseqid']].append((row['sseqid'],matchstart,matchend,row['strand']))
            temp = []

    if len(matchFL.keys()) == 0:
        print('''

ToxCodAn-Genome stopped to run, due to the following issue:
    No matching region containing a full-length toxin CDS was detected between the genome and the toxin database used as input.
    ''')
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> ToxCodAn-Genome suddenly finished, please check the issue above!")
        quit()
        quit()

    #handle overlap matches
    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> processing overlapped regions...")
    final = {}
    overlap = {}
    added = set([])
    MATCHREPORT = open(outF+"matched_regions.gtf","w")
    for k in matchFL.keys():
        for i in matchFL[k]:
            seqid = i[0]
            if i[1] < i[2]:
                MATCHREPORT.write(seqid+"\tTG:FL_hit\tmatch\t"+str(i[1])+"\t"+str(i[2])+"\t.\t+\t.\tID="+k+"\n")
            if i[1] > i[2]:
                MATCHREPORT.write(seqid+"\tTG:FL_hit\tmatch\t"+str(i[2])+"\t"+str(i[1])+"\t.\t+\t.\tID="+k+"\n")
            pos = i[0]+"--"+str(i[1])+"-"+str(i[2])
            st = i[0]+"--"+str(i[1])
            en = i[0]+"--"+str(i[2])
            if i[1] > i[2]:
                pos = i[0]+"--"+str(i[2])+"-"+str(i[1])
                st = i[0]+"--"+str(i[2])
                en = i[0]+"--"+str(i[1])
            final.setdefault(pos,set([]))
            final[pos].add(k)
            if len(overlap.keys()) != 0 and pos not in added:
                temp = []
                for position in overlap.keys():
                    stcheck = "--"+str(i[1])+"-"
                    encheck = "-"+str(i[2])
                    if i[1] > i[2]:
                        stcheck = "--"+str(i[2])+"-"
                        encheck = "-"+str(i[1])
                    if ( i[0] in position ) and ( ( stcheck in position ) or ( encheck in position ) ):
                        overlap[position].add(pos)
                        added.add(pos)
                if pos not in added:
                    temp.append(pos)
                if temp != []:
                    overlap.setdefault(pos, set([]))
                    overlap[pos].add(pos)
                    added.add(pos)
            if len(overlap.keys()) == 0:
                overlap.setdefault(pos, set([]))
                overlap[pos].add(pos)
                added.add(pos)
    MATCHREPORT.close()

    if os.path.isdir(outF+"exonerate_out/") == False:
        os.mkdir(outF+"exonerate_out/")

    ToxContigs = {}
    ToxLociN = {}
    ToxLoci = {}
    SEL = {}
    for k in overlap.keys():
        if os.path.isdir(outF+"exonerate_out/"+k) == False:
            os.mkdir(outF+"exonerate_out/"+k)
        id = k.split("--")[0]
        start = k.split("--")[1].split("-")[0]
        end = k.split("--")[1].split("-")[1]
        region = dna[id].seq[int(start)-100:int(end)+100] #100 nts up/downstream
        OUTdna = open(outF+"exonerate_out/"+k+"/"+k+"_dna.fasta","w")
        OUTdna.write(">"+k+"\n"+str(region)+"\n")
        OUTdna.close()
        for i in overlap[k]:
            cdsid = list(final[i])[0]
            toxfam = cdsid.split("_")[-1]
            if "SVMP" in toxfam:
                toxfam = "SVMP"
            if "SVSP" in toxfam:
                toxfam = "SVSP"
            ToxContigs.setdefault(toxfam, set([]))
            ToxContigs[toxfam].add(id)
            count=1
            for cdsid in final[i]:
                OUTcds = open(outF+"exonerate_out/"+k+"/"+i+"_cds_"+str(count)+".fasta","w")
                OUTcds.write(">"+cdsid+"\n"+str(cds[cdsid].seq)+"\n")
                OUTcds.close()
                #run exonerate
                genomic = outF+"exonerate_out/"+k+"/"+k+"_dna.fasta"
                cdsIN = outF+"exonerate_out/"+k+"/"+i+"_cds_"+str(count)+".fasta"
                exoout = outF+"exonerate_out/"+k+"/"+i+"_cds_"+str(count)+"_exonerate.out"
                sp.call("exonerate --bestn 1 -c "+str(cpu)+" --revcomp --showtargetgff --percent "+str(percid)+" --model est2genome " + cdsIN + " " + genomic + " > " + exoout, shell=True)
                count += 1
        selected, gff = _ParseExonerate_(outF+"exonerate_out/"+k+"/")
        if selected != "missing":
            OUTgff = open(outF+"exonerate_out/"+k+"/"+k+".gtf","w")
            OUTgff.write("\n".join(gff)+"\n")
            OUTgff.close()
            ToxLociN.setdefault(toxfam, []) #count toxin loci number
            ToxLociN[toxfam].append(selected.split("_cds_")[0])
            ToxLoci[k] = toxfam
            SEL[selected.split("_cds_")[0]] = k
            _OutputMatchedGTF_(k, outF+"exonerate_out/"+k+"/"+k+".gtf", outF+"matched_GTFs/")
        #if selected == "missing" and transcripts != None:
            #transcript search
            #print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> searching for toxin regions in genome using the transcripts file...")
            #genomic = outF+"exonerate_out/"+k+"/"+k+"_dna.fasta"
            #exoout = outF+"exonerate_out/"+k+"/transcripts_exonerate.out"
            #minimapout = outF+"exonerate_out/"+k+"/transcripts_minimap.sam"
            #sp.call("exonerate --bestn 1 -c "+str(cpu)+" --revcomp --showtargetgff --percent "+str(percid)+" --model est2genome " + transcripts + " " + genomic + " > " + exoout, shell=True)
            #sp.call("minimap2 -t"+str(cpu)+" -ax splice:hq -uf " + genomic + " " + transcripts + " > " + minimapout, shell=True)
            #print("missing", k)
            #minimap2 -t20 -ax splice:hq -uf ../Bfonsecai_v1.p_ctg.polish.fasta.mod.MAKER.soft.masked SB0060_toxins_cds.fasta > ToxinsCodingAnnotation_Bfon.sam
            #parse the transcripts, if it has a size higher than X bytes
            #call exonerate with all transcripts file
            #call convert .sam into gtf
            #open the gtf to catch all annotated regions (to not be used in exonerate search)

    print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> generating final output...")

    _OutputFinalGTF_(outF, ToxLoci)

    sp.check_output("gffread -x "+ outF+"pre_toxin_annotation_cds.fasta" +" -g "+ genome +" "+ outF+"pre_toxin_annotation.gtf", shell=True)

    #remove partials and incomplete annotations and CDSs shorter than "length"
    ToxLociN = _CheckPartials_(outF+"pre_toxin_annotation.gtf", outF+"pre_toxin_annotation_cds.fasta", outF, length)

    print("\n\t>>> Number of toxin loci identified in the genome:")
    TotalTox = 0
    for toxfam in sorted(ToxLociN.keys()):
        print("\t\t"+toxfam+" -> "+str(ToxLociN[toxfam]))
        TotalTox += ToxLociN[toxfam]
    print("\n\t\tTOTAL -> "+str(TotalTox))

    sp.check_output("gffread -x "+ outF+"toxin_annotation_cds.fasta" +" -g "+ genome +" "+ outF+"toxin_annotation.gtf", shell=True)
    _TranslateCDS_(outF)

    if keeptemp == "False" or keeptemp == False:
        os.system("rm -r "+outF+"blastDB/")
        os.system("rm -r "+outF+"exonerate_out/")
        os.system("rm "+outF+"blast.out")
        os.system("rm "+outF+"pre_toxin*")

##Options
def __main__():
    parser = OptionParser()
    parser.add_option("-g", "--genome", dest="genome", help="Mandatory - genome sequence in FASTA format, /path/to/genome.fasta", metavar="fasta", default=None)
    parser.add_option("-d", "--database", dest="database", help="Mandatory - database with coding sequences (CDSs) of toxins in FASTA format, /path/to/cds.fasta", metavar="fasta", default=None)
    parser.add_option("-C", "--cds", dest="cds", help="Optional - toxin coding sequences (CDSs) of the individual/species previously annotated from de-novo/genome-guided assembly in FASTA format, /path/to/toxin_cds.fasta", metavar="fasta", default=None)
    parser.add_option("-t", "--transcripts", dest="transcripts", help="Optional - transcripts recovered from venom tissue RNA-seq data for the individual/species using de-novo/genome-guided assembly methods in FASTA format, /path/to/transcripts.fasta", metavar="fasta", default=None)
    parser.add_option("-r", "--reads", dest="reads", help="Optional - pre-processed reads (i.e., adapters trimmed and low-quality reads removed) obtained from the toxin tissue of the species in FASTQ(.GZ) format. If single-end (or merged reads), specify only one file (e.g., path/to/reads.fastq). If paired-end, specify both files in a comma-separated format (e.g., path/to/reads_1.fastq,path/to/reads_2.fastq). If you also set transcripts file in the \"-t\" parameter, this parameter/file will be ignored.", metavar="fastq(.gz)", default=None)
    parser.add_option("-o", "--output", dest="output", help="Optional - output folder, /path/to/output_folder; if not defined, the output folder will be set in the current directory with the following name [ToxCodAnGenome_output]", metavar="folder", default=None)
    parser.add_option("-p", "--percid", dest="percid", help="Optional - threshold value used as the minimum percent identity between match CDSs and genome [default=80]", metavar="int", default="80")
    parser.add_option("-s", "--mingenesize", dest="mingenesize", help="Optional - threshold value used as the minimum size of a gene [default=400]", metavar="int", default="400")
    parser.add_option("-S", "--maxgenesize", dest="maxgenesize", help="Optional - threshold value used as the maximum size of a gene [default=50000]", metavar="int", default="50000")
    parser.add_option("-l", "--length", dest="length", help="Optional - minimum size of a CDS; it will remove any annotated CDS shorter than the specified threshold [default=200]", metavar="int", default="200")
    parser.add_option("-k", "--keeptemp", dest="keeptemp", help="Optional - keep temporary files. Use True to keep all temporary files or False to remove them [default=False]", metavar="boolean value", default="False")
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used in each step [default=1]", metavar="int", default="1")

    (options, args) = parser.parse_args()

    if options.genome == None or options.database == None:
        print(
        """

▄▄▄█████▓ ▒█████  ▒██   ██▒ ▄████▄   ▒█████  ▓█████▄  ▄▄▄       ███▄    █
▓  ██▒ ▓▒▒██▒  ██▒▒▒ █ █ ▒░▒██▀ ▀█  ▒██▒  ██▒▒██▀ ██▌▒████▄     ██ ▀█   █
▒ ▓██░ ▒░▒██░  ██▒░░  █   ░▒▓█    ▄ ▒██░  ██▒░██   █▌▒██  ▀█▄  ▓██  ▀█ ██▒
░ ▓██▓ ░ ▒██   ██░ ░ █ █ ▒ ▒▓▓▄ ▄██▒▒██   ██░░▓█▄   ▌░██▄▄▄▄██ ▓██▒  ▐▌██▒
  ▒██▒ ░ ░ ████▓▒░▒██▒ ▒██▒▒ ▓███▀ ░░ ████▓▒░░▒████▓  ▓█   ▓██▒▒██░   ▓██░
  ▒ ░░   ░ ▒░▒░▒░ ▒▒ ░ ░▓ ░░ ░▒ ▒  ░░ ▒░▒░▒░  ▒▒▓  ▒  ▒▒   ▓▒█░░ ▒░   ▒ ▒
    ░      ░ ▒ ▒░ ░░   ░▒ ░  ░  ▒     ░ ▒ ▒░  ░ ▒  ▒   ▒   ▒▒ ░░ ░░   ░ ▒░
  ░      ░ ░ ░ ▒   ░    ░  ░        ░ ░ ░ ▒   ░ ░  ░   ░   ▒      ░   ░ ░
             ░ ░   ░    ░  ░ ░          ░ ░     ░          ░  ░         ░
                           ░                  ░

          ██████  ███████ ███    ██  ██████  ███    ███ ███████
         ██       ██      ████   ██ ██    ██ ████  ████ ██
         ██   ███ █████   ██ ██  ██ ██    ██ ██ ████ ██ █████
         ██    ██ ██      ██  ██ ██ ██    ██ ██  ██  ██ ██
          ██████  ███████ ██   ████  ██████  ██      ██ ███████


>>>> ToxCodAn-Genome v1.0 January 2023 <<<<
      ****Use -h for help!****

BASIC USAGE:
toxcodan-genome.py -g genome.fasta -d CDS_database.fasta
        """)
        quit()

    if options.database != None and options.database.split(".")[-1] not in ["fasta", "fa", "fnn"]:
        print('''ToxCodAn-Genome was not able to run due to the ERROR below:
        You must set a CDS toxin file in FASTA format in the database parameter ("-d").

        If it is a compressed FASTA file in ZIP or GZ format, please decompress it accordingly.
            e.g., 'unzip file.fasta.zip' or 'gzip -d file.fasta.gz'
        ''')
        quit()

    if options.output != None:
        if not options.output.endswith("/"):
            options.output += "/"
    if options.output == None:
        CWD = os.getcwd()
        options.output = CWD+"/ToxCodAnGenome_output/"

    if os.path.isdir(options.output) == False:
        os.mkdir(options.output)

    if str(options.mingenesize)[-1] in ["K","k","M","m","G","g","B","b"] or str(options.maxgenesize)[-1] in ["K","k","M","m","G","g","B","b"]:
        print('''ToxCodAn-Genome was not able to run due to the ERROR below:
        The gene size parameters [-s/--mingenesize and -S/--maxgenesize] must be set using the nucleotide size with no abbreviations (such as K, M, G). Please, adjust it accordingly.
        For instance, if the user wants to set the minimum gene size to 1Kb, the parameter must be specified as "--mingenesize 1000"; if the user wants to set the maximum gene size to 150Kb, the parameter must be specified like "--maxgenesize 150000".
        ''')
        quit()

    if int(options.mingenesize) >= int(options.maxgenesize):
        print('''ToxCodAn-Genome was not able to run due to the ERROR below:
        The minimum gene size [-s/--mingenesize] must be shorter than the maximum gene size parameters [-S/--maxgenesize] (i.e., "--mingenesize" < "--maxgenesize"). If you set similar values, please consider using distinct values to have a proper gene size, which allows ToxCodAn-Genome to perform a robust search for matched regions and refinement of exon/intron boundaries.
        ''')
        quit()

    if options.genome != None and options.database != None:

        #if CDS file is set in the command line
        if options.cds != None:
            DB = _ParseFasta_(options.database)
            TR = _ParseFasta_(options.cds)
            S = open(options.output+"toxin_database.fasta", "w")
            for k in DB.keys():
                S.write(">"+k+"\n"+DB[k]+"\n")
            for k in TR.keys():
                if "_" not in k:
                    S.write(">"+k+"_TOXIN\n"+TR[k]+"\n")
                if "_" in k:
                    S.write(">"+k+"\n"+TR[k]+"\n")
            S.close()
            options.database = options.output+"toxin_database.fasta"


        print("""

▄▄▄█████▓ ▒█████  ▒██   ██▒ ▄████▄   ▒█████  ▓█████▄  ▄▄▄       ███▄    █
▓  ██▒ ▓▒▒██▒  ██▒▒▒ █ █ ▒░▒██▀ ▀█  ▒██▒  ██▒▒██▀ ██▌▒████▄     ██ ▀█   █
▒ ▓██░ ▒░▒██░  ██▒░░  █   ░▒▓█    ▄ ▒██░  ██▒░██   █▌▒██  ▀█▄  ▓██  ▀█ ██▒
░ ▓██▓ ░ ▒██   ██░ ░ █ █ ▒ ▒▓▓▄ ▄██▒▒██   ██░░▓█▄   ▌░██▄▄▄▄██ ▓██▒  ▐▌██▒
  ▒██▒ ░ ░ ████▓▒░▒██▒ ▒██▒▒ ▓███▀ ░░ ████▓▒░░▒████▓  ▓█   ▓██▒▒██░   ▓██░
  ▒ ░░   ░ ▒░▒░▒░ ▒▒ ░ ░▓ ░░ ░▒ ▒  ░░ ▒░▒░▒░  ▒▒▓  ▒  ▒▒   ▓▒█░░ ▒░   ▒ ▒
    ░      ░ ▒ ▒░ ░░   ░▒ ░  ░  ▒     ░ ▒ ▒░  ░ ▒  ▒   ▒   ▒▒ ░░ ░░   ░ ▒░
  ░      ░ ░ ░ ▒   ░    ░  ░        ░ ░ ░ ▒   ░ ░  ░   ░   ▒      ░   ░ ░
             ░ ░   ░    ░  ░ ░          ░ ░     ░          ░  ░         ░
                           ░                  ░

          ██████  ███████ ███    ██  ██████  ███    ███ ███████
         ██       ██      ████   ██ ██    ██ ████  ████ ██
         ██   ███ █████   ██ ██  ██ ██    ██ ██ ████ ██ █████
         ██    ██ ██      ██  ██ ██ ██    ██ ██  ██  ██ ██
          ██████  ███████ ██   ████  ██████  ██      ██ ███████

        """)
        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> starting ToxCodAn-Genome (v1.0 January 2023)...")
        CWD = os.getcwd()
        print("\tGenome file ->", options.genome)
        print("\tDatabase file ->", options.database)
        print("\tOutput folder ->", options.output)
        print("\tMinimum percent identity ->", options.percid)
        print("\tMinimum gene size ->", options.mingenesize)
        print("\tMaximum gene size ->", options.maxgenesize)
        print("\tMinimum size of a CDS ->", options.length)
        print("\tNumber of threads ->", options.cpu)

        #processing transcripts to detect CDSs
        if options.transcripts != None:
            print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> processing transcripts to detect toxin CDSs...")
            print("\tTranscripts file ->", options.transcripts)
            #rename transcripts header to avoid errors due to special characters
            a = open(options.transcripts,"r")
            renamed = open(options.output+"transcripts.fasta","w")
            countT = 1
            for line in a:
                if line.startswith(">"):
                    renamed.write(">transcript_"+str(countT)+"\n")
                    countT += 1
                if not line.startswith(">"):
                    renamed.write(line)
            a.close()
            renamed.close()
            #run CDSscreening
            sp.check_output("CDSscreening.py -t "+options.output+"transcripts.fasta -d "+options.database+" -o "+options.output+"screening/ -p "+str(options.percid)+" -c "+str(options.cpu) , shell=True)
            #concatenate CDS secreening and database
            DB = _ParseFasta_(options.database)
            TR = _ParseFasta_(options.output+"screening/cds_screening.fasta")
            S = open(options.output+"toxin_database.fasta", "w")
            for k in DB.keys():
                S.write(">"+k+"\n"+DB[k]+"\n")
            for k in TR.keys():
                if "_" not in k:
                    S.write(">"+k+"_TOXIN\n"+TR[k]+"\n")
                if "_" in k:
                    S.write(">"+k+"\n"+TR[k]+"\n")
            S.close()
            options.database = options.output+"toxin_database.fasta"
            print("\tUpdated database ->", options.database)
            os.system("rm -r "+options.output+"transcripts.fasta")
            #if user also set reads, the reads will be ignored
            if options.reads != None:
                print("""
    ***Warning***
        -> The user set both parameters \"--transcripts\" and \"--reads\".
        -> Only the transcripts will be considered.
        -> Ignoring reads file(s).
        -> If the user wants to run the transcriptome assembly module, please specify only reads file(s) within the "--reads" parameter or run the "TRassembly.py" script following its intructions.
                """)

        #processing reads
        if options.transcripts == None and options.reads != None:
            print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> assembling transcripts using the reads file...")
            print("\tReads file(s) ->", options.reads)
            #run GGassembly.py
            sp.check_output("TRassembly.py -g "+options.genome+" -r "+options.reads+" -o "+options.output+"assembly/ -c "+str(options.cpu) , shell=True)
            #rename transcripts header to avoid errors due to special characters
            a = open(options.output+"/assembly/transcripts.fasta", "r")
            renamed = open(options.output+"transcripts.fasta","w")
            countT = 1
            for line in a:
                if line.startswith(">"):
                    renamed.write(">transcript_"+str(countT)+"\n")
                    countT += 1
                if not line.startswith(">"):
                    renamed.write(line)
            a.close()
            renamed.close()
            #run CDSscreening.py
            sp.check_output("CDSscreening.py -t "+options.output+"transcripts.fasta -d "+options.database+" -o "+options.output+"/screening/ -p "+str(options.percid)+" -c "+str(options.cpu) , shell=True)
            #concatenate CDS secreening and database
            DB = _ParseFasta_(options.database)
            TR = _ParseFasta_(options.output+"screening/cds_screening.fasta")
            S = open(options.output+"toxin_database.fasta", "w")
            for k in DB.keys():
                S.write(">"+k+"\n"+DB[k]+"\n")
            for k in TR.keys():
                if "_" not in k:
                    S.write(">"+k+"_TOXIN\n"+TR[k]+"\n")
                if "_" in k:
                    S.write(">"+k+"\n"+TR[k]+"\n")
            S.close()
            options.database = options.output+"toxin_database.fasta"
            print("\tUpdated database ->", options.database)
            os.system("rm -r "+options.output+"transcripts.fasta")

        #run toxin annotation
        _ToxCodAn_(options.genome,
                    options.database,
                    options.output,
                    str(options.cpu),
                    str(options.percid),
                    str(options.mingenesize),
                    str(options.maxgenesize),
                    int(options.length),
                    str(options.keeptemp))

        print("\n"+dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S")+" >>>> ToxCodAn-Genome finished!")
        print("\n\t>>> Check the final annotation output: "+options.output+"toxin_annotation.gtf\n")

        #dependencies: python biopython pandas blast exonerate gffread
         #to the assembly module: hisat2 samtools stringtie trinity spades
        #conda create -n ToxcodanGenome -c bioconda python biopython pandas blast exonerate gffread hisat2 samtools stringtie trinity spades

if __name__ == '__main__':
    __main__()

#END