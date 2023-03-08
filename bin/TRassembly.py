#!/usr/bin/env python3
'''
TrnascriptomeAssembly (TRassembly) - Automated pipeline to perform transcriptome assembly using short-reads
Author: Pedro G. Nachtigall (pedronachtigall@gmail.com)
'''

##Modules
import os
import subprocess as sp
from optparse import OptionParser

##Functions
def _RunAssembly_(genome, reads, outF, method, cpu, M):
    ##GenomeGuided mode
    if method in ["gg","GG","both","BOTH"]:
        print(">>> Running Genome Guided assembly <<<")
        outFGG = outF+"GG/"
        if os.path.isdir(outFGG) == False:
            os.mkdir(outFGG)
        #index genome
        sp.check_output("hisat2-build -p "+str(cpu)+" "+genome+" "+outFGG+"GINDEX", shell=True)
        #map reads in the genome (no --dta)
        if "," in reads:
            #paired-end mode
            reads1 = reads.split(",")[0]
            reads2 = reads.split(",")[1]
            sp.check_output("hisat2 -p "+str(cpu)+" --rg-id SAMPLE --rg SM:SAMPLE --summary-file "+outFGG+"hisat2_summary.txt -x "+outFGG+"GINDEX -1 "+reads1+" -2 "+reads2+" -S "+outFGG+"mapped.sam", shell=True)
        if "," not in reads:
            #single-end mode
            sp.check_output("hisat2 -p "+str(cpu)+" --rg-id SAMPLE --rg SM:SAMPLE --summary-file "+outFGG+"hisat2_summary.txt -x "+outFGG+"GINDEX -U "+reads+" -S "+outFGG+"mapped.sam", shell=True)

        #convert SAM to BAM (samtools) and delete temporary files
        sp.check_output("samtools view -@ "+str(cpu)+" -b -S -o "+outFGG+"mapped.bam "+outFGG+"mapped.sam", shell=True)
        sp.check_output("rm "+outFGG+"mapped.sam",shell=True)
        #extract only mapped reads
        sp.check_output("samtools view -@ "+str(cpu)+" -b -F 4 "+outFGG+"mapped.bam > "+outFGG+"mapped_filtered.bam", shell=True)
        sp.check_output("rm "+outFGG+"mapped.bam",shell=True)
        #sort bam file
        sp.check_output("samtools sort -@ "+str(cpu)+" "+outFGG+"mapped_filtered.bam -o "+outFGG+"mapped_sorted.bam", shell=True)
        sp.check_output("rm "+outFGG+"mapped_filtered.bam",shell=True)
        #index bam file
        sp.check_output("samtools index "+outFGG+"mapped_sorted.bam", shell=True)

        #assembly transcripts - StringTie
        sp.check_output("stringtie -p "+str(cpu)+" -o "+outFGG+"transcripts.gtf "+outFGG+"mapped_sorted.bam", shell=True)
        #recover transcripts - StringTie
        sp.check_output("gffread -w "+outFGG+"st_transcripts.fasta -g "+genome+" "+outFGG+"transcripts.gtf", shell=True)

        #assembly transcripts - trinitGG
        sp.check_output("Trinity --genome_guided_bam "+outFGG+"mapped_sorted.bam --max_memory "+M+" --genome_guided_max_intron 1000 --CPU "+str(cpu)+" --output "+outFGG+"trinityGG/ --full_cleanup", shell=True)

        #remove genome index
        sp.check_output("rm "+outFGG+"GINDEX*",shell=True)

        #concatenate assemblies to a final file
        os.system("cat "+outFGG+"st_transcripts.fasta "+outFGG+"trinityGG/Trinity-GG.fasta >> "+outF+"transcripts.fasta")

    ##De novo mode
    if method in ["dn","DN","both","BOTH"]:
        print(">>> Running De novo assembly <<<")
        outFDN = outF+"DN/"
        if os.path.isdir(outFDN) == False:
            os.mkdir(outFDN)
        if "," in reads:
            #paired-end mode
            reads1 = reads.split(",")[0]
            reads2 = reads.split(",")[1]
            sp.check_output("Trinity --seqType fq --left "+reads1+" --right "+reads2+" --max_memory "+M+" --CPU "+str(cpu)+" --output "+outFDN+"trinitydn/ --full_cleanup", shell=True)
            sp.check_output("rnaspades.py --threads "+str(cpu)+" --phred-offset 33 -1 "+reads1+" -2 "+reads2+" -o "+outFDN+"spades/", shell=True)
        if "," not in reads:
            #single-end mode
            sp.check_output("Trinity --seqType fq --single "+reads+" --max_memory "+M+" --CPU "+str(cpu)+" --output "+outFDN+"trinitydn/ --full_cleanup", shell=True)
            sp.check_output("rnaspades.py --threads "+str(cpu)+" --phred-offset 33 -s "+reads+" -o "+outFDN+"spades/", shell=True)

        #concatenate assemblies to a final file
        os.system("cat "+outFDN+"trinitydn.Trinity.fasta "+outFDN+"spades/transcripts.fasta >> "+outF+"transcripts.fasta")

    print(">>> Final assembly: "+outF+"transcripts.fasta")

##Options
def __main__():
    parser = OptionParser()
    parser.add_option("-g", "--genome", dest="genome", help="Mandatory - genome sequence in FASTA format, /path/to/genome.fasta", metavar="fasta", default=None)
    parser.add_option("-r", "--reads", dest="reads", help="Mandatory - pre-processed short-reads (i.e., adapters trimmed and low-quality reads removed) obtained from the toxin tissue of the species in FASTQ(.GZ) format. If single-end (or merged reads), specify only one file (e.g., path/to/reads.fastq). If paired-end, specify both files in a comma-separated format (e.g., path/to/reads_1.fastq,path/to/reads_2.fastq). the files can be compressed in GZ format.", metavar="fastq(.gz)", default=None)
    parser.add_option("-o", "--output", dest="output", help="Optional - output folder, /path/to/output_folder; if not defined, the output folder will be set in the current directory with the following name [assembly]", metavar="folder", default=None)
    parser.add_option("-m", "--method", dest="method", help="Optional - method used to assemble the transcriptome [gg/dn/both]. Choice between genome-guided (\"-m gg\"; which considers only genome-guided assembly), de novo (\"-m dn\"; which considers only de novo assembly) or both (\"-m both\"; which considers both methods). [default = both]", metavar="str", default="both")
    parser.add_option("-c", "--cpu", dest="cpu", help="Optional - number of threads to be used in each step [default=1]", metavar="int", default="1")
    parser.add_option("-M", "--memory", dest="M", help="Optional - Max memory usage to be passed to Trinity assembler [default=8G], use the same format as stated by Trinity assembler (e.g., 40G, 4G, 400M, etc)", metavar="string", default="8G")

    (options, args) = parser.parse_args()

    if options.genome == None or options.reads == None:
        print(
        """
>>>> Automated pipeline to peform Transcriptome Assembly <<<<
      ****Use -h for help!****

PAIRED-END READS USAGE:
TRassembly.py -g genome.fasta -r reads_1.fastq(.gz),reads_2.fastq(.gz)

SINGLE-END READS USAGE:
TRassembly.py -g genome.fasta -r reads.fastq(.gz)
        """)
        quit()

    if options.method not in ["gg","GG","dn","DN","both","BOTH"]:
        print("""
Error: the \"-m|--method\" pamater was wrongly assigned.
    Choose one of the following:
        gg -> to run genome-guided assembly
        dn -> to run de novo assembly
        both -> to run both methods [default]
        """)
        quit()

    if options.output != None:
        if not options.output.endswith("/"):
            options.output += "/"
    if options.output == None:
        CWD = os.getcwd()
        options.output = CWD+"/assembly/"

    if os.path.isdir(options.output) == False:
        os.mkdir(options.output)

    if options.genome != None and options.reads != None:

        _RunAssembly_(options.genome,
                        options.reads,
                        options.output,
                        options.method,
                        str(options.cpu),
                        options.M)

if __name__ == '__main__':
    __main__()

#END
