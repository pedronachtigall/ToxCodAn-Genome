![ToxcodanGenome_logo](/ToxcodanGenome_logo.png)

# ToxCodAn-Genome
[![Published in Giga Science](https://img.shields.io/badge/published%20in-Giga%20Science-blue)](https://doi.org/10.1093/gigascience/giad116)

**ToxCodAn-Genome** is a computational tool designed to annotate toxin genes in genomes of venomous lineages.

The guide for toxin annotation in genomes is [here](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide).

# Getting Started

# Installation

Clone the ToxCodAn-Genome respository and add the bin folder into your PATH:
```
git clone https://github.com/pedronachtigall/ToxCodAn-Genome.git
export PATH=$PATH:$PWD/ToxCodAn-Genome/bin/
```
 - ```echo "export PATH=$PATH:$PWD/ToxCodAn-Genome/bin/" >> ~/.bashrc``` to add it permanently.
    - Replace ```~/.bashrc``` to ```~/.bash_profile``` if needed.
 - Apply "execution permission" to executables if needed: ```chmod +x $PWD/ToxCodAn-Genome/bin/*```

## Requirements
 - [Python](https://www.python.org/), [biopython](https://biopython.org/), and [pandas](https://pandas.pydata.org/)
 - [NCBI-BLAST](https://www.ncbi.nlm.nih.gov/books/NBK279671/) (v2.9.0 or above)
 - [Exonerate](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
 - [Miniprot](https://github.com/lh3/miniprot)
 - [GffRead](https://github.com/gpertea/gffread)
 - [Hisat2](http://daehwankimlab.github.io/hisat2/) - Optional (used in Transcriptome assembly)
 - [Samtools](http://www.htslib.org/) - Optional (used in Transcriptome assembly)
 - [StringTie](https://ccb.jhu.edu/software/stringtie/) - Optional (used in Transcriptome assembly)
 - [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) - Optional (used in Transcriptome assembly)
 - [SPAdes](https://github.com/ablab/spades) - Optional (used in Transcriptome assembly)

Ensure that all requirements are working properly.

:warning: If the user wants to install ToxCodAn-Genome and all dependencies using [Conda environment](https://docs.conda.io/projects/conda/en/latest/user-guide/concepts/environments.html), follow the steps below:

- **Tip:** First, ensure that you have added all conda channels properly in the following order:
    - ```conda config --add channels defaults```
    - ```conda config --add channels bioconda```
    - ```conda config --add channels conda-forge```
- Create the environment:
    - ```conda create -n ToxcodanGenome -c bioconda python biopython pandas blast exonerate miniprot gffread hisat2 samtools stringtie trinity spades```
- Git clone the ToxCodAn-Genome repository and add the bin to your PATH:
    - ```git clone https://github.com/pedronachtigall/ToxCodAn-Genome.git```
    - ```echo "export PATH=$PATH:$PWD/ToxCodAn-Genome/bin/" >> ~/.bashrc```
        - Replace ```~/.bashrc``` to ```~/.bash_profile``` if needed.
- It may be needed to apply "execution permission" to all bin executables in "ToxCodAn-Genome/bin/":
    - ```chmod +x ToxCodAn-Genome/bin/*```
- Then, run ToxCodAn-Genome as described in the ["Usage"](https://github.com/pedronachtigall/ToxCodAn-Genome/edit/main/README.md#usage) section.
- To activate the environment to run ToxCodAn-Genome just use the command: ```conda activate ToxcodanGenome```
- To deactivate the environment just use the command: ```conda deactivate```

# Toxin Database

We have designed toxin databases for some venomous lineages that were tested and working well.

The toxin databases with a brief descriptions are available in the ["Databases"](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Databases) folder.

You can check the [guide](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide) to learn how to design a custom toxin database (when venom tissue transcriptomic data is available or not).

The toxin database must be set with the ```-d``` option (see below).

# Usage
```
Usage: toxcodan-genome.py [options]

Options:
  -h, --help            show this help message and exit
  -g fasta, --genome=fasta
                        Mandatory - genome sequence in FASTA format,
                        /path/to/genome.fasta
  -d fasta, --database=fasta
                        Mandatory - database with coding sequences (CDSs) of
                        toxins in FASTA format, /path/to/cds.fasta
  -C fasta, --cds=fasta
                        Optional - toxin coding sequences (CDSs) of the
                        individual/species previously annotated from de-
                        novo/genome-guided assembly in FASTA format,
                        /path/to/toxin_cds.fasta
  -t fasta, --transcripts=fasta
                        Optional - transcripts recovered from venom tissue
                        RNA-seq data for the individual/species using de-
                        novo/genome-guided assembly methods in FASTA format,
                        /path/to/transcripts.fasta
  -r fastq(.gz), --reads=fastq(.gz)
                        Optional - pre-processed reads (i.e., adapters trimmed
                        and low-quality reads removed) obtained from the toxin
                        tissue of the species in FASTQ(.GZ) format. If single-
                        end (or merged reads), specify only one file (e.g.,
                        path/to/reads.fastq). If paired-end, specify both
                        files in a comma-separated format (e.g.,
                        path/to/reads_1.fastq,path/to/reads_2.fastq). If you
                        also set transcripts file in the "-t" parameter, this
                        parameter/file will be ignored.
  -u fasta, --uniprot=fasta
                        Optional - path to uniprot/toxprot database to perform
                        an extra step of similarity search to compare the
                        annotated toxins to the uniprot/toxprot database. It
                        will output a report file in TXT format containing the
                        annotated toxin and the best hit in the
                        uniprot/toxprot database, alongside their percent
                        identity, alignment length, toxin length, and
                        uniprot/toxprot entry length.
  -o folder, --output=folder
                        Optional - output folder, /path/to/output_folder; if
                        not defined, the output folder will be set in the
                        current directory with the following name
                        [ToxCodAnGenome_output]
  -p int, --percid=int  Optional - threshold value used as the minimum percent
                        identity between match CDSs and genome [default=80]
  -s int, --mingenesize=int
                        Optional - threshold value used as the minimum size of
                        a gene [default=400]
  -S int, --maxgenesize=int
                        Optional - threshold value used as the maximum size of
                        a gene [default=50000]
  -l int, --length=int  Optional - minimum size of a CDS; it will remove any
                        annotated CDS shorter than the specified threshold
                        [default=200]
  -k boolean value, --keeptemp=boolean value
                        Optional - keep temporary files. Use True to keep all
                        temporary files or False to remove them
                        [default=False]
  -c int, --cpu=int     Optional - number of threads to be used in each step
                        [default=1]
```
Basic usage:
```
toxcodan-genome.py -g genome.fasta -d toxin_database.fasta
```

Check our [tutorial](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Tutorial) to learn how to use ToxCodAn-Genome.

Check our [guide](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide) to have full details about running ToxCodAn-Genome and how to perform toxin annotation in genomes.

# Inputs

ToxCodAn-Genome has the following inputs as mandatory:
 - The genome in FASTA format through the ```-g``` option.
 - The toxin database in FASTA format through the ```-d``` option.

# Outputs

By default, ToxCodAn-Genome outputs the following files:
```
ToxCodAnGenome_output/
├── annotation_removed.txt
├── annotation_warning.txt
├── matched_GTFs/
│   ├── contig_1--1234-5678.gtf
│   ├── contig_2--11234-15678.gtf
│   ├── ...
│   └── contig_N--111234-115678.gtf
├── matched_regions.gtf
├── toxin_annotation_cds.fasta
├── toxin_annotation.gtf
└── toxin_annotation_pep.fasta
```

Description of the output files:
```
toxin_annotation -> final toxin annotation files (including a gtf and two fasta files of CDSs and peptides)
annotation_warning.txt -> list of annotations in the final toxin annotation file that need manual inspection (may represent truncated paralogs, pseudogenes, and/or erroneous annotations)
annotation_removed.txt -> annotations that were removed from the final toxin annotation file (may represent erroneous/incomplete annotations)
matched_regions.gtf -> regions of genome matching to full-length toxin CDSs in the database (that returned or not a toxin annotation)
matched_GTFs/ -> folder containing all matched regions that returned a toxin annotation (each file is named by the genome position as follows "contig--start-end")
```

If you want to keep all temporary files, run ToxCodAn-Genome with the parameter ```-k True```.

# Citation

If you use or discuss **ToxCodAn-Genome**, its guide, or any script available at this repository, please cite:

Nachtigall et al. (2024) ToxCodAn-Genome: an automated pipeline for toxin-gene annotation in genome assembly of venomous lineages. Giga Science. DOI: [https://doi.org/10.1093/gigascience/giad116](https://doi.org/10.1093/gigascience/giad116)

# License

[GNU GPLv3](https://www.gnu.org/licenses/gpl-3.0.html)

# Contact
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

# Frequently Asked Questions (FAQ)

**[Q1]** What Operation System (OS) do I need to use ToxCodAn-Genome?
  - We tested ToxCodAn-Genome in Linux Ubuntu 16, 18 and 20. However, we believe that ToxCodAn-Genome should work on any UNIX OS able to have all dependencies installed.

**[Q2]** How long will take to ToxCodAn-Genome finish the analysis?
 - We tested ToxCodAn-Genome using a personal computer (6-Core i7 with 16Gb memory) and 6 threads (```-c 6```). It took less than 2 minutes to finish the analysis using a genome with an average size of 1.6Gb. If the user has more threads available for use, the running time will decrease.

**[Q3]** ToxCodAn-Genome is returning an error in the "generating final output" step similar to ```subprocess.CalledProcessError``` and ```Segmentation fault (core dumped)```. What should I do?
 - This error can be caused by one or more lines containing a huge sequence. Some tools and packages, like Bio::DB::Fasta used by GffRead, can't process a fasta file with lines containing more than 65,536 characters. So, if you have any large sequence in one unique line, do the following:
    - download the script [BreakLines.py](https://github.com/pedronachtigall/CodAn/blob/master/scripts/BreakLines.py) to "break" the genomic sequences into 100 nts per line
    - ```wget https://raw.githubusercontent.com/pedronachtigall/CodAn/master/scripts/BreakLines.py```
    - run BreakLines script: ```python BreakLines.py genome.fasta genome_breaklines.fasta```
    - use the "genome_breaklines.fasta" to run ToxCodAn-Genome. It will solve this issue.
