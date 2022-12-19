![ToxcodanGenome_logo](/ToxcodanGenome_logo.png)

# ToxCodAn-Genome
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
 - [GffRead](https://github.com/gpertea/gffread)
 - [Hisat2](http://daehwankimlab.github.io/hisat2/) - Optional (used in Transcriptome assembly)
 - [Samtools](http://www.htslib.org/) - Optional (used in Transcriptome assembly)
 - [StringTie](https://ccb.jhu.edu/software/stringtie/) - Optional (used in Transcriptome assembly)
 - [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) - Optional (used in Transcriptome assembly)
 - [SPAdes](https://github.com/ablab/spades) - Optional (used in Transcriptome assembly)

Ensure that all requirements are working properly.


# Toxin Database

We have designed toxin databases for some venomous lineages that were tested and working well.

The toxin databases with a brief descriptions are available in the ["Databases"](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Databases) folder.

You can check the [guide](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide) to learn how to design a custom toxin database (when venom tissue transcriptomic data is available or not).

The toxin database must be set with the ```-d``` option (see below).

# Usage
```
```
Basic usage:
```
toxcodan-genome.py -g genome.fasta -d toxin_database.fasta
```

Check our [tutorial](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Tutorial) to learn how to use ToxCodAn-Genome.

Check our [guide](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Guide) to have full details about ToxCodAn-Genome and how to perform toxin annotation in genomes.

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

Nachtigall et al. *under review*

# Contact
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

# Frequently Asked Questions (FAQ)

**[Q1]** What Operation System (OS) do I need to use ToxCodAn-Genome?
  - We tested ToxCodAn-Genome in Linux Ubuntu 16, 18 and 20. However, we believe that ToxCodAn-Genome should work on any UNIX OS able to have all dependencies installed.

**[Q2]** ToxCodAn-Genome is returning an error in the "generating final output" step similar to ```subprocess.CalledProcessError``` and ```Segmentation fault (core dumped)```. What should I do?
 - This error can be caused by one or more lines containing a huge sequence. Some tools and packages, like Bio::DB::Fasta used by GffRead, can't process a fasta file with lines containing more than 65,536 characters. So, if you have any large sequence in one unique line, do the following:
    - download the script [BreakLines.py](https://github.com/pedronachtigall/CodAn/blob/master/scripts/BreakLines.py)
    - run BreakLines script: ```python BreakLines.py genome.fasta genome_breaklines.fasta```
    - use the "genome_breaklines.fasta" to run ToxCodAn-Genome
