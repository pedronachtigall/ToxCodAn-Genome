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
 - [GffRead](https://github.com/gpertea/gffread)
 - [Hisat2](http://daehwankimlab.github.io/hisat2/) - Optional (used in Transcriptome assembly)
 - [Samtools](http://www.htslib.org/) - Optional (used in Transcriptome assembly)
 - [StringTie](https://ccb.jhu.edu/software/stringtie/) - Optional (used in Transcriptome assembly)
 - [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) - Optional (used in Transcriptome assembly)
 - [SPAdes](https://github.com/ablab/spades) - Optional (used in Transcriptome assembly)

Ensure that all requirements are working properly.


# Toxin Database

# Usage

Check our [tutorial](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Tutorial) to learn how to use ToxCodAn-Genome.

# Inputs

# Outputs

# Citation

If you use or discuss **ToxCodAn-Genome**, its guide, or any script available at this repository, please cite:

Nachtigall et al. *under review*

# Contact
:bug::sos::speech_balloon:

To report bugs, to ask for help and to give any feedback, please contact **Pedro G. Nachtigall**: pedronachtigall@gmail.com

# Frequently Asked Questions (FAQ)
