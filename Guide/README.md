
<div align="center">
<center>

![ToxcodanGenome_logo](figures/ToxcodanGenome_logo.png)

# The guide to annotate toxin genes in genomes
### Pedro G. Nachtigall

</center>
</div>

### Summary

- [Introduction](#introduction)
- [Toxin annotation](#toxin-annotation)
	- [ToxCodAnGenome](#toxcodangenome)
	- [Checking annotations](#checking-annotations)
- [NonToxin annotation](#nontoxin-annotation)
- [Extra: Plotting toxin loci](#plotting-toxin-loci)

# Introduction
The Guide to performing toxin gene annotation in genomes is part of [**ToxCodAn-Genome**](https://github.com/pedronachtigall/ToxCodAn-Genome) and designed to walk you through our toxin annotation pipeline.

Before start and walking through the guide, it will be good to have some basic knowledge about bioinformatics.

If you don't have much experience in using command lines, running programs, or bioinformatics in general, then you can follow the ["Basic Bioinformatics"](https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide#basic-bioinformatics) section in the ToxCodAn's guide to venom gland transcriptomes. This section has some general resources and information that may help you to better understand bioinformatics and the rest of this Guide. The Guide is not designed to teach you *everything*, so we highly recommend working through the training resources first!

# Toxin annotation

## ToxCodAnGenome

### Toxin database

We have designed toxin database with curated toxin CDSs for ..., which can be downloaded [here](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Databases).

If you are working with some venomous lineage that is still not on our toxin database list, you can design a custom toxin database as decribed below.

### Custom toxin database

If you and your research group have been extensively working within the venomous taxa of the species being analyzed and have a well-curated set of toxin CDSs available in [FASTA](https://en.wikipedia.org/wiki/FASTA_format) format, it can be used as the toxin database or integrate some of the available toxin databases.

If using a custom toxin database as the only toxin database, just set the path to the custom database file with the parameter ```-d```. In this case, ensure you have the toxin family annotated in the header of each toxin CDS after an ```_``` ("underscore") symbol (e.g., ```>Sequence1_TOXIN```).
 - ***Tip:*** If your database only contains one toxin family, it can be easily done with [Perl](https://www.perl.org/):
     - ```perl -p -e 's/^(>.*)$/$1_TOXIN/g' custom_toxin_db.fasta > custom_toxin_db_renamed.fasta```
     - Replace ```"_TOXIN"``` to the target toxin family (e.g., for snake venom metalloproteinase use ```"_SVMP"``` and so on).

If you want to integrate some of the ToxCodAn-Genome's toxin database, just indicate the path to the custom toxin database file with the parameter ```-C```. In this case, you will not need to specify the toxin family in the header (as previously described); however, if it is not specified, the ToxCodAn-Genome will consider the sequence as a "generic" toxin, by adding a string ```_TOXIN``` at the end of header of each sequence (e.g., ```>Sequence1_TOXIN```). We strongly recommend to perform the annotation and add the toxin family in each sequence header to facilitate downstream analysis.

#### Curated toxin CDS set

To design a set of curated toxin CDSs, you may survey available data in sequence databases, such as [TSA](https://www.ncbi.nlm.nih.gov/genbank/TSA) and [Genbank](https://www.ncbi.nlm.nih.gov/genbank/) from NCBI, [EMBL-EBI](https://www.ebi.ac.uk/), [ENSEMBL](https://www.ensembl.org/index.html), [CNGBdb](https://db.cngb.org/) (China National GeneBank DataBase), and repositories and databases designed to specific genomes/purposes. You can take advantage of the well-designed search engine of these databases by using "keywords" that may help you to retrieve a high-quality set of toxin CDSs.

Here, we briefly described how to search for complete CDSs of the toxin SVMP (snake venom metalloproteinase) of vipers available in NCBI as an example, but you can easily modify the search parameters and use other databases to design a specific set of sequences for your purposes.

 - Access the GenBank website: https://www.ncbi.nlm.nih.gov/genbank/
 - Use the following set of keywords in the: ```Viperidae AND snake venom metalloproteinase AND "complete cds"```
 - Download the result in GenBank format:
	- "send to:" -> "Complete Record" -> "File" -> "Format: GenBank" -> "Creat File"
	- You can download the data in any other format, like downloading the "coding sequence" in FASTA format. But you may adjust the header of sequences accordingly (see ["Custom toxin database"](https://github.com/pedronachtigall/ToxCodAn-Genome/edit/main/Guide/README.md#custom-toxin-database) above).
 - Run our script to retrieve the toxin CDSs and have it ready to use with ToxCodAn-Genome:
	- ```ParseGenBank.py inut.gb output_cds.fasta TOXIN```
	- Replace ```"TOXIN"``` to the target toxin family (e.g., ```"SVMP"``` for snake venom metalloproteinase).

You can also search for RNA-seq data from venom tissue to be re-analyzed. These datasets can be downloaded, the toxin transcripts assembled and its toxin CDSs annotated and retrieved to be used as a Toxin Database. You can follow any available pipeline, [ToxCodAn's guide](https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide#the-guide) to venom gland transcriptomes, or the ["Venom tissue transcriptome"](https://github.com/pedronachtigall/ToxCodAn-Genome/edit/main/Guide/README.md#venom-tissue-transcriptome) section below.

Additionally, you can also survey for accession numbers and links to the datasets of toxin CDSs in published manuscripts (in the "data availability" sections) or within its supplementary files (some authors usually keep its curated sequences in supplementary files or deposited in other databases rather than NCBI's GenBank/TSA, such as [figshare](https://figshare.com/) and [dryad](https://datadryad.org/stash)). You can also contact authors to request such toxin CDSs when not available in any specific database and/or supplementary file.

 - ***Tip:*** To better curate a toxin set and ensure that you have designed a high-quality set, it is good to acquire knowledge about the toxins, their protein domains, and also about the venomous lineage being studied.
	- You can find details about venomous lineages and their toxins in [VenomZone](https://venomzone.expasy.org/), which is a good resource to retrieve information on venoms from six taxa (snakes, scorpions, spiders, cone snails, sea anemones, and insects), as well as on their biological targets and effects. It also has a description of the protein domains of each toxin family and links to retrieve its peptide sequences from [Uniprot](https://legacy.uniprot.org/).
	- You can also use the protein id from specific toxins to access its information in Uniprot and better understand the structure of domains of these toxins (e.g., https://legacy.uniprot.org/uniprot/Q9W6M5).
	- To identify protein domains and signal peptide in the identified toxins and check if it presents the expected toxin structure, you can use some web-servers like [HMMER](https://www.ebi.ac.uk/Tools/hmmer/), [InterProScan](https://www.ebi.ac.uk/interpro/search/sequence/), and [SignalP](https://services.healthtech.dtu.dk/service.php?SignalP).
	- You can also survey the literature to find reviews about specific toxin families and venomous lineages.

#### Venom tissue transcriptome

If you have venom tissue transcriptomic data available for the species being annotated, you should consider using this data in the toxin annotation step to improve the final set generated by ToxCodAn-Genome. In this sense, you can detect the toxin CDSs and annotate them by following your own pipeline, or follow the guide to venom gland transcriptomics available within [ToxCodAn](https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide#the-guide) tool and published in [Briefings in Bioinformatics](https://doi.org/10.1093/bib/bbab095).

If you and your research group do not have a strong background in performing such transcriptome annotation, we designed two scripts to help on this task that can be run separately or integrated into the main ToxCodAn-Genome pipeline.

<details>
<summary>Expand "Transcriptome Assembly" Section</summary>

<div align="center">
<center>

![TRassembly_workflow](figures/TRassembly_workflow.png)

</center>
</div>

**Pre-processing of reads**

First, ensure that the adapters are trimmed and low-quality reads filtered. It can be performed by using any available tool, such as [trim_galore!](https://github.com/FelixKrueger/TrimGalore), [fastp](https://github.com/OpenGene/fastp), and many others. Here, we just set a simple command to run ```trim_galore```:
 - ```trim_galore --paired --phred33 --length 75 -q 25 --stringency 1 -e 0.1 -o sample_trimmed sample_r1.fastq.gz sample_r2.fastq.gz```

After removing adapters and low-quality reads, you can move to the transcriptome assembly step.

**Transcriptome assembly**

We designed a script to run transcriptome assembly to run the genome-guided methods of [Trinity](https://github.com/trinityrnaseq/trinityrnaseq/wiki) and [StringTie](https://ccb.jhu.edu/software/stringtie/) and the *de novo* method of Trinity and [rnaSPAdes](https://cab.spbu.ru/software/spades/) to assemble most of the toxin transcripts in the dataset ([Holding et al., 2018](https://doi.org/10.3390/toxins10060249)).

If using paired-end reads:
```
TRassembly.py -g genome.fasta -r sample_r1.fastq.gz,sample_r2.fastq.gz -c 20 -M 20G
```

If using single-end reads (or merged reads):
```
TRassembly.py -g genome.fasta -r sample_reads.fastq.gz -c 20 -M 20G
```

 - The final transcriptome assembly can be found at ```assembly/transcripts.fasta```. The output directory can be changed by using the parameter ```-o```.
 - Please adjust the number of threads ```-c``` and memory usage ```-M``` accordingly to your system.
 - Run ```TRassembly.py -h``` to print the help message.
 - ```TRassembly.py``` handles both type of files: ```.fastq``` and ```.fastq.gz```.
 - It may take a while to finish!

You may also run each assembler and method separately and also consider using other tools (e.g., extender, NGEN, velvet, etc) to increase the probability of retrieving most of the toxin transcripts. We strongly recommend you to test several assemblers and paremeters to ensure a high-quality toxin transcripts assembly, because it was not thoroughly tested to all venomous lineage and it is a step that still needs improvement. You can follow the instructions to run other assemblers [here](https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide#transcriptome-assembly).

- ***Tip about transcriptome assembly:*** if using a paired-end read dataset, you may consider merging reads to improve the *de novo* assembly by using [PEAR](https://www.h-its.org/software/pear-paired-end-read-merger/). It will increase the size of the reads by merging pairs with high-quality overlap.
	- ```pear -k -j 20 -f sample_r1.fastq.gz -r sample_r1.fastq.gz -o sample_merged_reads```
	- ```-j``` is the number of threads, adjust accordingly.

 - ***Tip about Trinity:*** if using a paired-end read dataset, ensure the number of reads in both files match. If the number of reads does not match, Trinity may return an error and stop the assembling step. To avoid this issue, you can merge reads as described above or you can take advantage of [seqtk](https://github.com/lh3/seqtk) to subsample your dataset to a specific number of reads in both files.
	- ```seqtk sample -s100 sample_r1.fastq 5000000 > sample_sub1.fastq```
	- ```seqtk sample -s100 sample_r2.fastq 5000000 > sample_sub2.fastq```
	- Always use the same random seed to keep a proper pair in both files (parameters ```-s```, which is set to ```100``` in our example).
	- Here we set the number of reads to 5 million (```5000000```), but it must be adjusted accordingly to your dataset.


</details>
<br>

<details>
<summary>Expand "Toxin CDS annotation" Section</summary>

<div align="center">
<center>

![CDSscreening_workflow](figures/CDSscreening_workflow.png)

</center>
</div>

**Toxin CDS annotation**

We designed a script to screen a transcriptome assembly (assembled by our script ```TRassembly.py``` or any other pipeline) using any specific set of toxin CDS sequences. To identify toxin CDSs in the transcripome assembly, you can set one of the available Toxin databases or use any set of sequences with curated toxin CDSs designed by you.

The script ```CDSscreening.py``` can be easily run as follows:

```
CDSscreening.py -t transcripts.fasta -d CDS_database.fasta -c 20
```
 - The toxin CDS output can be found at ```screening/cds_screening.fasta```.
	- The output directory can be changed by using the parameter ```-o```.
 - Set the path to your transcriptome assembly accordingly.
	- e.g., if using the output from ```TRassembly.py``` just set ```-t assembly/transcripts.fasta```.
 - Set the path to your toxin CDS database accordingly (e.g., ```$PATH/to/Viperidae_db.fasta```).
 - Adjust the number of threads ```-c``` accordingly to your system.
 - Run ```CDSscreening.py -h``` to print the help message.

If you feel that some toxins may not being properly annotated by ```CDSscreening.py``` pipeline or wants to ensure that all toxins are being correctly screened, you can consider running [ToxCodAn](https://github.com/pedronachtigall/ToxCodAn) and follow its [guide](https://github.com/pedronachtigall/ToxCodAn/tree/master/Guide) to venom gland transcriptomics to perform a manual curation of the toxins present in the transcriptome being analyzed. You may also consider running other annotation tools, like [Venomix](https://bitbucket.org/JasonMacrander/venomix/src/master/), [Trinotate](https://github.com/Trinotate/Trinotate), and [Dammit](https://github.com/dib-lab/dammit).

</details>
<br>

### Running ToxCodAn-Genome

## Checking annotations

### Checking reliable annotations

### Checking warning annotations

### Checking matched regions with no annotation

# NonToxin annotation

To perform NonToxin annotation (which is a general gene annotation), you can use any available tool. Here, we describe how to use [funannotate](https://funannotate.readthedocs.io/en/latest/#). It is a package that performs an automated gene prediction and a functional annotation, by integrating several methods and tools. This program is well documented, and you can modify the commands in this guide to reach your purposes.

In this guide, we follow 4 easy-to-run steps:
 1. Training (using RNA-seq data)
 2. Prediction (using training data and additional protein databases)
 3. Updating (using RNA-seq data to back over predictions and add UTR annotations)
 4. Functional Annotation (annotate the final set of predictions with gene names, GO terms, and functional information)

Install funannotate and all dependencies to ensure it is properly working (you can follow the instructions [here](https://funannotate.readthedocs.io/en/latest/install.html) and run a test to ensure it is properly working (```funannotate test -t all --cpus 10```).

Download the database specific to the species being analyzed.
 - e.g., for snakes we use the tetrapoda database: ```funannotate setup -d funannotate_db -b tetrapoda```

Ensure to use a trimmed and filtered paired-end reads (see "Pre-processing of reads" subsection in "Transcriptome Assembly").

## Repetitive region annotation (*Masking*)
First, ensure you have a soft-masked version of the genome being analyzed. Soft-masked genomes use lower-case letters to indicate repetitive regions.

If you don't have a soft masked version of the genome being analyzed, just expand this section and follow its instructions.

<details>
<summary>Expand "Masking the genome" Section</summary>

**Masking the genome**

Here, we describe how to use the [Extensive De-novo TE Annotator (EDTA)](https://github.com/oushujun/EDTA) to identify and annotate transposable elements, long-terminal repeats (LTRs), and other repetitive regions in the genome.

```
perl EDTA.pl --genome genome.fasta --threads 20 --sensitive 1 --anno 1
```
 - Adjust the number of threads accordingly to your system (```--threads```)
 - It is important to note that this step is slow and will take a while (up to a few days).
	- This is slow due to the use of parameter ```--sensitive 1``` that also runs [RepeatModeler/RepeatMasker](https://www.repeatmasker.org/
) to improve its final annotation.

**Soft-masking**

EDTA automatically only hard-masks the genome (which consists of replacing repetitive regions from nucleotide to "N"). Then, follow the steps below to soft-mask the genome:

```
#rename the hard-masked genome
mv genome.fasta.mod.MAKER.masked genome.fasta.mod.MAKER.hard.masked

#create a new directory to perform soft-masking
mkdir softmask && cd softmask
ln -s ../genome.fasta .
ln -s ../genome.fasta.mod.EDTA.anno/genome.fasta.mod.EDTA.RM.out .

#mimic how EDTA performs hard-masking to soft-mask the genome
perl ~/.conda/envs/EDTA/share/EDTA/util/make_masked.pl -genome genome.fasta -rmout genome.fasta.mod.EDTA.RM.out -maxdiv 30 -minscore 1000 -minlen 1000 -hardmask 0 -misschar N -threads 20
mv genome.fasta.new.masked ../genome.fasta.mod.MAKER.soft.masked
cd ..
```

 - Use the soft-masked version of the genome (```genome.fasta.mod.MAKER.soft.masked```) in further steps.

</details>
<br>


## Training

```
funannotate train -i genome.fasta.mod.MAKER.soft.masked -o annotate \
	--left sample_tissue1_r1.fastq.gz sample_tissue2_r1.fastq.gz sample_tissue3_r1.fastq.gz \
	--right sample_tissue1_r2.fastq.gz sample_tissue2_r2.fastq.gz sample_tissue3_r2.fastq.gz \
	--no_trimmomatic --max_intronlen 30000 \
	--cpus 20 --species "Genus species"
```
 - ```--cpus``` - set number of threads accordingly to your system.
 - ```--left``` and ```--right``` - set any available RNA-seq data, but it must follow similar positions in "left" and "right" parameters.
 - ```--species``` - set the species accordingly.
 - This step is slow and will take a while (up to a few days).

## Prediction

```
funannotate predict -i genome.fasta.mod.MAKER.soft.masked -o annotate \
	--transcript_evidence annotate/training/funannotate_train.stringtie.gtf.fasta \
	--protein_evidence $FUNANNOTATE_DB/uniprot_sprot.fasta \
	--busco_db tetrapoda --busco_seed_species chicken --organism other --max_intronlen 30000 \
	--cpus 20 --species "Genus species"
```
 - ```--cpus``` - set number of threads accordingly to your system.
 - ```--busco_db``` - set this parameter accordingly.
 - ```--busco_seed_species``` - set this parameter accordingly.
 - ```--species``` - set the species accordingly.
 - ```$FUNANNOTATE_DB/uniprot_sprot.fasta``` - if you set funannotate correctly, it must work. Otherwise, set the path to ```uniprot_sprot.fasta```.
 - This step is slow and will take a while (up to a few days).

## Updating

```
funannotate update -i annotate --cpus 20
```
 - ```--cpus``` - set number of threads accordingly to your system.
 - This step is slow and will take a while.

## Functional Annotation

```
funannotate iprscan -i annotate -c 20 -m local --iprscan_path $PATH/to/IPRSCAN/interproscan.sh
funannotate annotate -i annotate --busco_db tetrapoda --cpus 20
```
 - ```--cpus``` - set number of threads accordingly to your system.
 - ```--busco_db``` - set this parameter accordingly.
 - ```$PATH/to/IPRSCAN/``` - set it accordingly to the installed path of ```interproscan.sh``` in your system.
 - Both steps are slow and will take a while.

# Plotting toxin loci
