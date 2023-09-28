## Tutorial

A quick tutorial to use ToxCodAn-Genome.

- Download the [Balt_test.fasta](https://github.com/pedronachtigall/ToxCodAn-Genome/blob/main/Tutorial/Balt_test.fasta) file. It is composed of 17 genomic loci of toxins identified in the *Bothrops alternatus* genome.

- Download the Viperidae toxin database [here](https://github.com/pedronachtigall/ToxCodAn-Genome/tree/main/Databases) or using the command below.
    - ```wget https://raw.githubusercontent.com/pedronachtigall/ToxCodAn-Genome/main/Databases/Viperidae_db.fasta```

## Running ToxCodAn-Genome

```
toxcodan-genome.py -g Balt_test.fasta -d Viperidae_db.fasta
```

The user can set several optional parameters to run ToxCodAn-Genome as desired:
 - ```-o folder``` -> output folder ```/path/to/output_folder```; [default="ToxCodAnGenome_output"]
 - ```-p int``` -> threshold value used as the minimum percenty identity between match CDSs and genome [default=80]
 - ```-s int``` -> threshold value used as the minimum size of a gene [default=400]
 - ```-S int``` -> threshold value used as the maximum size of a gene [default=50000]
 - ```-l int``` -> minimum size of an annotated CDS; it will remove any annotated CDS shorter than the specified threshold [default=200]
 - ```-k boolean``` -> keep temporary files. Use True to keep all temporary files or False to remove them [default=False]
 - ```-c int``` -> number of threads to be used in each step [default=1]

:warning:**Warning**:warning:

We strongly recommend to use the default options but paying attention to the ```-c``` option, which will decrease the running time of ToxCodAn-Genome proportionally as the number of threads being used.

## Expected results

If ToxCodAn-Genome is running properly, you should see a similar output for the number of toxin loci annotated in the tutorial file.

```
	>>> Number of toxin loci identified in the genome:
		CTL -> 1
		HYAL -> 1
		NGF -> 1
		NUC -> 1
		PLA2 -> 1
		SVMP -> 5
		SVSP -> 4
		VEGF -> 1
		Vespryn -> 1
		Waprin -> 1
```
