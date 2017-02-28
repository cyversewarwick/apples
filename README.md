[![Build Status](https://travis-ci.org/cyversewarwick/apples.svg?branch=master)](https://travis-ci.org/cyversewarwick/apples)
[![](https://images.microbadger.com/badges/image/cyversewarwick/apples_conservation.svg)](https://microbadger.com/images/cyversewarwick/apples_conservation)
[![](https://images.microbadger.com/badges/version/cyversewarwick/apples_conservation.svg)](https://microbadger.com/images/cyversewarwick/apples_conservation)

# APPLES (Analysis of Plant Promoter-Linked Elements)
![APPLES logo](https://github.com/cyversewarwick/apples/blob/master/files/APPLES.png)

APPLES is a set of tools to analyse promoter sequences on a genome-wide scale. In this CyVerse-compatible version, two **main modules** are provided: 

 -  APPLES\_rbh: Find Orthologs as Reciprocal Best Hits
 -  APPLES\_conservation: Find Non-Coding Conserved Regions

In addition, the following tools are also exposed to the user:

 - APPLES\_utr: Extract sequences based on FASTA and GFF3 files

The following diagram illustrates the structure of these modules:

![APPLES workflow](https://github.com/cyversewarwick/apples/blob/master/files/APPLES_workflow.png)

## Background
The original APPLES package is described at [this address](http://www2.warwick.ac.uk/fac/sci/dcs/people/sascha_ott/tools_and_software/apples)

## Publications

* Nathaniel J. Davies, Peter Krusche, Eran Tauber and Sascha Ott, **Analysis of 5’ gene regions reveals extraordinary conservation of novel non-coding sequences in a wide range of animals**, BMC Evolutionary Biology, 2015, [doi: 10.1186/s12862-015-0499-6](http://dx.doi.org/10.1186/s12862-015-0499-6)

* Laura Baxter, Aleksey Jironkin, Richard Hickman, Jay Moore, Christopher Barrington, Peter Krusche, Nigel P. Dyer, Vicky Buchanan-Wollaston, Alexander Tiskin, Jim Beynon, Katherine Denby, and Sascha Ott, **Conserved Noncoding Sequences Highlight Shared Components of Regulatory Networks in Dicotyledonous Plants**, Plant Cell, 2012, [doi:10.1105/tpc.112.103010](http://dx.doi.org/10.1105/tpc.112.103010)

## Modules

#### APPLES\_conservation_multiple

##### Inputs

:whale: Species String
This is a string of Species names separated by ",".
 
 - Note that there is no "," behind the last species;
 - The first species is the central species;

Example: Species\_1,Species\_2,Species_3

:whale: Sequence Database Folder

With Species_1 being the central species, you will have the following folder structure:
```
<input_folder>
	+-- Species_1
	|   +-- PlantA.fa
	|   +-- PlantA.bed
	|   +-- PlantA_utr5.bed
	|   +-- PlantA_utr3.bed
	+-- Species_2
	|   +-- PlantA.fa
	|   +-- PlantA.bed
	|   +-- PlantA_utr5.bed
	|   +-- PlantA_utr3.bed
	|   +-- rbhSearch_result.txt
	+-- Species_3
	|   +-- PlantA.fa
	|   +-- PlantA.bed
	|   +-- PlantA_utr5.bed
	|   +-- PlantA_utr3.bed
	|   +-- rbhSearch_result.txt
	.
	.
```
See `/cyverseZone/home/shared/cyverseuk/apples_testdata/apples_conservation_multiple/app_short` for an example.

![Screenshot of APPLES_conservation_multiple on CyVerse DE](https://github.com/cyversewarwick/apples/blob/master/files/screenshot_conservation_multiple.png)

#### APPLES\_rbh
The APPLES\_rbh module finds Orthologs as Reciprocal Best Hits

[Run APPLES_rbh on CyVerse](https://de.cyverse.org/de/?type=apps&app-id=d99ca952-dbe2-11e6-9e37-0242ac120003)

##### Version History
 - 1.0

##### Inputs

 - `Protein FASTA` of Species A
 - `Protein FASTA` of Species B

#### UTR Tool
The APPLES\_utr module extracts sequences based on FASTA and GFF3 files of **a** species

![Screenshot of APPLES_utr on CyVerse DE](https://github.com/cyversewarwick/apples/blob/master/files/screenshot_utr.png)

[Run APPLES_utr on CyVerse](https://de.cyverse.org/de/?type=apps&app-id=d99ca952-dbe2-11e6-9e37-0242ac120003)

##### Version History

 - 1.1-stable Added parallelisation option [fa9ebdd]
 - 1.0 Simple version adopted from Grannysmith

##### Inputs
For a Species X:

 - `Gene FASTA*` - This is the file from which you wish to extract your sequences from. Provided that you have the matching GFF3 annotation, this file may be genome, scaffold or others based.
 - `GFF3*` - This is the file which annotates the FASTA file.
 - `Gene ID Identifier Text**` - This is the text which prefixes the Gene ID in the 9th column of the GFF3 file. Check your GFF3 to see what goes here.
 - `Sequence Length` - The number of bases which you wish to extract upstream.
 - `Stop at Neighbouring Gene` - Check this if you wish the sequence extraction to stop at neighbouring gene.
 - `Include the 5-prime UTR region` - Check to start the upstream at TSS so that the sequence include the UTR region. Otherwise start at 5-prime.

```
* - Sequences of a species are queried from a pair of FASTA and GFF3 files. This requires that the Sequence IDs in both files to match. In the FASTA file, this is the ID following the `>` charactor in the description lines; in the GFF3 file, this is the value stored in the first column of the gene lines (i.e. lines that says "gene" in the 3rd column).
```

```
** - To understand the `Gene ID Identifier Text` works, here are a couple of examples:

Use "ID=" if your `gff3` file looks like this:
`Niben101Scf00059        maker   gene    513034  528469  .       +       .       ID=Niben101Scf00059g04019;Alias=maker-Niben101Scf00059-snap-gene-4.18`

Use "ID=gene:" if your `gff3` file looks like this:
`1       tair    gene    31170   33153   .       -       .       ID=gene:AT1G01050;Name=PPA1;biotype=protein_coding;description=Soluble inorganic pyrophosphatase 1 [Source:UniProtKB/Swiss-Prot%3BAcc:Q93V56];gene_id=AT1G01050;logic_name=tair`
```

#### Conservation Module
The APPLES_conservation module finds Non-Coding Conserved Regions

![Screenshot of APPLES_conservation on CyVerse DE](https://github.com/cyversewarwick/apples/blob/master/files/screenshot_conservation.png)

[Run APPLES_conservation on CyVerse](https://de.cyverse.org/de/?type=apps&app-id=d99ca952-dbe2-11e6-9e37-0242ac120003)

##### Inputs
There are three sections of inputs for the conservation module. The first two are identical to that of the utr module with each one being for one of the two species. In the third section:

 - `Orthologs` - A total of 4 columns (tab-separated) are required in this file. Column 1: Species A's protein ID; Column 2: Species B's protein ID; Column 3: Species B's gene ID; Column 4: Species A's gene ID. i.e. "SpeciesA_proteinID SpeciesB_proteinID SpeciesB_geneID SpeciesA_geneID". This is the format in which results from the APPLES_rbh module are produced.
 - `Orthologs Mode` - Results from the Pseudo-Orthologs option is used as a controlled result which is only useful when compared with the result produced by using the correspoinding (proper) orthologs. If you don't know what it means, please use the default mode.
 - `Window Size` - The Seaweed algorithm aligns substrings of the given sequences (the length of which are specified in each species's "Sequence Length" argument) at a time. The length of this substring is called the "Window Size". It is recommended to use one of these values: 30 / 60 (default) / 80 / 100


##### Parallelisation
Use this following command to split the orthologs file:
`split -d --number=l/$(nproc) rbhSearch_result_PlantA_PlantB.txt rbhSearch_result_PlantA_PlantB.txt`

## Accessibility

Similar to all of the CyVerse UK applications developed at Warwick. There are 3 options when it comes to using our applications:

1. Via the [CyVerse Discovery Environment](https://de.cyverse.org/de/). This is the recommended approach to a new user. This is the easiest option since a full user interface is provided to the user.
2. Using the Docker images that are available on our [Docker Hub repository :whale:](https://hub.docker.com/u/cyversewarwick/). Each application/tool has a corresponding image.
3. With the source codes that are hosted on our [Github repository :octocat:](https://github.com/cyversewarwick). This approach will give you more information of how the application actually works. We are always looking to improve our code, so feel free to send us a pull request.

The modules related to APPLES can be searched on the [CyVerse Discovery Environment](https://de.cyverse.org/)  using the "**apples**" keyword in the application search box as shown in this screenshot:

![Search for APPLES on CyVerse DE](https://github.com/cyversewarwick/apples/blob/master/files/screenshot_search.png)






