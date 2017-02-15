# APPLES (Analysis of Plant Promoter-Linked Elements) #

APPLES is a set of tools to analyse promoter sequences on a genome-wide scale.

In this CyVerse-compatible version, two **main modules** are provided: 

 -  APPLES\_rbh: Find Orthologs as Reciprocal Best Hits
 -  APPLES\_conservation: Find Non-Coding Conserved Regions

In addition, the following tools are also exposed to the user:

 - APPLES_utr: Extract sequences based on FASTA and GFF3 files

## Background
The original APPLES package is described at [this address](http://www2.warwick.ac.uk/fac/sci/dcs/people/sascha_ott/tools_and_software/apples)

## Publications

* Nathaniel J. Davies, Peter Krusche, Eran Tauber and Sascha Ott, **Analysis of 5’ gene regions reveals extraordinary conservation of novel non-coding sequences in a wide range of animals**, BMC Evolutionary Biology, 2015, [doi: 10.1186/s12862-015-0499-6](http://dx.doi.org/10.1186/s12862-015-0499-6)

* Laura Baxter, Aleksey Jironkin, Richard Hickman, Jay Moore, Christopher Barrington, Peter Krusche, Nigel P. Dyer, Vicky Buchanan-Wollaston, Alexander Tiskin, Jim Beynon, Katherine Denby, and Sascha Ott, **Conserved Noncoding Sequences Highlight Shared Components of Regulatory Networks in Dicotyledonous Plants**, Plant Cell, 2012, [doi:10.1105/tpc.112.103010](http://dx.doi.org/10.1105/tpc.112.103010)

## Modules

#### APPLES\_rbh
1.0

#### UTR Tool

##### Version History
1.1-stable Added parallelisation option [fa9ebdd]
1.0 Simple version adopted from Grannysmith

"ID=" if your `gff3` file looks like this:
`Niben101Scf00059        maker   gene    513034  528469  .       +       .       ID=Niben101Scf00059g04019;Alias=maker-Niben101Scf00059-snap-gene-4.18`

"ID=gene:" if your `gff3` file looks like this:
`1       tair    gene    31170   33153   .       -       .       ID=gene:AT1G01050;Name=PPA1;biotype=protein_coding;description=Soluble inorganic pyrophosphatase 1 [Source:UniProtKB/Swiss-Prot%3BAcc:Q93V56];gene_id=AT1G01050;logic_name=tair`

#### Conservation Module
Inputs:
Sequences of a species are queried from a pair of FASTA and GFF3 files. This requires that the Sequence IDs in both files to match. In the FASTA file, this is the ID following the `>` charactor in the description lines; in the GFF3 file, this is the value stored in the first column of the gene lines (i.e. lines that says "gene" in the 3rd column).

##### Parallelisation

`split -d --number=l/$(nproc) rbhSearch_result_PlantA_PlantB.txt rbhSearch_result_PlantA_PlantB.txt`

The Seaweed algorithm aligns substrings of the given sequences (the length of which are specified in each species's "Sequence Length" argument) at a time. The length of this substring is called the "Window Size". It is recommended to use one of these values: 30 / 60 (default) / 80 / 100

## Accessibility

Similar to all of the CyVerse UK applications developed at Warwick. There are 3 options when it comes to using our applications:

1. Via the [CyVerse Discovery Environment](https://de.cyverse.org/de/). This is the recommended approach to a new user. This is the easiest option since a full user interface is provided to the user.
2. Using the Docker images that are available on our [Docker Hub repository :whale:](https://hub.docker.com/u/cyversewarwick/). Each application/tool has a corresponding image.
3. With the source codes that are hosted on our [Github repository :octocat:](https://github.com/cyversewarwick). This approach will give you more information of how the application actually works. We are always looking to improve our code, so feel free to send us a pull request.

The modules related to APPLES can be searched on the [CyVerse Discovery Environment](https://de.cyverse.org/)  using the "**apples**" keyword in the application search box as shown in this screenshot:
![Search for APPLES on CyVerse DE](https://github.com/cyversewarwick/apples/blob/master/files/screenshot_search.png)






