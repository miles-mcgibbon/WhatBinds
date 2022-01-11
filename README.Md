# WhatBinds
---
This script is designed to fetch all protein sequences from the NCBI Protein database that match a specified "organism" and "protein" search query. It then performs several analyses:

1. Performs a multiple sequence alignment using `clustalo`
2. Plots the level of conservation between the fetched sequences and saves this as a .svg image using `plotcon`
3. Analyses the results of the multiple sequence alignment using `infoalign`
4. Scans all fetched protein sequences for any motifs present in the PROSITE motif database
5. Queries [BindingDB](https://www.bindingdb.org/) protein targets for similarity to any fetched sequences, and gets known ligands of good hits
6. Queries [The Protein Data Bank](https://www.rcsb.org/) protein structures for similarity to any fetched sequences, and gets known ligands of good hits
7. Clusters fetched ligands into groups of similar ligands based on their [Tanimoto similarity scores](https://www.surechembl.org/knowledgebase/84207-tanimoto-coefficient-and-fingerprint-generation)

# Installation
---
Currently, the only method to install the python enviroment and dependencies is via **conda**. I'm happy to create **virtualenv** and **docker** based alternatives if there is sufficient interest, but for now **conda** it is!

### Dependencies
- Conda
- clustalo
- emboss tools
- entrez-direct tools
- Linux OS

To install and set up the dependencies and environment for the WhatBinds script, simply run to following commands:

```
git clone https://github.com/miles-mcgibbon/WhatBinds.git

cd WhatBinds

chmod 755 setup.sh

./setup.sh
```

This will install all dependencies and create a conda environment with all the necessary python packages installed.

When the program is run for the first time, it will download the BindingDB database for local storage which is a few hundred Mb. Everytime the program is run, it checks for an updated version of BindingDB, and downloads and updates the local copy if it finds one.

# Running the Program
---

### Essential Arguments

The program has two **essential arguments** that it needs to work, **-organism** and **-protein**

```bash
python3 whatbinds.py -organism "TAXONOMY IDENTIFIER" -protein "PROTEIN NAME"
```

These arguments **must each be wrapped in double quotes** when supplied to the program. For example, if searching for glucose-6-phosphatase proteins in birds (*Aves*), the search would be:

```bash
python3 whatbinds.py -organism "Aves" -protein "glucose-6-phosphatase"
```

By default, the program looks for an **exact match** for the protein name. If you wish to search for partial matches, an \* can be added to either end of the protein search query to glob for more results (make sure this is **inside the double quotation marks!**):

```bash
# protein names ending in glucose-6-phosphatase
python3 whatbinds.py -organism "Aves" -protein "*glucose-6-phosphatase"

# protein names starting with glucose-6-phosphatase
python3 whatbinds.py -organism "Aves" -protein "glucose-6-phosphatase*"

# protein names containing glucose-6-phosphatase
python3 whatbinds.py -organism "Aves" -protein "*glucose-6-phosphatase*"
```

If one of these arguments is missing, you will be presented with an error telling you which is missing:
```bash
python3 whatbinds.py -organism "Aves"
>>>
ValueError: '-protein' is not in list
```

When the program is run, a directory is created to store the files and analysis results from the search. This directory is named based on the supplied search terms. so for the *Aves* example above a directory would be created called `Aves_glucose-6-phosphatase` to store working files and analysis output. \* symbols are substituted for the word "glob" in filenames to prevent confusion with bash commands. If the directory already exists this means an identical search has already been done, and the you will get a prompt:

```bash
WARNING - data from an identical search exists already!
Continuing will delete the existing data - continue? (y/N)
```

**Answering yes will delete the existing search folder and everything inside it**. Answering no will exit the program leaving the existing data untouched - you will need to rename the folder or move it before running the analysis again.

### Optional Arguments

By default, the program excludes isoforms, subdomains and partial sequences from the sequence results, and runs on one thread. These can be changed via the optional arguments below.

The four optional arguments for the program are as follows:

1. -threads : The number of threads you wish to use to run the program, default 1 if not supplied
2. -partial : If supplied, the search will return partial sequences that are normally excluded
2. -isoforms : If supplied, the search will return isoform sequences that are normally excluded
2. -subdomains : If supplied, the search will return subdomain sequences that are normally excluded

So to find glucose-6-phosphatase proteins from *Aves* including partial sequences, subdomains and isoforms, the command would be:

```bash
python3 whatbinds.py -organism "Aves" -protein "glucose-6-phosphatase" \
                           -partial -subdomains -isoforms
```

To run this search using 8 threads to speed up the analysis, the command would be:

```bash
python3 whatbinds.py -organism "Aves" -protein "glucose-6-phosphatase" \
                           -partial -subdomains -isoforms -threads 8
```

### Excess Sequences

When excess sequences (>1000) are returned from the search, you will see a prompt asking for a solution. For example, if searching:

```bash
python3 whatbinds.py -organism "rodent" -protein "*amyloid*" \
                           -threads 4 -partial -isoforms -subdomains

```

You will see the following prompt:

```bash
Found 1275 relevant sequences...
**********************
WARNING: Large number of sequences detected! (1275)
It is recommended not to exceed 1000 sequences in a single analysis
Available filtering methods: 

***COMPLEX METHODS***
0 - Equal number of sequences from each species will be fetched (recommended)


***SIMPLE METHODS***
Only the first 1000 sequences will be fetched - please pick the sort method:
1 - Accession
2 - Date Modified
3 - Date Released
4 - Default Order
5 - Organism Name
6 - Taxonomy ID


Or type EXIT to end the analysis
Method choice:
```

Here, simply enter the number (0-6) of the solution you would like to use to reduce the number of sequences that are fetched, or enter "EXIT" to abort the analysis.

# Analysis Results/Program Output
---
Throughout the running of the program, the analysis results are displayed to the screen. The results are also all saved as csv files (or an svg file in the case of the sequence conservation plot) in a `summary_files` folder inside the search specific directory. 

Again using the example search above (sped up for efficiency with multiple threads):


```bash
python3 whatbinds.py -organism "Aves" -protein "glucose-6-phosphatase" -threads 8
```

The analysis results and progess will be printed to screen as below:

1. The multiple sequence alignment is performed

```bash
Found 69 relevant sequences...
Sequences fetched!
Performing MSA...
Using 8 threads
Read 69 sequences (type: Protein) from Aves_glucose-6-phosphatase/Aves_glucose-6-phosphatase.prot.fa
not more sequences (69) than cluster-size (100), turn off mBed
Calculating pairwise ktuple-distances...
Ktuple-distance calculation progress done. CPU time: 0.40u 0.01s 00:00:00.41 Elapsed: 00:00:00
Guide-tree computation done.
Progressive alignment progress done. CPU time: 2.80u 0.10s 00:00:02.90 Elapsed: 00:00:01
Alignment written to Aves_glucose-6-phosphatase/Aves_glucose-6-phosphatase.aln
Sequences aligned!
```

2. The conservation between sequences is plotted and analysed:

```bash
Checking conservation:
Mean SeqLen : 360.507
Mean AlignLen : 400.536
Mean Gaps : 6.942
Mean GapLen : 40.029
Mean Ident : 315.986
Mean Similar : 12.203
Mean Differ : 32.319
Mean % Change : 20.793
```

The conservation plot is written to:

`Aves_glucose-6-phosphatase/summary_files/Aves_glucose-6-phosphatase_conservation_plot.svg`

![](aves_cons_plot.png)

and a summary of the alignment is written to:

`Aves_glucose-6-phosphatase/summary_files/Aves_glucose-6-phosphatase_alignment_info.csv`

With the following information:

|Sequence|# SeqLen|AlignLen|Gaps  |GapLen |Ident  |Similar|Differ |% Change          |
|------|--------|--------|------|-------|-------|-------|-------|------------------|
|0     |358.0   |488.0   |20.0  |130.0  |263.0  |14.0   |81.0   |46.106556         |
|1     |358.0   |488.0   |20.0  |130.0  |263.0  |14.0   |81.0   |46.106556         |
|2     |226.0   |311.0   |14.0  |85.0   |84.0   |43.0   |99.0   |72.990356         |
|3     |358.0   |488.0   |20.0  |130.0  |268.0  |12.0   |78.0   |45.081965999999994|
|4     |358.0   |488.0   |20.0  |130.0  |267.0  |12.0   |79.0   |45.286884         |
...
|169   |358.0   |488.0   |20.0  |130.0  |269.0  |11.0   |78.0   |44.877047999999995|
|170   |299.0   |405.0   |16.0  |106.0  |179.0  |46.0   |74.0   |55.80246700000001 |
|mean  |334.064 |461.099 |19.789|127.035|195.322|32.181 |106.561|57.877            |


3. The PROSITE motifs for each sequence are counted:

```bash
Counting PROSITE db motifs...
100%|-----------------------------------------------| 69/69 [00:03<00:00, 19.61it/s]


Total motif counts:
Motif AMIDATION found 62 times


Relative motif frequencies: 
Motif: AMIDATION present 0.898550 times per sequence
```

A summary of the motif analysis is written to:

`Aves_glucose-6-phosphatase/summary_files/Aves_glucose-6-phosphatase_PROSITE_motif_summary.csv`

with the following information:

|Motif|Count  |Count Per Sequence|Found in Sequences|Sequences Queried|Taxonomy Search Key|Protein Search Key|
|-----|-------|------------------|------------------|-----------------|-------------------|------------------|
|AMIDATION|78     |0.45614035087719296|['pku41280.1s', ... , 'xp_009899158.1s']|171              |Aves               |*glucose-6-phosphatase*|
|TNF_1|1      |0.005847953216374269|['xp_010199970.1s']|171              |Aves               |*glucose-6-phosphatase*|
|THIOL_PROTEASE_HIS|1      |0.005847953216374269|['opj74608.1s']   |171              |Aves               |*glucose-6-phosphatase*|



4. The PDB and BindingDB are queried with the sequences using `blastp` to attempt to find known ligands:

```bash
Querying BindingDB for known ligands...
100%|-----------------------------------------------| 69/69 [00:07<00:00,  9.82it/s]
Querying the Protein Data Bank for known ligands...
100%|-----------------------------------------------| 69/69 [00:16<00:00,  4.31it/s]
100%|-----------------------------------------------| 69/69 [00:00<00:00, 276.06it/s]
No structures or ligands found!
```
The criteria for a match used by the program is >90% identity and <0.01 evalue.

For the *Aves* glucose-6-phosphatase example, no ligands have been found. To demonstrate how the analysis is presented, below is the ligand analysis output for the search:

```bash
python3 whatbinds.py -organism "mammalia" -protein "ABC transporter" -threads 8
```

Which gives the printed output:

```bash
Ligand Cluster 1: (7 members)
CC(=O)Oc1ccc2n(C)c(SSc3c(C(=O)Nc4ccccc4)c4cc(OC(C)=O)ccc4n3C)c(C(=O)Nc3ccccc3)c2c1
CCOC(=O)C1=C(NC(C)=C(C1c1ccccc1Cl)C(=O)Nc1ccccn1)c1ccc(cc1)-c1c(C)nc2sccn12
CCC(=O)c1ccc(cc1NC(=O)c2ccc3ccccc3n2)c4cn(nn4)c5ccc(CCN6CCc7cc(OC)c(OCCOCCOCCOC)cc7C6)cc5
CC(C)C[CH]1N2[CH](Cc3c1[nH]c4cc(OC5CCCC5)ccc34)C(=O)N[CH](CCC(=O)OC(C)(C)C)C2=O
CCc1c2CN3C(=O)C4=C(C=C3c2nc5ccc(O)cc15)[C](O)(CC)C(=O)OC4
COc1cc2CCN(CCc3ccc(NC(=O)c4cc(OC)c(OC)cc4NC(=O)c5cnc6ccccc6c5)cc3)Cc2cc1OC
CC[C]1(O)C(=O)OCC2=C1C=C3N(Cc4cc5c(CN(C)C)c(O)ccc5nc34)C2=O
*************
Ligand Cluster 2: (2 members)
Nc1ncnc2n(cnc12)[CH]3O[CH](CO[P](O)(=O)O[P](O)(=O)O[P](O)(O)=O)[CH](O)[CH]3O
Nc1ncnc2n(cnc12)[CH]3O[CH](CO[P](O)(=O)O[P](O)(O)=O)[CH](O)[CH]3O
*************
Ligand Cluster 3: (1 members)
OCCS
*************
Ligand Cluster 4: (1 members)
[Mg++]
*************
Ligand Cluster 5: (1 members)
[O-][P]([O-])([O-])=O
*************
Ligand Cluster 6: (1 members)
OC[CH]1O[CH](O)[CH](O)[CH](O)[CH]1O
*************
Ligand Cluster 7: (1 members)
CC(=O)N[CH]1[CH](O)O[CH](CO)[CH](O)[CH]1O
*************
Ligand Cluster 8: (1 members)
OCCNCCNc1ccc(NCCNCCO)c2C(=O)c3c(O)ccc(O)c3C(=O)c12
*************
Ligand Cluster 9: (3 members)
CCCCCCCCC=CCCCCCCCC(=O)OC[CH](CO[P](O)(=O)OCCN)OC(=O)CCCCCCCC=CCCCCCCCC
CCCCCCCCCCCCCCCCCC(=O)OC[CH](CO[P](O)(=O)OCCN)OC(=O)CCCCCCCCC
CCCCCCCCCCCC(=O)OC[CH](CO[P](O)(=O)OCC[N+](C)(C)C)OC(=O)CCCCCCCCCCC
*************
Ligand Cluster 10: (1 members)
CN1CCN(CC1)Cc2ccc(cc2)C(=O)Nc3ccc(C)c(Nc4nccc(n4)c5cccnc5)c3
*************
Ligand Cluster 11: (2 members)
CC(C)CCC[CH](C)[CH]1CC[CH]2[CH]3CC=C4C[CH](O)CC[C]4(C)[CH]3CC[C]12C
CC(C)CCC[CH](C)[CH]1CC[CH]2[CH]3CC=C4C[CH](CC[C]4(C)[CH]3CC[C]12C)OC(=O)CCC(O)=O
*************
Ligand Cluster 12: (1 members)
C[C]12CC[CH]3[CH](CCc4cc(O[S](O)(=O)=O)ccc34)[CH]1CCC2=O
*************
```

A summary of the known ligand analysis is saved to:

`mammalia_ABC_transporter/summary_files/mammalia_ABC_transporter_known_ligand_structural_summary.csv`

with the following information contained:

|Ligand                                                                            |Cluster Group|ZINC Database ID|PubChem ID|PDB Ligand ID|Ki (nM)|Homologous PDB Structures                                                                                                                                                                         |
|----------------------------------------------------------------------------------|-------------|----------------|----------|-------------|-------|--------------------------|
|CC(=O)Oc1ccc2n(C)c(SSc3c(C(=O)Nc4ccccc4)c4cc(OC(C)=O)ccc4n3C)c(C(=O)Nc3ccccc3)c2c1|1            |ZINC26161005    |5328602.0 |             |       |1IS0,1NZL,...,7NER|
|CCOC(=O)C1=C(NC(C)=C(C1c1ccccc1Cl)C(=O)Nc1ccccn1)c1ccc(cc1)-c1c(C)nc2sccn12       |1            |ZINC26838738    |15005345.0|             |       |                                                                                                                                                                                                  |


Ligands are clustered if their pairwise Tanimoto similarity is >~0.8.

### Program Filestructure

The filestructure created by a single search is as follows:

```bash
# search query
python3 whatbinds.py -organism "Aves" -protein "glucose-6-phosphatase"

# output
|--- Aves_glucose-6-phosphatase
|   |--- blast_results # contains one blastp output csv file for each sequence file
|   |--- split_seqs # contains individual sequence files
|   |--- summary_files # contains summary files for the user
|   |   |--- Aves_glucose-6-phosphatase_alignment_info.csv 
|   |   |--- Aves_glucose-6-phosphatase_conservation_plot.svg
|   |   |--- Aves_glucose-6-phosphatase_PROSITE_motif_summary.csv
|   |   |--- # would also be a ligand summary file here if any had been found
|   |--- working_files # contains files produced by emboss tools and blast databases
|   |   |--- Aves_glucose-6-phosphatase_alignment.info
|   |   |--- Aves_glucose-6-phosphatase.aln
|   |   |--- Aves_glucose-6-phosphatase.prot.fa
|   |   |--- bindingdb.pdb
|   |   |--- bindingdb.phr
|   |   |--- bindingdb.pin
|   |   |--- bindingdb.pot
|   |   |--- bindingdb.psq
|   |   |--- bindingdb.ptf
|   |   |--- bindingdb.pto

 ```
 
Files contained in `working_files` are initially written to the parent directory (e.g. `Aves_glucose-6-phosphatase`) and moved into the `working_files` directory at the end of the script.
