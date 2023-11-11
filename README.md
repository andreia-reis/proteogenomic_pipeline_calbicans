# Proteogenomic pipeline for Candida albicans
Source code accompanying the manuscript entitled "A proteogenomic pipeline for the analysis of protein biosynthesis errors in the human pathogen Candida albicans.", published in BioRXiv (https://doi.org/10.1101/2023.10.31.564356).

## Author information
Inês Correia<sup>1a</sup>, Carla Oliveira<sup>1</sup>, Andreia Reis<sup>1</sup>, Ana R Guimarães<sup>1</sup>, Susana Aveiro<sup>2</sup>, Pedro Domingues<sup>2</sup>, Ana R Bezerra<sup>1</sup>, Rui Vitorino<sup>1</sup>, Gabriela Moura<sup>1</sup>, Manuel A. S. Santos<sup>1,3a</sup>

<sup>1</sup> Institute of Biomedical Sciences (iBiMED) and Department of Medical Sciences (DCM)
<sup>2</sup> Department of Chemistry, University of Aveiro, 3810-193, Aveiro, Portugal.
<sup>3</sup> Multidisciplinary Institute of Ageing (MIA-Portugal), University of Coimbra, 3004-504 Coimbra, Portugal

<sup>a</sup> Corresponding authors: Manuel A. S. Santos (mansilvasantos@uc.pt) and Inês Correia (inescorreia@ua.pt)

## Files description

<b>1. CaPeaksMainFiles.R</b>
   
<i>- Input file(s):</i>

protein-peptides.csv

peptides.cvs

Ca_A22-s07-m01-r149_AT_PROTEIN_29032022.fasta 
   
   <i>- Output file(s):</i>

protein-peptides_NoSub.txt

protein-peptides_NoSubForAccept.txt

protein-peptides_Sub.txt

   <i>- Description:</i>

   This script generates the main files for validation (in directed searches) and codon assignment. In the protein-peptides file from PEAKS, the duplicates are removed (due to allelic duplication) and the file is separated into peptides with substitutions (Sub) and without mutations (noSubForAccept). The noSubForAccept file is further filtered to remove peptides whose PTMs do not pass the Ascore filter (20) and peptides with insertions and/or deletions that do not pass the ion intensity filter (5%) by matching with peptidesInDels.cvs file. The Subfile from protein-peptides.cvs is matched to the peptides.cvs file to remove entries with mutations not passing the ion intensity filter (5%). The script also generates a mutated proteome database, which will be used in a second round of PEAKS DB search for directed searches: it takes the mutated peptides without duplicates, searches the corresponding protein, and changes de protein sequence according to the mutated peptide.
   
<b>2. CaPeaksValidation.R</b>
   
<i>- Input file(s):</i>

protein-peptides_Sub.txt

protein-peptides_NoSubForAccept.txt

protein-peptides_NoSub.txt

 <i>- Output files(s):</i>

protein-peptides_SubAccepted.txt

protein-peptidesAllCodons.txt

<i>-  Description:</i>

   This script filters mut peptides by accepting only amino acid misincorporations whose base peptides (peptides without the mutation) are present in the mainfile from PEAKS (original one with all indels and PTMs before their validation – protein-peptides_NoSubForAccept). A file is created to analyse the correspondence (useful for areas comparison, for instance). This script also prepares a file with all peptides (the ones with accepted substitutions plus the filtered NoSub peptides from PEAKS – protein-peptides_NoSub).

<b>3. CaSubCodonsAnalysis.R</b>
   
 <i>- Input file(s):</i>

protein-peptides_SubAccepted.txt

Ca_A22-s07-m01-r149_AT_CODING_29032022.fasta

  <i>- Output files(s):</i>
   
tablefinalproteomics.txt

TOPcodonsSubst.txt

<i>-  Description:</i>
   
   This script assigns the top most frequent codons to amino acids and calculates codon/sub frequencies.

<b>4. CaAllCodonsAnalysisI.R</b>

<i>- Input file(s):</i>

Ca_A22-s07-m01-r149_AT_CODING_29032022.fasta

Ca_A22-s07-m01-r149_AT_PROTEIN_29032022.fasta

protein-peptidesAllCodons.txt

<i>- Output files(s):</i>

tablefinalproteomics.txt

<i>-  Description:</i>

   This script assigns all codons to amino acids and calculates codon/sub frequencies.

<b>5. CaAllCodonAnalysisII.R</b>

<i>-  Input file(s):</i>

tablefinalproteomics.txt

<i>-  Output files(s):</i>

ALLcounts_final_X.txt (X is the aa of your choice!) 

<i>- Description:</i>

For a specific aa (of your choice!) calculates the codon/sub frequencies
   
## Full workflow (described in the methods section)

![image](https://github.com/andreia-reis/proteogenomic_pipeline_calbicans/assets/19263451/63bbd7eb-6772-4485-9cb4-71b036e00e3c)

Figure 1. Workflow for mistranslation analysis. Protein extracts from actively growing cultures were resolved in SDS-PAGE and divided into 8 fractions. Bands were manually excised from the gel, digested, and injected separately in the MS system. Raw MS/MS data was analysed by PEAKS Xpro software. The list of identified peptides was filtered to remove duplicates and low-quality hits, and codon assignment and mistranslation frequency was calculated using R scripts. Numbers indicate the number of peptides identified or validated in each step of the pipeline from T0 (wild-type) strain. 
