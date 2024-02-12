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

peptides.csv

peptides_InDels.cvs

Ca_A22-s07-m01-r149_AT_PROTEIN_29032022.fasta 
   
   <i>- Output file(s):</i>

protein-peptides_NoSub.txt

protein-peptides_NoSubForAccept.txt

protein-peptides_Sub.txt

   <i>- Description:</i>

   This script generates the main files for peptide validation and codon assignment. In the protein-peptides file from PEAKS, the duplicates are removed (due to allelic duplication) and the file is divided into two new files containing: peptides with substitutions (Sub) and without substitutions (noSubForAccept). The noSubForAccept file is further filtered to remove peptides whose PTMs do not pass the Ascore filter (20), and therefore are not depicted in the PTM column at protein-peptides file, and peptides with insertions and/or deletions that do not pass the ion intensity filter (5%) by matching with peptidesInDels.cvs file. The peptidesInDels.cvs is the peptides.cvs file extracted from PEAKS that contains only peptides with insertions and deletions. The Sub file is matched to the peptides.cvs file containing only mutations found by SPIDER to remove entries with mutations not passing the ion intensity filter (5%).
   
<b>2. CaPepValidation.R</b>
   
<i>- Input file(s):</i>

protein-peptides_Sub.txt

protein-peptides_NoSubForAccept.txt

protein-peptides_NoSub.txt

 <i>- Output files(s):</i>

protein-peptides_SubAccepted.txt

protein-peptidesAllCodons.txt

protein-peptides_SubAcceptedF1.txt

protein-peptidesAllCodonsF1.txt

protein-peptides_SubAcceptedF2.txt

protein-peptidesAllCodonsF2.txt

<i>-  Description:</i>

   This script filters sub peptides by accepting only amino acid misincorporations whose base peptides (peptides without the mutation) are present in the mainfile from PEAKS (original one with all indels and PTMs before their validation – protein - peptides_NoSubForAccept). An intermediate file is created to analyse the correspondence. This script also prepares a final file with all filtered peptides (the ones with validated substitutions and accepted NoSub peptides): protein-peptides_AllCodons. 
   Other files with different filtering parameters are also created (F1 and F2). 
   F1 removes duplicated peptides (with and without substitutions) whose CGDid, peptide sequence and mass are the same (considers duplicates those peptides that contain the same PTMs in different locations).
   F2 removes duplicated peptides (with and without substitutions) whose CGDid and peptide sequence are the same (considers duplicates those peptides with same base sequence, independently of having PTMs).

<b>3. CaSubCodonsAnalysis.R</b>
   
 <i>- Input file(s):</i>

protein-peptides_SubAccepted.txt

Ca_A22-s07-m01-r149_AT_CODING_29032022.fasta

  <i>- Output files(s):</i>
   
pepsubst.fasta

pepsubst_nosymbol.fasta

originalprotein.fasta

DupliPeptide.txt

DupliPepDifferentCodon.txt

tablefinalproteomics.txt/.cvs

TOPcodonsSubst.txt

<i>-  Description:</i>
   
   This script assigns codons to amino acids and calculates codon/sub frequencies.
Protein-peptide correspondence is analysed to identify putative errors on codon assignment due multiple peptide matching within the same protein. Duplicates are removed: same peptide matched to different proteins. As DNA sequence may be different, the information is kept so the putative error on codon assignment can be analysed. The Protein-DNA correspondence is analysed to analyse putative discrepancies between the files.

<b>4. CaAllCodonsAnalysisI.R</b>

<i>- Input file(s):</i>

Ca_A22-s07-m01-r149_AT_CODING_29032022.fasta

Ca_A22-s07-m01-r149_AT_PROTEIN_29032022.fasta

protein-peptidesAllCodons.txt

<i>- Output files(s):</i>

tablefinalproteomics.txt

<i>-  Description:</i>

   This script assigns codons to amino acids from all validated and accepted peptides. This is necessary to obtain the total number of codons encoding a specific amino acid and calculate the mistranslation frequency. Total codon frequency is obtained from the file protein-peptidesAllCodons. Protein-peptide correspondence is analysed to identify putative errors on codon assignment due multiple peptide matching within the same protein. The Protein-DNA correspondence is analysed to analyse putative discrepancies between the files.

<b>5. CaAllCodonAnalysisII.R</b>

<i>-  Input file(s):</i>

tablefinalproteomics.txt

<i>-  Output files(s):</i>

tablefinalproteomicsDupli_with_X.txt

DupliPeptide.txt

DupliPepDifferentCodon.txt

tablefinalproteomics_with_X.txt (X is the aa of your analysis!)

ALLcounts_final_X.txt

<i>- Description:</i>

This script analyses which codons, and how many, are assigning each specific amino acid (ex. Histidine is being codified by 8621 CAU and 4832 CAC codons). This is necessary to obtain the total number of codons encoding a specific amino acid so mistranslation frequency can be calculated. For each amino acid, the script checks if all matches are possible (absence of NA or NNN) and removes duplicates: same peptide matched to different proteins. When the assigned codon is different, the information is kept so the putative error on codon assignment can be analysed.
   
## Full workflow (described in the methods section)

![image](https://github.com/andreia-reis/proteogenomic_pipeline_calbicans/assets/19263451/63bbd7eb-6772-4485-9cb4-71b036e00e3c)

Figure 1. Workflow for mistranslation analysis. Protein extracts from actively growing cultures were resolved in SDS-PAGE and divided into 8 fractions. Bands were manually excised from the gel, digested, and injected separately in the MS system. Raw MS/MS data was analysed by PEAKS Xpro software. The list of identified peptides was filtered to remove duplicates and low-quality hits, and codon assignment and mistranslation frequency was calculated using R scripts. Numbers indicate the number of peptides identified or validated in each step of the pipeline from T0 (wild-type) strain. 

[![DOI](https://zenodo.org/badge/715792884.svg)](https://zenodo.org/doi/10.5281/zenodo.10651621)
