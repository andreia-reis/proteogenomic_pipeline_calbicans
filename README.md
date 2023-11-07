# Proteogenomic pipeline for Candida albicans
Source code accompanying the manuscript entitled "A proteogenomic pipeline for the analysis of protein biosynthesis errors in the human pathogen Candida albicans.", published in BioRXiv (https://doi.org/10.1101/2023.10.31.564356).

## Author information
Inês Correia<sup>1a</sup>, Carla Oliveira<sup>1</sup>, Andreia Reis<sup>1</sup>, Ana R Guimarães<sup>1</sup>, Susana Aveiro<sup>2</sup>, Pedro Domingues<sup>2</sup>, Ana R Bezerra<sup>1</sup>, Rui Vitorino<sup>1</sup>, Gabriela Moura<sup>1</sup>, Manuel A. S. Santos<sup>1,3a</sup>

<sup>1</sup> Institute of Biomedical Sciences (iBiMED) and Department of Medical Sciences (DCM)
<sup>2</sup> Department of Chemistry, University of Aveiro, 3810-193, Aveiro, Portugal.
<sup>3</sup> Multidisciplinary Institute of Ageing (MIA-Portugal), University of Coimbra, 3004-504 Coimbra, Portugal

<sup>a</sup> Corresponding authors: Manuel A. S. Santos (mansilvasantos@uc.pt) and Inês Correia (inescorreia@ua.pt)

## Files description

1. CaPeaksMainFiles.R
   
   Input file(s):
   
   Output files(s):

   Description
   
2. CaPeaksValidation.R
   
   Input file(s):

   Output files(s):

   Description

3. CaSubCodonsAnalysis.R
   
   Input file(s):

   Output files(s):

   Description

4. CaAllCodonsAnalysisI.R

   Input file(s):

   Output files(s):

   Description

5. CaAllCodonAnalysisII.R

   Input file(s):

   Output files(s):

   Description

   
## Full workflow (described in the methods section)

![image](https://github.com/andreia-reis/proteogenomic_pipeline_calbicans/assets/19263451/63bbd7eb-6772-4485-9cb4-71b036e00e3c)

Figure 1. Workflow for mistranslation analysis. Protein extracts from actively growing cultures were resolved in SDS-PAGE and divided into 8 fractions. Bands were manually excised from the gel, digested, and injected separately in the MS system. Raw MS/MS data was analysed by PEAKS Xpro software. The list of identified peptides was filtered to remove duplicates and low-quality hits, and codon assignment and mistranslation frequency was calculated using R scripts. Numbers indicate the number of peptides identified or validated in each step of the pipeline from T0 (wild-type) strain. 
