
library(data.table)
library(stringr)
library(seqRFLP)
library(IRanges)
library(seqRFLP)
library(dplyr)
library(Biostrings)
library(pracma)
library(annotate)
library(plyr)
library(filesstrings)
options(scipen=999)

# create AllCodons folder and work within
mainDir<-getwd()
subDir1<-"AllCodons"
file<-""
ifelse(!dir.exists(file.path(mainDir, subDir1)), dir.create(file.path(mainDir, subDir1)), FALSE)
file.copy(c("Ca_A22-s07-m01-r149_AT_CODING_29032022.fasta", "Ca_A22-s07-m01-r149_AT_PROTEIN_29032022.fasta", paste0("protein-peptidesAllCodons",file,".txt")), "AllCodons")
setwd(paste0(getwd(), "/AllCodons"))
getwd()

listpep <- read.table("protein-peptidesAllCodons.txt", sep=",", header=TRUE)
#listpep <- read.table("protein-peptidesAllCodonsF1.txt", sep=",", header=TRUE)
#listpep <- read.table("protein-peptidesAllCodonsF2.txt", sep=",", header=TRUE)

listpep<-listpep[,-c(1,25,30:32,34:56)] ###(CGD,ID,pepseq4:seqmut,seqmut2:end)
listpep<-listpep[order(listpep$CGDid),]

for (i in 1:nrow(listpep)){
  if(is.na(listpep$seqOriginal2[i])==TRUE){
    listpep$seqOriginal2[i]=listpep$pepSeq3[i]
  }
}

countINS<-c()
countDEL<-c()
for(i in 1:length(listpep$Peptide)){
  countINS[[i]]<-grep("\\(ins\\)", listpep$Peptide[i])
  countDEL[[i]]<-grep("\\(del\\s\\w+\\)", listpep$Peptide[i])
}
countINS[sapply(countINS, function(x) length(x)==0)] <- NA
countINS<- countINS[!is.na(as.vector(countINS))]
print(paste("There are", length(countINS), "insertions in this table"))
countDEL[sapply(countDEL, function(x) length(x)==0)] <- NA
countDEL<- countDEL[!is.na(as.vector(countDEL))]
print(paste("There are", length(countDEL), "deletions in this table"))

######UPLOAD ORIGINAL PROTEIN SEQUENCE######
ALLproteinseq<-readAAStringSet(file="Ca_A22-s07-m01-r149_AT_PROTEIN_29032022.fasta")
seq_name = names(ALLproteinseq)
sequence = paste(ALLproteinseq)
proteins <- data.frame(seq_name, sequence)
names(proteins)[1]<-"cgd"
names(proteins)[2]<-"proteinseq"
proteins$cgd<-gsub("(.+\\_.+\\_B)\\s+Allele\\s+of\\s+.+\\_.+A", "\\1", proteins$cgd)
proteins$proteinseq<-gsub("\\*$", "", proteins$proteinseq)
proteins$proteinseq<-gsub("\\*+", "X", proteins$proteinseq)
proteins<-proteins[order(proteins$cgd),]
finalproteins<-merge(listpep,proteins, by.x="CGDid", by.y="cgd", all.x=TRUE)

# match pattern between peptide and protein seq

sink("resultsfromMATCH.txt")
for (i in 1:length(finalproteins$CGDid)) {
  print(as.character(finalproteins$CGDid[i]))
  print(paste("line:",i))
  print(matchPattern(as.character(finalproteins$seqOriginal2[i]), as.character(finalproteins$proteinseq[i])))
}
closeAllConnections()

###################################
header<-c()
position<-c()
resultsFromMatch<-readLines("resultsfromMATCH.txt")
len<-length(resultsFromMatch)
cnone<-grep("NONE", resultsFromMatch)
print(paste("There are", length(cnone), "IDs with no match. Before continue, check the resultsfromMATCH.txt"))
c<-c()
for (i in 1:length(cnone)){
  c[i]<-cnone[i]-3
}
cnone<-append(cnone, c)
resultsFromMatch[cnone]<-0
for (i in 1:len) {
  if ((grepl("\\[[2-9]+\\]", resultsFromMatch[i])) == TRUE) {
    print("You are getting more than one match from the same peptide in the same protein")
    print("Check it manually, remove it and run again from line 74")
  }
}
###################################

ALLdnaseq<-readDNAStringSet(file="Ca_A22-s07-m01-r149_AT_CODING_29032022.fasta")
seq_name = names(ALLdnaseq)
sequence = paste(ALLdnaseq)
dna <- data.frame(seq_name, sequence)
names(dna)[1]<-"cgd"
names(dna)[2]<-"dnaseq"
dna$cgd<-gsub("(.+\\_.+\\_B)\\s+Allele\\s+of\\s+.+\\_.+A", "\\1", dna$cgd)
dna<-dna[order(dna$cgd),]
finaldna<-merge(finalproteins,dna, by.x="CGDid", by.y="cgd", all.x=TRUE)
finaldna$dnaLen<-nchar(as.character(finaldna$dnaseq))
finaldna$proteinLen<-nchar(as.character(finaldna$proteinseq))

write.table(finaldna, file="tablefinalproteomics.txt", sep="\t", na="N/A", row.names = FALSE)
write.csv(finaldna, file="tablefinalproteomics.csv", na="N/A", row.names = FALSE)
closeAllConnections()

# check match protein-DNA by length
table<-read.delim("tablefinalproteomics.txt", header=TRUE, sep="\t")
protein<-nchar(as.character(table$proteinseq))
DNA<-table$dnaLen
check<-(DNA==3*protein+3)
names(check)=table$CGDid
table(check)
as.character(table$CGDid[!check])
##if needed, change the original CGD DNA file and run again from line 94##
