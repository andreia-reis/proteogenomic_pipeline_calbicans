
library(data.table)
library(stringr)
library(seqRFLP)
library(IRanges)
library(protr)
library(seqRFLP)
library(dplyr)
library(Biostrings)
library(pracma)
library(annotate)

listpep <- read.csv("protein-peptides_SubAccepted.txt")
#listpep <- read.csv("protein-peptides_SubAcceptedF1.txt")
#listpep <- read.csv("protein-peptides_SubAcceptedF2.txt")

listpep<-listpep[order(listpep$CGDid),]
listpep<-listpep[!sapply(listpep, function(x) all(is.na(x)))]
pepsubst <- data.frame(listpep$CGDid,listpep$pepSeq2)
pepsubst.fasta = dataframe2fas(pepsubst, file="pepsubst.fasta")
pepsubst_no <- data.frame(listpep$CGDid,listpep$seqOriginal2)
pepsubst_no.fasta = dataframe2fas(pepsubst_no, file="pepsubst_nosymbol.fasta")

finalproteins<-listpep[,c(1,48)] ##change here for proteinSeq column if needed
originalprotein.fasta = dataframe2fas(finalproteins, file="originalprotein.fasta")

# match pattern between peptide and protein seq

sink("resultsfromMATCH.txt")
for (i in 1:length(listpep$CGDid)) {
  print(as.character(listpep$CGD[i]))
  print(matchPattern(as.character(listpep$seqOriginal2[i]), as.character(listpep$proteinseq[i])))
}
closeAllConnections()

############################################
header<-c()
position<-c()
resultsFromMatch<-readLines("resultsfromMATCH.txt")
len<-length(resultsFromMatch)
cnone<-grep("NONE", resultsFromMatch)
print(paste("There are", length(cnone), "IDs with no match. Before continue, please check the resultsfromMATCH.txt"))
c<-c()
for (i in 1:length(cnone)){
  c[i]<-cnone[i]-3
}
cnone<-append(cnone, c)
resultsFromMatch[cnone]<-0
for (i in 1:len) {
  if ((grepl("\\[[2-9]+\\]", resultsFromMatch[i])) == TRUE) {
    print("You are getting more than one match from the same peptide in the same protein")
    print("Check it manually, remove it and run again from line 37")
  }
}
############################################

fileDNA<-readFASTA(file="Ca_A22-s07-m01-r149_AT_CODING_29032022.fasta")
n <- 3
ID=c()
codon<-c()
colsProt <- as.vector(grep("protPos.+", colnames(listpep)))
for(i in 1:length(listpep$CGDid)) {
  for (j in 1:length(fileDNA)) {
  if (names(fileDNA[j]) == listpep$CGDid[i]) {
    #print(paste0("i:", i))
    #print(paste0("j:", j))
    codonstemp<-sapply(seq(1,nchar(fileDNA[[j]]),by=n), function(x) substr(fileDNA[[j]], x, x+n-1))
    ID[i]<-names(fileDNA[j])
    codon[[i]]<-codonstemp[as.numeric(listpep[i,colsProt])]
    listpep$dnaLen[i]<-nchar(fileDNA[[j]])
  }
}
}
listpep$posDNA<-as.character(codon)
listpep$posDNA<-gsub("c\\(", "" ,listpep$posDNA, perl=TRUE)
listpep$posDNA<-gsub("\\)", "" ,listpep$posDNA, perl=TRUE)
listpep$posDNA<-gsub("\"", "" ,listpep$posDNA, perl=TRUE)
setDT(listpep)[, paste0("posDNA", 1:max(listpep$countSub)) := tstrsplit(listpep$posDNA, "\\, ")]
change<-colnames(listpep)
colnames(listpep)<-gsub("X.", "" ,change, perl=TRUE)
listpep<-do.call(cbind.data.frame, listpep)
listpep$proteinLen<-nchar(as.character(listpep$proteinseq))

## REMOVE duplicates: same peptide matched to different proteins even if different DNA
# the info regarding different codon assignment is kept

nrow(listpep[duplicated( listpep$Peptide),])
write.table(listpep[duplicated( listpep$Peptide),], "DupliPeptide.txt")

listpep_nodupli<-distinct(listpep, Peptide, Mass, RT, posDNA, .keep_all = TRUE)
nrow(listpep_nodupli[duplicated( listpep_nodupli$Peptide),])
write.table(listpep_nodupli[duplicated( listpep_nodupli$Peptide),], "DupliPepDifferentCodon.txt")

listpep_nodupli<-distinct( listpep, Peptide, Mass, RT, .keep_all = TRUE)
write.table(listpep_nodupli, file="tablefinalproteomics.txt", sep="\t", na="N/A", row.names = FALSE)
write.csv(listpep_nodupli, file="tablefinalproteomics.csv", na="N/A", row.names = FALSE)
closeAllConnections()

# check match protein-DNA by length
table<-read.delim("tablefinalproteomics.txt", header=TRUE, sep="\t")
table$proteinLen<-nchar(as.character(table$proteinseq))
DNA<-table$dnaLen
check<-(DNA==3*table$proteinLen+3)
names(check)=table$CGDid
table(check)
table$CGDid[!check]
##if needed, change the original CGD DNA file and run again from line 57##

# Count Sub codons
table<-read.table("tablefinalproteomics.txt", header = TRUE)
setDT(table)[, paste0("aaSub", 1:max(table$countSub)) := tstrsplit(table$aaSub, ",")]
table<-table[,table[,-colnames(table)[grepl("aaSub$", colnames(table))], with=FALSE],]
table<-table[,table[,-colnames(table)[grepl("posDNA$", colnames(table))], with=FALSE],]
t<-lapply(table[,colnames(table)[grepl("aaSub", colnames(table))], with=FALSE], gsub, pattern = "(\\w)\\(sub\\s(\\w)\\)", replacement = "\\1(\\2)")
t<-lapply(X = t, FUN = function(t) gsub(pattern = " ", replacement = "", x = t, fixed = TRUE))
table<-table[,table[,-colnames(table)[grepl("aaSub", colnames(table))], with=FALSE],]
t<-do.call(cbind,t)
table<-cbind(table , t)

col<-grep("aaSub", colnames(table))
col2<-grep("posDNA", colnames(table))
list<-c()
table$x<-1
for(i in 1:max(table$countSub)){
  save<-aggregate(table$x ~ table[,col[i], with=FALSE][[1]] + table[,col2[i], with=FALSE][[1]], data=table, FUN=sum)
  names(save)<-c("aaSub","posDNA","freq")
  #save$pos<-i
  list[[i]]<-save
  names(list)[[i]]<-paste0("n",i)
}
final_list<-do.call(rbind,list)
setDT(final_list)
final_list[, freq:=sum(freq), .(aaSub, posDNA)] 
final_list = unique(final_list)
final_list<-final_list[order(final_list$freq, decreasing = T),]
write.table(final_list, "TOPcodonsSubst.txt", row.names = FALSE)
