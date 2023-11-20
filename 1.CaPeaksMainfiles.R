
library(data.table)
library(stringr)
library(plyr)
library(Biostrings)
library(reshape2)
library(dplyr)
library(seqRFLP)
library(pracma)
library(gsubfn)

mainDir<-getwd()
subDir1<-"OtherFasta"
ifelse(!dir.exists(file.path(mainDir, subDir1)), dir.create(file.path(mainDir, subDir1)), FALSE)

# Original file .csv from PEAKS
mainfile <- read.csv(file="protein-peptides.csv")

countINS<-c()
countDEL<-c()
for(i in 1:length(mainfile$Peptide)){
  countINS[[i]]<-grep("\\(ins\\)", mainfile$Peptide[i])
  countDEL[[i]]<-grep("\\(del\\s\\w+\\)", mainfile$Peptide[i])
}
countINS[sapply(countINS, function(x) length(x)==0)] <- NA
countINS<- countINS[!is.na(as.vector(countINS))]
print(paste("There are", length(countINS), "insertions in this table"))
countDEL[sapply(countDEL, function(x) length(x)==0)] <- NA
countDEL<- countDEL[!is.na(as.vector(countDEL))]
print(paste("There are", length(countDEL), "deletions in this table"))

# Order by CGDid and Peptide
mainfile<-mainfile[order(mainfile$Protein.Accession, mainfile$Peptide),]
names(mainfile)[3]<-"CGDid"
mainfile$ID<-gsub("\\_[A-B]", "", mainfile$CGDid)
mainfile$ID<-type.convert(mainfile$ID)

# separate files
mainfile$countSub <- str_count(mainfile$Peptide, "\\(sub ")
extractSub<-subset(mainfile, mainfile$countSub>=1)
extractNoSub<-subset(mainfile, mainfile$countSub<1)

######GENERATE NoSub FILE FOR PEPTIDE VALIDATION######

write.table(extractNoSub, "protein-peptides_NoSubForAccept.txt", row.names = FALSE)
extractNoSubForAccept <- read.delim(file="protein-peptides_NoSubForAccept.txt", sep=" ")

# remove duplicates from different alleles
extractNoSubForAccept<-distinct(extractNoSubForAccept, ID, Peptide, Mass, X.10lgP, Length, .keep_all = TRUE)

countINS<-c()
countDEL<-c()
for(i in 1:length(extractNoSubForAccept$Peptide)){
  countINS[[i]]<-grep("\\(ins\\)", extractNoSubForAccept$Peptide[i])
  countDEL[[i]]<-grep("\\(del\\s\\w+\\)", extractNoSubForAccept$Peptide[i])
}
countINS[sapply(countINS, function(x) length(x)==0)] <- NA
countINS<- countINS[!is.na(as.vector(countINS))]
print(paste("There are", length(countINS), "insertions in this table"))
countDEL[sapply(countDEL, function(x) length(x)==0)] <- NA
countDEL<- countDEL[!is.na(as.vector(countDEL))]
print(paste("There are", length(countDEL), "deletions in this table"))

extractNoSubForAccept$pepSeq<-gsub("\\(\\+\\d+\\)", "", extractNoSubForAccept$Peptide, perl=TRUE)
extractNoSubForAccept$pepSeq<-gsub("\\(\\-\\d+\\)", "", extractNoSubForAccept$pepSeq, perl=TRUE)
extractNoSubForAccept$pepSeq<-gsub("\\(\\+\\d+\\.\\d+\\)", "", extractNoSubForAccept$pepSeq, perl=TRUE)
extractNoSubForAccept$pepSeq<-gsub("\\(\\-\\d+\\.\\d+\\)", "", extractNoSubForAccept$pepSeq, perl=TRUE)
extractNoSubForAccept$pepSeq<-gsub("\\(\\+\\.\\d+\\)", "", extractNoSubForAccept$pepSeq, perl=TRUE)
extractNoSubForAccept$pepSeq<-gsub("\\(\\-\\.\\d+\\)", "", extractNoSubForAccept$pepSeq, perl=TRUE)

extractNoSubForAccept$pepSeq2<-gsub("[A-Z]\\.(.+)\\.[A-Z]", "\\1", extractNoSubForAccept$pepSeq, perl=TRUE)
extractNoSubForAccept$pepSeq2<-gsub("^\\-\\.(.+)\\.[A-Z]", "\\1", extractNoSubForAccept$pepSeq2, perl=TRUE)
extractNoSubForAccept$pepSeq2<-gsub("[A-Z]\\.(.+)\\.\\-$", "\\1", extractNoSubForAccept$pepSeq2, perl=TRUE)
extractNoSubForAccept$pepSeq2<-gsub("^\\-\\.(.+)\\.\\-$", "\\1", extractNoSubForAccept$pepSeq2, perl=TRUE)
extractNoSubForAccept$pepSeq2<-gsub("^\\-\\.(.+)", "\\1", extractNoSubForAccept$pepSeq2, perl=TRUE)
extractNoSubForAccept$pepSeq2<-gsub("^[A-Z]\\.(.+)", "\\1", extractNoSubForAccept$pepSeq2, perl=TRUE)
extractNoSubForAccept$pepSeq2<-gsub("(.+)\\.[A-Z]$", "\\1", extractNoSubForAccept$pepSeq2, perl=TRUE)
extractNoSubForAccept$pepSeq2<-gsub("(.+)\\.\\-", "\\1", extractNoSubForAccept$pepSeq2, perl=TRUE)

extractNoSubForAccept$pepSeq3<-gsub("[A-Z]\\(ins\\)", "", extractNoSubForAccept$pepSeq2, perl=TRUE)
extractNoSubForAccept$pepSeq3<-gsub("\\(del\\s(\\w+)\\)", "\\1", extractNoSubForAccept$pepSeq3, perl=TRUE)

extractNoSubForAccept$pepSeq4<-gsub("[A-Z]\\.(.+)\\.[A-Z]", "\\1", extractNoSubForAccept$Peptide, perl=TRUE)
extractNoSubForAccept$pepSeq4<-gsub("^\\-\\.(.+)\\.[A-Z]", "\\1", extractNoSubForAccept$pepSeq4, perl=TRUE)
extractNoSubForAccept$pepSeq4<-gsub("[A-Z]\\.(.+)\\.\\-$", "\\1", extractNoSubForAccept$pepSeq4, perl=TRUE)
extractNoSubForAccept$pepSeq4<-gsub("^\\-\\.(.+)\\.\\-$", "\\1", extractNoSubForAccept$pepSeq4, perl=TRUE)
extractNoSubForAccept$pepSeq4<-gsub("^\\-\\.(.+)", "\\1", extractNoSubForAccept$pepSeq4, perl=TRUE)
extractNoSubForAccept$pepSeq4<-gsub("^[A-Z]\\.(.+)", "\\1", extractNoSubForAccept$pepSeq4, perl=TRUE)
extractNoSubForAccept$pepSeq4<-gsub("(.+)\\.[A-Z]$", "\\1", extractNoSubForAccept$pepSeq4, perl=TRUE)
extractNoSubForAccept$pepSeq4<-gsub("(.+)\\.\\-", "\\1", extractNoSubForAccept$pepSeq4, perl=TRUE)

write.table(extractNoSubForAccept, "protein-peptides_NoSubForAccept.txt", row.names = FALSE)

######GENERATE NoSub FILE FOR MISTRANSLATION FREQUENCY######

# validate InDels _ peptidesInDels.csv 
extractNoSub_noInDels<-extractNoSubForAccept[grep("\\(ins\\)|\\(del\\s(\\w+)\\)", extractNoSubForAccept$Peptide, invert = T),]
extractNoSub_InDels<-extractNoSubForAccept[grep("\\(ins\\)|\\(del\\s(\\w+)\\)", extractNoSubForAccept$Peptide),]

mainfile_pepInDels <- read.csv("peptide_InDels.csv")
extractNoSub_InDels$matchMass <- match(extractNoSub_InDels$Mass, mainfile_pepInDels$Mass, nomatch = NA_integer_)
extractNoSub_InDels$matchPep <- match(extractNoSub_InDels$pepSeq4, mainfile_pepInDels$Peptide, nomatch = NA_integer_)

extractNoSub_InDels<-filter_at(extractNoSub_InDels, .vars = vars(matchMass, matchPep), .vars_predicate = all_vars(!is.na(.)))
extractNoSub_InDels<-as.data.frame(extractNoSub_InDels)
extractNoSub_InDels<-extractNoSub_InDels[,grep("match", colnames(extractNoSub_InDels), invert = T)]

# remove low Ascore PTM peptides 
extractNoSub_noPTM<-extractNoSub_noInDels[grep("\\(", extractNoSub_noInDels$Peptide, invert = T),]
extractNoSub_PTM<-extractNoSub_noInDels[grep("\\(", extractNoSub_noInDels$Peptide),]
extractNoSub_PTM<-extractNoSub_PTM[grep(".+", extractNoSub_PTM$PTM),]

extractNoSub<-rbind(extractNoSub_noPTM, extractNoSub_PTM, extractNoSub_InDels)

countINS<-c()
countDEL<-c()
for(i in 1:length(extractNoSub$Peptide)){
  countINS[[i]]<-grep("\\(ins\\)", extractNoSub$Peptide[i])
  countDEL[[i]]<-grep("\\(del\\s\\w+\\)", extractNoSub$Peptide[i])
}
countINS[sapply(countINS, function(x) length(x)==0)] <- NA
countINS<- countINS[!is.na(as.vector(countINS))]
print(paste("There are", length(countINS), "insertions in this table"))
countDEL[sapply(countDEL, function(x) length(x)==0)] <- NA
countDEL<- countDEL[!is.na(as.vector(countDEL))]
print(paste("There are", length(countDEL), "deletions in this table"))

write.table(extractNoSub, "protein-peptides_NoSub.txt", row.names = FALSE)

######GENERATE Sub FILE######
extractSub$pepSeq<-gsub("\\(\\+\\d+\\)", "", extractSub$Peptide, perl=TRUE)
extractSub$pepSeq<-gsub("\\(\\-\\d+\\)", "", extractSub$pepSeq, perl=TRUE)
extractSub$pepSeq<-gsub("\\(\\+\\d+\\.\\d+\\)", "", extractSub$pepSeq, perl=TRUE)
extractSub$pepSeq<-gsub("\\(\\-\\d+\\.\\d+\\)", "", extractSub$pepSeq, perl=TRUE)
extractSub$pepSeq<-gsub("\\(\\+\\.\\d+\\)", "", extractSub$pepSeq, perl=TRUE)
extractSub$pepSeq<-gsub("\\(\\-\\.\\d+\\)", "", extractSub$pepSeq, perl=TRUE)

extractSub$pepSeq2<-gsub("[A-Z]\\.(.+)\\.[A-Z]", "\\1", extractSub$pepSeq, perl=TRUE)
extractSub$pepSeq2<-gsub("^\\-\\.(.+)\\.[A-Z]", "\\1", extractSub$pepSeq2, perl=TRUE)
extractSub$pepSeq2<-gsub("[A-Z]\\.(.+)\\.\\-$", "\\1", extractSub$pepSeq2, perl=TRUE)
extractSub$pepSeq2<-gsub("^\\-\\.(.+)\\.\\-$", "\\1", extractSub$pepSeq2, perl=TRUE)
extractSub$pepSeq2<-gsub("^\\-\\.(.+)", "\\1", extractSub$pepSeq2, perl=TRUE)
extractSub$pepSeq2<-gsub("^[A-Z]\\.(.+)", "\\1", extractSub$pepSeq2, perl=TRUE)
extractSub$pepSeq2<-gsub("(.+)\\.[A-Z]$", "\\1", extractSub$pepSeq2, perl=TRUE)
extractSub$pepSeq2<-gsub("(.+)\\.\\-", "\\1", extractSub$pepSeq2, perl=TRUE)

extractSub$pepSeq3<-gsub("[A-Z]\\(ins\\)", "", extractSub$pepSeq2, perl=TRUE)
extractSub$pepSeq3<-gsub("\\(del\\s(\\w+)\\)", "\\1", extractSub$pepSeq3, perl=TRUE)

extractSub$pepSeq4<-gsub("[A-Z]\\.(.+)\\.[A-Z]", "\\1", extractSub$Peptide, perl=TRUE)
extractSub$pepSeq4<-gsub("^\\-\\.(.+)\\.[A-Z]", "\\1", extractSub$pepSeq4, perl=TRUE)
extractSub$pepSeq4<-gsub("[A-Z]\\.(.+)\\.\\-$", "\\1", extractSub$pepSeq4, perl=TRUE)
extractSub$pepSeq4<-gsub("^\\-\\.(.+)\\.\\-$", "\\1", extractSub$pepSeq4, perl=TRUE)
extractSub$pepSeq4<-gsub("^\\-\\.(.+)", "\\1", extractSub$pepSeq4, perl=TRUE)
extractSub$pepSeq4<-gsub("^[A-Z]\\.(.+)", "\\1", extractSub$pepSeq4, perl=TRUE)
extractSub$pepSeq4<-gsub("(.+)\\.[A-Z]$", "\\1", extractSub$pepSeq4, perl=TRUE)
extractSub$pepSeq4<-gsub("(.+)\\.\\-", "\\1", extractSub$pepSeq4, perl=TRUE)

CGD<-c()
for (i in 1:length(extractSub$CGDid)) {
  x<-i
  CGD[i]<-paste(extractSub$CGDid[i],x,sep="|")
}
extractSub<-cbind.data.frame(CGD, extractSub)

extractSub$seqOriginal<-gsub("[A-Z]\\(ins\\)", "", extractSub$pepSeq, perl=TRUE)
extractSub$seqOriginal<-gsub("\\(del\\s(\\w+)\\)", "\\1", extractSub$seqOriginal, perl=TRUE)
extractSub$seqOriginal<-gsub("[A-Z]\\(sub\\s([A-Z])\\)", "\\1" ,extractSub$seqOriginal, perl=TRUE)
extractSub$seqMut<-gsub("\\(ins\\)", "", extractSub$pepSeq, perl=TRUE)
extractSub$seqMut<-gsub("\\(del\\s\\w+\\)", "", extractSub$seqMut, perl=TRUE)
extractSub$seqMut<-gsub("([A-Z])\\(sub\\s[A-Z]\\)", "\\1" ,extractSub$seqMut, perl=TRUE)

extractSub$seqOriginal2<-gsub("[A-Z]\\(ins\\)", "", extractSub$pepSeq2, perl=TRUE)
extractSub$seqOriginal2<-gsub("\\(del\\s(\\w+)\\)", "\\1", extractSub$seqOriginal2, perl=TRUE)
extractSub$seqOriginal2<-gsub("[A-Z]\\(sub\\s([A-Z])\\)", "\\1" ,extractSub$seqOriginal2, perl=TRUE)
extractSub$seqMut2<-gsub("\\(ins\\)", "", extractSub$pepSeq2, perl=TRUE)
extractSub$seqMut2<-gsub("\\(del\\s\\w+\\)", "", extractSub$seqMut2, perl=TRUE)
extractSub$seqMut2<-gsub("([A-Z])\\(sub\\s[A-Z]\\)", "\\1" ,extractSub$seqMut2, perl=TRUE)

extractSub$seqOriginal_PTM<-gsub("[A-Z]\\(ins\\)", "", extractSub$Peptide, perl=TRUE)
extractSub$seqOriginal_PTM<-gsub("\\(del\\s(\\w+)\\)", "\\1", extractSub$seqOriginal_PTM, perl=TRUE)
extractSub$seqOriginal_PTM<-gsub("[A-Z]\\(sub\\s([A-Z])\\)", "\\1" ,extractSub$seqOriginal_PTM, perl=TRUE)
extractSub$seqMut_PTM<-gsub("\\(ins\\)", "", extractSub$Peptide, perl=TRUE)
extractSub$seqMut_PTM<-gsub("\\(del\\s\\w+\\)", "", extractSub$seqMut_PTM, perl=TRUE)
extractSub$seqMut_PTM<-gsub("([A-Z])\\(sub\\s[A-Z]\\)", "\\1" ,extractSub$seqMut_PTM, perl=TRUE)

extractSub$aaSub<-str_extract_all(extractSub$pepSeq2, "\\w\\(sub [^()]+\\)")
extractSub$aaSub<-gsub("c\\(\"", "" ,extractSub$aaSub, perl=TRUE)
extractSub$aaSub<-gsub("\"\\)", "" ,extractSub$aaSub, perl=TRUE)
extractSub$aaSub<-gsub("\"", "" ,extractSub$aaSub, perl=TRUE)
extractSub$temp<-gsub("[A-Z]\\(sub [A-Z]\\)", "_" ,extractSub$pepSeq3, perl=TRUE)
flag <- gregexpr("_", extractSub$temp)
extractSub$posSub<-as.character(flag)
extractSub$posSub<-gsub("c\\(", "" ,extractSub$posSub, perl=TRUE)
extractSub$posSub<-gsub("\\)", "" ,extractSub$posSub, perl=TRUE)

vectorpos<-c()
for (i in 1:length(extractSub$CGDid)) {
  if (extractSub$posSub[i] %like% ":" == TRUE) {
    vectorpos<-unlist(strsplit(as.character(extractSub$posSub[i]), "[:]"))
    start <- as.numeric(vectorpos[1])
    end<-as.numeric(vectorpos[length(vectorpos)])
    a<-as.vector(seq(start, end, by=1))
    extractSub$posSub[i]<-paste(as.character(a),collapse=", ")
  }
}

setDT(extractSub)[, paste0("aaOriginal", 1:max(extractSub$countSub)) := tstrsplit(extractSub$aaSub, ",")]
t<-lapply(extractSub[,colnames(extractSub)[grepl("aaOriginal", colnames(extractSub))], with=FALSE], gsub, pattern = "\\w\\(sub\\s(\\w)\\)", replacement = "\\1")
t<-lapply(X = t, FUN = function(t) gsub(pattern = " ", replacement = "", x = t, fixed = TRUE))
t<-do.call(cbind,t)
extractSub<-extractSub[,extractSub[,-colnames(extractSub)[grepl("aaOriginal", colnames(extractSub))], with=FALSE],]
extractSub<-cbind(extractSub , t)
setDT(extractSub)[, paste0("info", 1:max(extractSub$countSub)) := tstrsplit(extractSub$posSub, ",")]
temp<-sapply(extractSub[,colnames(extractSub)[grepl("info", colnames(extractSub))], with=FALSE], as.numeric)
temp<-(temp + extractSub$Start) - 1
colnames(temp)<-gsub("info","protPos",colnames(temp))
extractSub<-cbind(extractSub,temp)

# select mutations by ion intensity_peptides.csv 
mainfile_pep <- read.csv("peptide.csv")

extractSub$matchMass <- match(extractSub$Mass, mainfile_pep$Mass, nomatch = NA_integer_)
extractSub$matchPep <- match(extractSub$pepSeq4, mainfile_pep$Peptide, nomatch = NA_integer_)

# remove rows with NAs in matchMass and matchPep
extractSub_valPep<-filter_at(extractSub, .vars = vars(matchMass, matchPep), .vars_predicate = all_vars(!is.na(.)))
extractSub_valPep<-as.data.frame(extractSub_valPep)

# remove duplicates from different alleles
extractSub_nodupli<-distinct(extractSub_valPep, ID, Peptide, Mass, X.10lgP, Length, .keep_all = TRUE)

countINS<-c()
countDEL<-c()
for(i in 1:length(extractSub_nodupli$Peptide)){
  countINS[[i]]<-grep("\\(ins\\)", extractSub_nodupli$Peptide[i])
  countDEL[[i]]<-grep("\\(del\\s\\w+\\)", extractSub_nodupli$Peptide[i])
}
countINS[sapply(countINS, function(x) length(x)==0)] <- NA
countINS<- countINS[!is.na(as.vector(countINS))]
print(paste("There are", length(countINS), "insertions in this table!"))
countDEL[sapply(countDEL, function(x) length(x)==0)] <- NA
countDEL<- countDEL[!is.na(as.vector(countDEL))]
print(paste("There are", length(countDEL), "deletions in this table!"))

pep_original <- data.frame(extractSub_nodupli$CGD,extractSub_nodupli$seqOriginal)
pep_original.fasta = dataframe2fas(pep_original, file="OtherFasta/pep_original.fasta")
pep_mutate <- data.frame(extractSub_nodupli$CGD,extractSub_nodupli$seqMut)
pep_mutate.fasta = dataframe2fas(pep_mutate, file="OtherFasta/pep_mutate.fasta")

######UPLOAD ORIGINAL PROTEIN SEQUENCE######
ALLproteinseq<-readAAStringSet(file="Ca_A22-s07-m01-r149_AT_PROTEIN_29032022.fasta")
seq_name = names(ALLproteinseq)
sequence = paste(ALLproteinseq)
proteins <- data.frame(seq_name, sequence)
names(proteins)[1]<-"CGDid"
names(proteins)[2]<-"proteinseq"
proteins<-proteins[order(proteins$CGDid),]
finalproteins<-merge(extractSub_nodupli,proteins, by="CGDid", all.x=TRUE)
CGD<-c()
for (i in 1:length(finalproteins$CGDid)) {
  x<-i
  CGD[i]<-paste(finalproteins$CGDid[i],x,sep="|")
}
finalproteins<-cbind.data.frame(CGD, finalproteins)
finalproteins<-finalproteins[,c(1,55)] ##check proteinseq column number
originalprotein.fasta = dataframe2fas(finalproteins, file="OtherFasta/originalprotein.fasta")
ALL_merged<-merge(extractSub_nodupli,proteins, by="CGDid", all.x=TRUE)
ALL_merged<-ALL_merged[,c(2,3,4,1,5:54)] ##check total column number
write.table(ALL_merged, "protein-peptides_Sub.txt", row.names = FALSE)
