
library(seqinr)
library(seqRFLP)
library(stringr)
library(data.table)
library(plyr)
library(dplyr)

mainDir<-getwd()
subDir1<-"F1"
subDir2<-"F2"
ifelse(!dir.exists(file.path(mainDir, subDir1)), dir.create(file.path(mainDir, subDir1)), FALSE)
ifelse(!dir.exists(file.path(mainDir, subDir2)), dir.create(file.path(mainDir, subDir2)), FALSE)

# open mainfiles
pepsubst <- read.delim("protein-peptides_Sub.txt", sep=" ")
pepsubst<-pepsubst[order(pepsubst$CGDid),]

mainfile_noSub <- read.delim(file="protein-peptides_NoSubForAccept.txt", sep=" ")
mainfile_noSub$indexing<-rownames(mainfile_noSub)
col<-grep("protPos", colnames(pepsubst))
col2<-grep("aaOriginal", colnames(pepsubst))

# Check allele issue (peptides from different alleles or paralogs, misassigned as mutated)
tempi<-c()
for(i in 1:nrow(pepsubst)){
  #print(i)
  tempi<-mainfile_noSub[grepl(pepsubst$seqMut[i], mainfile_noSub$pepSeq),]
  if(nrow(tempi)>0) {
    pepsubst$found[i]<-"DECLINE"
    pepsubst$index[i]<-paste(tempi$indexing, collapse="|")
  }
  else{
    pepsubst$found[i]<-"ACCEPT"
    pepsubst$index[i]<-NA
  }
}
pepsubst_check<-pepsubst[pepsubst$found=="ACCEPT",]
table(pepsubst_check$found)

pepsubst_check$found<-NULL

# Check wild type peptides for the corresponding mutations (same wt peptide sequence but may have different PTMs)
tempi<-c()
for(i in 1:nrow(pepsubst_check)){
  #print(i)
  tempi<-mainfile_noSub[grepl(pepsubst_check$seqOriginal[i], mainfile_noSub$pepSeq) & pepsubst_check$CGDid[i] == mainfile_noSub$CGDid,]
  if(nrow(tempi)>0) {
    pepsubst_check$found[i]<-"ACCEPT"
    pepsubst_check$index[i]<-paste(tempi$indexing, collapse="|")
  }
  else{
    pepsubst_check$found[i]<-"DECLINE"
    pepsubst_check$index[i]<-NA
  }
}
pepsubst_Accepted<-pepsubst_check[pepsubst_check$found=="ACCEPT",]
table(pepsubst_Accepted$found)

# Generate complete mainfile: pepsubst accepted + mainfile nosub (with filtered indels and PTMs)
write.csv(pepsubst_Accepted, file="protein-peptides_SubAccepted.txt", row.names = FALSE)
mainfile_noSub <- read.delim(file="protein-peptides_NoSub.txt", sep=" ")
allcodons_prep<-rbind.fill(pepsubst_Accepted, mainfile_noSub)
write.csv(allcodons_prep, file="protein-peptidesAllCodons.txt", row.names = FALSE)

### EXTRA FILTERS

# F1 filtering =pepseq3 =CGDid =mass
pepsubstF1<-distinct(pepsubst, CGDid, pepSeq3, Mass, .keep_all = TRUE)
pepsubst_AcceptedF1<-distinct(pepsubst_Accepted, CGDid, pepSeq3, Mass, .keep_all = TRUE)
mainfile_noSubF1<-distinct(mainfile_noSub, CGDid, pepSeq3, Mass, .keep_all = TRUE)

# F2 filtering =pepseq3 =CGDid
pepsubstF2<-distinct(pepsubst, CGDid, pepSeq3, .keep_all = TRUE)
pepsubst_AcceptedF2<-distinct(pepsubst_Accepted, CGDid, pepSeq3, .keep_all = TRUE)
mainfile_noSubF2<-distinct(mainfile_noSub, CGDid, pepSeq3, .keep_all = TRUE)

# pepsubst accepted F1 + mainfile nosub F1
write.csv(pepsubst_AcceptedF1, file="F1/protein-peptides_SubAcceptedF1.txt", row.names = FALSE)
allcodons_prepF1<-rbind.fill(pepsubst_AcceptedF1, mainfile_noSubF1)
write.csv(allcodons_prepF1, file="F1/protein-peptidesAllCodonsF1.txt", row.names = FALSE)

# pepsubst accepted F2 + mainfile nosub F2
write.csv(pepsubst_AcceptedF2, file="F2/protein-peptides_SubAcceptedF2.txt", row.names = FALSE)
allcodons_prepF2<-rbind.fill(pepsubst_AcceptedF2, mainfile_noSubF2)
write.csv(allcodons_prepF2, file="F2/protein-peptidesAllCodonsF2.txt", row.names = FALSE)
