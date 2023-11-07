
library(stringr)
library(data.table)
library(foreach)
library(iterators)
library(protr)
library(pracma)
library(Biostrings)
library(data.table)
library(seqRFLP)
library(IRanges)
library(plyr)
library(dplyr)
library(annotate)
library(filesstrings)

table<-read.table("tablefinalproteomics.txt", header = TRUE)
mainDir<-getwd()

###WHICH AA YOU ARE LOOKING FOR?###
aa<-c("S") #change here for the aa of interest
#aa<-c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "T", "V", "W", "Y") #all other amino acids
for (x in 1:length(aa)){
  print(x)
  subDir<-aa[x]
  ifelse(!dir.exists(file.path(mainDir, subDir)), dir.create(file.path(mainDir, subDir)), FALSE)
  table$countAA <- str_count(table$seqOriginal2, aa[x])
  if(sum(table$countAA)==0){
    next
  }
  else{
    extractAA<-subset(table, table$countAA>=1)
    n <- 3
    ID=c()
    codon<-c()
    for (i in 1:length(extractAA$CGDid)) {
      pos<-as.integer(gregexpr(aa[x], extractAA$seqOriginal2[i], perl = TRUE)[[1]])
      protPos<-(extractAA$Start[i] + pos) - 1
      codonstemp<-sapply(seq(1,nchar(extractAA$dnaseq[i]),by=n), function(x) substr(extractAA$dnaseq[i], x, x+n-1))
      codon[[i]]<-codonstemp[protPos]
    }
    extractAA$posDNA<-as.character(codon)
    extractAA$posDNA<-gsub("c\\(", "" ,extractAA$posDNA, perl=TRUE)
    extractAA$posDNA<-gsub("\\)", "" ,extractAA$posDNA, perl=TRUE)
    extractAA$posDNA<-gsub("\"", "" ,extractAA$posDNA, perl=TRUE)
    max<-max(extractAA$countAA)
    setDT(extractAA)[, paste0("posDNA", 1:max):=tstrsplit(extractAA$posDNA, "\\, ")]
    
    # CHECKING NAs or NNNs 
    check<-c()
    for (i in 1:length(extractAA$CGDid)) {
      if(extractAA$posDNA1[i]=="NA" || extractAA$posDNA1[i]=="NNN") {
        print(i)
        check[i]<-i
      }
    }
    check<-check[!is.na(check)]
    print(paste0("There are ", length(check), " peptides with NAs OR NNNs in the DNA first codon (posDNA1)"))
    print(as.character(extractAA$CGDid[as.integer(check)]))
    
    write.table(extractAA, file=paste0(aa[x],"/","tablefinalproteomicsDupli_with_",aa[x],".txt"), sep="\t", na="N/A", row.names = FALSE)
    
    # REMOVE duplicates: same peptide matched to different proteins
    ##the info regarding different codon assignment is kept
    
    print(nrow(extractAA[duplicated( extractAA$Peptide),]))
    write.table(extractAA[duplicated( extractAA$Peptide),], file=paste0(aa[x],"/", "DupliPeptide.txt"), row.names = FALSE)
    
    extractAA_nodupli<-distinct(extractAA, Peptide, Mass, RT, posDNA, .keep_all = TRUE)
    print(nrow(extractAA_nodupli[duplicated( extractAA_nodupli$Peptide),]))
    write.table(extractAA_nodupli[duplicated( extractAA_nodupli$Peptide),], file=paste0(aa[x],"/", "DupliPepDifferentCodon.txt"), row.names = FALSE)
    
    extractAA_nodupli<-distinct( extractAA, Peptide, Mass, RT, .keep_all = TRUE)
    write.table(extractAA_nodupli, file=paste0(aa[x],"/","tablefinalproteomics_with_",aa[x],".txt"), sep="\t", na="N/A", row.names = FALSE)
  }
}

rm(list=ls())

##### Create frequency table

###WHICH AA YOU ARE LOOKING FOR?###
aa<-c("S") #change here for the aa of interest
#aa<-c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "T", "V", "W", "Y") #all other amino acids
for (x in 1:length(aa)){
  list<-c()
  table<-read.table(file=paste0(aa[x],"/","tablefinalproteomics_with_",aa[x],".txt"), header = TRUE, sep="\t")
  table<-table[colnames(table)[-grep("posDNA$", colnames(table))]]
  table$c<-1
  for(i in 1:max(table$countAA)){
    col<-grep("posDNA", colnames(table))
    save<-aggregate(table$c ~ table[,col[i]], data=table, FUN=sum)
    names(save)<-c("posDNA","freq")
    list[[i]]<-save
  }
  final_list<-do.call(rbind,list)
  setDT(final_list)
  final_list[, freq:=sum(freq), .(posDNA)] 
  final_list = unique(final_list)
  final_list<-final_list[!grepl("N\\/A", final_list$posDNA),]
  final_list<-final_list[order(final_list$freq, decreasing = T),]
  write.table(final_list, paste0(aa[x],"/","ALLcounts_final_",aa[x],".txt"), row.names = FALSE)
}
