
rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
wd = as.character(args[1])
infile = as.character(args[2]) #reads file, pre-grepped for main pattern match
outfile = as.character(args[3]) #file to print results to
bc.length = as.numeric(args[4]) #barcode length = 20
n.match = as.numeric(args[5]) #number of bp for main match
n.match.opp = as.numeric(args[6]) #number of bp for match on other side
fl.R = as.character(args[7]) #flanking sequence R
fl.L = as.character(args[8]) #flanking sequence L
full.length = as.numeric(args[9]) #1=take only full length barcodes

setwd(wd)

library(Biostrings)
library(seqinr)

#read in file with cell barcode, UMI, read name, sequence etc
nr = length(count.fields(infile, sep="\t"))
maxfields = max(count.fields(infile, sep="\t"))
reads = read.table(infile, nrows = nr, comment.char="", colClasses = rep("character",maxfields), 
                   fill=T, header=F, sep="\t", stringsAsFactors = F)

#make sequence
fl.s = s2c(fl.L)

#find rows which have a match to flanking sequence
#one mismatch allowed, no indels, can change below
grep.list = lapply(reads[,2], 
       function(x) matchPattern(fl.L, x, min.mismatch=0, 
                                max.mismatch=1, with.indels=F))
logi = sapply(grep.list, function(x) length(x)>0)

if(sum(logi)>0){
  #extract reads which have a match
  reads.to.keep = reads[logi,]
  grep.to.keep = grep.list[logi]
  #get start and end position of match
  starts = sapply(grep.to.keep, function(x) start(x)[1])
  ends = sapply(grep.to.keep, function(x) end(x)[1])
  reads.to.keep = cbind(reads.to.keep,starts,ends)
  #check that start and end base are same as reference i.e. that edit is in the middle somewhere if it exists
  logi2 = apply(reads.to.keep, 1, function(x) {
    seq.s = s2c(as.character(x[2]))
    if(as.numeric(x[5])==0){
      return(T)
    }else if(as.numeric(x[6])==(length(seq.s)+1)){
      return(F)
    }else{
      return( fl.s[1]==seq.s[as.numeric(x[5])] && fl.s[length(fl.s)]==seq.s[as.numeric(x[6])] )
    }
  })
  reads.to.keep = reads.to.keep[logi2,]
  
  #extract part of read which goes into barcode
  if(nrow(reads.to.keep)>0){
    vbcs = character(nrow(reads.to.keep))
    for(i in 1:nrow(reads.to.keep)){
      if((reads.to.keep[i,ncol(reads.to.keep)]+bc.length)<=nchar(reads.to.keep[i,2])){
        vbcs[i] = substr(reads.to.keep[i,2],(reads.to.keep[i,ncol(reads.to.keep)]+1),
                         (reads.to.keep[i,ncol(reads.to.keep)]+bc.length))
      }else{
        vbcs[i] = substr(reads.to.keep[i,2],(reads.to.keep[i,ncol(reads.to.keep)]+1),
                         nchar(reads.to.keep[i,2]))
      }
    }
    bcs = cbind(reads.to.keep[,1],reads.to.keep[,4],vbcs,reads.to.keep[,c(5,6,3)],
                reads.to.keep[,2], stringsAsFactors = F)
    
    if(n.match.opp>0){
      #barcodes must be full length
      #remove barcodes which do not have 4 bp exact match on rhs of barcode
      bcs = bcs[(bcs[,5] + bc.length + 1 + n.match.opp)<=nchar(bcs[,7]),]
      st.read = bcs[,5] + bc.length + 2
      en.read = bcs[,5] + bc.length + 1 + n.match.opp
      bcs = bcs[substr(bcs[,7],st.read,en.read)==substr(fl.R,1,n.match.opp),]
    }else{
      if(full.length==1){
        #remove barcodes which are not 20 base pairs
        bcs = bcs[nchar(bcs[,3])>=bc.length,]
      }
    }
    
    #write output of viral barcodes, 10x barcode, UMI etc
    if(nrow(bcs)==1){
      write.table(bcs, file=outfile, quote=FALSE, row.names=FALSE, 
                  col.names=FALSE, na="NA", sep="\t")
    }else{
      #sort by cell barcode, then by UMI
      bcs = bcs[order(bcs[,1],bcs[,2]),]
      write.table(bcs, file=outfile, quote=FALSE, row.names=FALSE, 
                  col.names=FALSE, na="NA", sep="\t")
    }
  }
}




