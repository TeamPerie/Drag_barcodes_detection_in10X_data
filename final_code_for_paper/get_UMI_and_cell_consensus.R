rm(list = ls())
options(echo=TRUE)
args = commandArgs(trailingOnly=TRUE)
print(args)
wd = as.character(args[1])
infile = as.character(args[2])
outfile = as.character(args[3])
outfile2 = as.character(args[4])
bc.length = as.numeric(args[5])
setwd(wd)

library(seqinr)

#script to read in all barcodes and find consensus
#Assumes barcode is bc.length bp long
vbcs = read.table(infile, header = F, sep = "\t", stringsAsFactors = F,
                  colClasses = c(rep("character",3),rep("numeric",2),rep("character",3)))

#get unique 10x cell barcode and UMI combinations
BC.UMI.uni = unique(vbcs[,1:2])

con.seq = character(nrow(BC.UMI.uni)) #consensus seq for each UMI
prop.con = character(nrow(BC.UMI.uni)) #proportion of supporting reads for each base
obs = character(nrow(BC.UMI.uni)) #number of reads covering each base
side = character(nrow(BC.UMI.uni)) #side of main flanking region match
for(i in 1:nrow(BC.UMI.uni)){
  
  #extract relevant rows
  rel.rows = vbcs[which(vbcs[,1]==BC.UMI.uni[i,1] & vbcs[,2]==BC.UMI.uni[i,2]),]
  
  #take one sequence for each read 
  #(reads may be there twice if grep match was found on both sides of barcode)
  reads.uni = unique(rel.rows[,6])
  rel.rows2 = do.call("rbind",lapply(unique(rel.rows[,6]), function(x) {
    if(nrow(rel.rows[rel.rows[,6]==x,])==1){
      return(rel.rows[rel.rows[,6]==x,])
    }else{
      rel.rows2 = rel.rows[rel.rows[,6]==x,]
      #arbitrarily return first row 
      return(rel.rows2[1,])
    }
  }))
  
  #organise sequences into columns with correct bases in vertical columns
  vbc.mat = apply(rel.rows2, 1, function(x) {
    seq.c = s2c(x[3])
    seq.l = length(seq.c)
    if(x[8]=="R"){
      return(c(rep(NA,bc.length-seq.l),seq.c))
    }else{
      return(c(seq.c,rep(NA,bc.length-seq.l)))
    }
  })
  
  #get consensus for each row
  con.seq.c = apply(vbc.mat, 1, function(x) {
    ux = unique(x[!is.na(x)])
    return(ux[which.max(tabulate(match(x, ux)))])
  })
  con.seq[i] = paste0(con.seq.c[!is.na(con.seq.c)], collapse = "")
  
  #get proportion of reads which support consensus
  con.prop.c = apply(vbc.mat, 1, function(x) {
    ux = unique(x[!is.na(x)])
    return(max(tabulate(match(x, ux)))/length(x[!is.na(x)]))
  })
  prop.con[i] = paste0(signif(con.prop.c[!is.na(con.prop.c)],digits = 3), collapse = ":")
  
  #get number of reads each base observed with
  obs.c = apply(vbc.mat, 1, function(x) {
    return(length(x[!is.na(x)]))
  })
  obs[i] = paste0(obs.c, collapse = ":")
  
  #get side of barcode
  if(nchar(con.seq[i])==bc.length){
    side[i] = "B"
  }else{
    side[i] = rel.rows[1,8]
  }
}

output.tab = cbind(BC.UMI.uni, con.seq, obs, prop.con, side, stringsAsFactors=F)

#filter UMIs which do not have good consensus across reads
kept = apply(output.tab, 1, function(x) {
  if(sum(as.numeric(unlist(strsplit(x[5],":",fixed=T)))<=0.5)<=3){
    return(1)
  }else{
    return(0)
  }
})
output.tab = cbind(output.tab, kept)

write.table(output.tab, file = outfile, sep = ",", quote = F,
            col.names = c("BC_10x","UMI_10x","cons_BC_lenti","Nread_per_base",
                          "prop_match_cons","side","kept"), row.names = F)

#keep only UMIs with good consensus
output.tab = output.tab[output.tab[,7]==1,]

#same approach as above, but now get cell consensus across UMIs
BC.uni = unique(output.tab[,1])

con.seq = character(length(BC.uni))
prop.con = character(length(BC.uni))
obs = character(length(BC.uni))
side = character(length(BC.uni))
for(i in 1:length(BC.uni)){
  
  #extract relevant rows
  rel.rows = output.tab[which(output.tab[,1]==BC.uni[i]),]
  
  #organise sequences into columns with correct bases in vertical columns
  vbc.mat = apply(rel.rows, 1, function(x) {
    seq.c = s2c(x[3])
    seq.l = length(seq.c)
    if(x[6]=="R"){
      return(c(rep(NA,bc.length-seq.l),seq.c))
    }else if(x[6]=="B"){
      return(seq.c)
    }else{
      return(c(seq.c,rep(NA,bc.length-seq.l)))
    }
  })
  
  #get consensus for each row
  con.seq.c = apply(vbc.mat, 1, function(x) {
    ux = unique(x[!is.na(x)])
    return(ux[which.max(tabulate(match(x, ux)))])
  })
  con.seq[i] = paste0(con.seq.c[!is.na(con.seq.c)], collapse = "")
  
  #get proportion of reads which support consensus
  con.prop.c = apply(vbc.mat, 1, function(x) {
    ux = unique(x[!is.na(x)])
    return(max(tabulate(match(x, ux)))/length(x[!is.na(x)]))
  })
  prop.con[i] = paste0(signif(con.prop.c[!is.na(con.prop.c)],digits=3), collapse = ":")
  
  #get number of reads each base observed with
  obs.c = apply(vbc.mat, 1, function(x) {
    return(length(x[!is.na(x)]))
  })
  obs[i] = paste0(obs.c, collapse = ":")
  
  #get side of barcode
  if(nchar(con.seq[i])==bc.length){
    side[i] = "B"
  }else{
    side[i] = rel.rows[1,6]
  }
}

output.tab2 = cbind(BC.uni, con.seq, obs, prop.con, side)

#add column indicating barcodes which are filtered
#Filter cells which don't have good consensus aross UMIs
kept = apply(output.tab2, 1, function(x) {
  if(sum(as.numeric(unlist(strsplit(x[4],":",fixed=T)))<=0.5)<=3){
    return(1)
  }else{
    return(0)
  }
})
output.tab2 = cbind(output.tab2, kept)

write.table(output.tab2, file = outfile2, sep = ",", quote = F,
            col.names = c("BC_10x","cons_BC_lenti","Numi_per_base",
                          "Prop_match_cons","side","kept"), row.names = F)


