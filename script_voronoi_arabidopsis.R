require(deldir)
require(nloptr)
require(DEoptim)
require(gpclib)
require(sp)
require(tripack)
require(geometry)
require(fdrtool)



setwd('/mnt/1800228E002272C4/Dropbox/Work/WangWenyi')



mapman = read.csv('mapping_arabidopsis/ath_Araport11_2017-03-14_mapping.txt', header=T, sep="\t",quote="") # Arabidopsis
tmp = grepl("not assigned.unknown",mapman[,2])
mapman = mapman[tmp==F,]
tmp_atg = mapman[,3]
mapman = cbind(mapman,tmp_atg)


mapman[,1] = paste(mapman[,1], ".", sep="")
mapman[,1] = paste("R.",mapman[,1], sep="")
mapman[,1] = paste("bin_",mapman[,1], sep="")
mapman = rbind(mapman[1,], mapman)
mapman$NAME = as.vector(mapman$NAME)
mapman[1,1] = "bin_R."
mapman[1,2] = "root"

rm(list=ls(pattern="tmp"))

### turn into a list

tmp = unique(mapman[,1])

tmp_f1= function(x, y) { grep(x,substr(y,1,nchar(as.character(x))), fixed=T)}


bins = sapply(tmp, function(x) mapman[tmp_f1(x,mapman[,1]),3] )
bins = lapply(bins,unique)
bins = lapply(bins, as.vector)
bins = lapply(bins, function(x) x[x != ""])


bins.description = as.vector(unlist(sapply(names(bins),function(x) mapman[mapman[,1] == x,2][1])))
names(bins.description) = names(bins)
rm(list=ls(pattern="tmp"))





### clean
rm(mapman)
rm(list=ls(pattern="tmp"))

### add hand edited short names
bins.description2 = bins.description
# load('shortnames')
# names(shortnames) = gsub("bin_", "bin_R.", names(shortnames), fixed = T)
# names(shortnames) = gsub("bin_R.R.", "bin_R.", names(shortnames), fixed = T)
# shortnames = shortnames[names(bins.description2)]
# 
# tmp = cbind(names(bins.description2), bins.description2, shortnames)
# write.table(tmp, "shortnames_edit.txt", sep = "\t")
#manual edit
tmp = read.csv("shortnames_edit.txt", header = T, row.names = 1,sep = "\t", quote = "")
shortnames = as.vector(tmp$shortnames)
names(shortnames) = names(bins.description2)


simpleCap <- function(x) {
  s <- substr(x, 1,1)
  s2 <- substr(x, 2,nchar(x))
  paste(toupper(s), s2, sep="", collapse=" ")
}

shortnames2 = sapply(shortnames, simpleCap)

f1 = function(x,y){
  tmp = strsplit(x,split=".", fixed=T)[[1]]
  tmp[length(tmp)] = y
  tmp = paste(tmp, collapse=".")
  return(tmp)
}
for(i in 1:length(bins.description2)){
  bins.description2[i] = f1(bins.description[i], shortnames2[i])
}



### create list of terms
# O1 <- c("bin_R." = list(unique(unlist(bins))),bins)
# #names(O1) <- c("R","T11","T12")
# O1 = O1[c(1,3:length(O1))] # remove the empty root
# O1 = lapply(O1, function(x) toupper(substr(x[grep("solyc", x)], 1,14)))
O1 = bins[lapply(O1,length) > 0]

### load the data


data1 = read.csv("data_arabidopsis/comparsion 1.csv", header = T, row.names = 1, sep = "\t")
data2 = read.csv("data_arabidopsis/comparsion 2.csv", header = T, row.names = 1, sep = "\t")
data3 = read.csv("data_arabidopsis/comparsion 3.csv", header = T, row.names = 1, sep = "\t")
data4 = read.csv("data_arabidopsis/comparsion 4.csv", header = T, row.names = 1, sep = "\t")

Y.names = unique(c(rownames(data1),rownames(data2),rownames(data3),rownames(data4)))
Y = matrix(0, nrow = length(Y.names), ncol = 4, dimnames = list(Y.names, paste("comparison",1:4, sep = "_")))
Y[rownames(data1),1] = as.numeric(as.vector(data1[,1]))
Y[rownames(data2),2] = as.numeric(as.vector(data2[,1]))
Y[rownames(data3),3] = as.numeric(as.vector(data3[,1]))
Y[rownames(data4),4] = as.numeric(as.vector(data4[,1]))

Y[is.na(Y)] = 0
Y[Y == Inf] = max(Y[Y != Inf]) + 1
Y[Y == -Inf] = min(Y[Y != -Inf]) - 1

### filter the hierarchy

O2 = lapply(O1, function(x) x[x %in% rownames(Y)]) 
O2 = O2[unlist(lapply(O2, length)) > 0]


### create vtm file
f.vtm = function(Y,Ycol, O2,bins.description){
  
  #nodeId
  nodeId = names(O2)
  
  #parentId
  tmp = sapply(names(O2), function(x) max(gregexpr(pattern ='.',x, fixed=T)[[1]][1:(length(gregexpr(pattern ='.',x, fixed=T)[[1]]) - 1)]   ))
  parentId = substr(names(O2), 1, tmp)
  parentId[1] = "project"
  parentId[parentId == ""] = "bin_R."
  parentId[nodeId == parentId] = "bin_R."
  
  #names
  namesId = bins.description[nodeId]
  
  #weights
  weight = unlist(lapply(O2, function(x) sum(Y[x,1])))
  weight = weight/max(weight)
  
  # input table
  input_table = data.frame("nodeId" = nodeId, "parentId" = parentId, "name" = namesId, "weight" = weight)
  
  rm(list=ls(pattern="tmp"))
  
  ### renormalize
  
  Yn = apply(Y,2,function(x) (x/sum(x)) * 10000)
  weightn =  unlist(lapply(O2, function(x) sum(Yn[x,Ycol])))
  
  
  # identify the level
  tmp = sapply(names(O2), function(x) length(gregexpr(pattern ='.',x, fixed=T)[[1]]))
  tmp2 = sapply(names(O2), function(x) min(gregexpr(pattern ='.',x, fixed=T)[[1]]))
  tmp[1] = 0
  
  #childId
  childNr = sapply(names(O2), function(x) length(grep(x,names(O2), fixed=T)))
  
  # generate new leaves in the vtm according to the weight of the category
  input_table2 = cbind(input_table, childNr)
  input_table2[,4] = weightn
  
  # give the root the child number
  input_table2["bin_R.",5] = nrow(input_table2)
  
  
  # add new leaves according to the weight
  
  f.newleaves = function(x){
    if(as.numeric(x[5]) == 1){
      tmp = round(as.numeric(x[4]))
      tmp = as.numeric(tmp)
      x2 = as.vector(t(x))
      if(tmp > 0){
        tmp_mat = data.frame("nodeId"=paste(rep(x2[1],tmp), "synth", 1:tmp,".", sep=""),"parendId"=rep(x2[1],tmp),"name"=rep(x2[3],tmp),"weight"=rep(1,tmp),"childNr"=rep(1,tmp))
      }else{
        tmp_mat = NULL
      } 
    }else{
      tmp_mat = NULL
    }
    return(tmp_mat)
  }
  
  #f.newleaves(input_table2['bin_8.2.9.',])	      
  
  new_leaves = apply(input_table2,1,f.newleaves)
  new_leaves = do.call(rbind, new_leaves)
  rownames(new_leaves) = new_leaves[["nodeId"]]
  colnames(new_leaves) = colnames(input_table2)
  
  
  
  input_table3 = rbind(input_table2, new_leaves)
  input_table3 = input_table3[,1:4]
  
  
  # reorder the table accordint to the hierarchy
  
  tmp0 = rownames(input_table3)
  tmp0 = gsub("bin_", "", tmp0)
  tmp0 = gsub("R", "", tmp0)
  tmp0 = gsub("synth", "", tmp0)
  tmp = sapply(tmp0, function(x)  strsplit(x, split=".",fixed=T)[[1]])
  tmp2 = max(unlist(lapply(tmp, length)))
  tmp3=c()
  for(i in 1:tmp2){
    tmp3 = cbind(tmp3, unlist(lapply(tmp, function(y) y[i])))
  }
  tmp3[is.na(tmp3)] = 0
  tmp4 = matrix(as.numeric(tmp3), ncol=ncol(tmp3))
  rownames(tmp4) = rownames(input_table3)
  for(i in ncol(tmp4):1){
    tmp4 = tmp4[order(tmp4[,i]),]
  }
  tmp4 = tmp4[c(nrow(tmp4), 1:(nrow(tmp4)-1)),]
  
  
  input_table3 = input_table3[rownames(tmp4),]
  
  
  
  
  # recalculate the weights
  tmp = sapply(rownames(input_table3), function(x) sum( grepl(x,rownames(input_table3), fixed=T) & grepl("synth",rownames(input_table3), fixed=T))) 	
  input_table3[["weight_approx"]] = tmp
  input_table3[["weight_approx"]][1] = sum(grepl("synth",rownames(input_table3), fixed=T))
  input_table3 = input_table3[as.vector(input_table3[["weight_approx"]]) > 0,]
  # write the result
  
  write.table(input_table3, 'vtm.txt', sep=";", row.names=F)
  
  # run the VTM
  system(paste("java -jar Voronoi/Voronoi-Treemap-Library/build/libs/JVoroTreemap.jar", paste(getwd(), '/vtm.txt', sep=""), sep=" "))
  
  
  # retrive the result 
  vtm = read.table('vtm-finished.txt', header=T, sep=";")
  
  #remove the file
  system(paste("rm", paste(getwd(), "/vtm-finished.txt", sep="")))
  
  # fill the names and bin info
  vtm = vtm[is.na(vtm[["nodeId"]]) == F & vtm[["nodeId"]] != 0,]
  vtm[["nodeId"]] = as.vector(input_table3[["nodeId"]])[as.vector(vtm[["nodeId"]])]
  vtm[["parentID"]] = as.vector(input_table3[as.vector(vtm[["nodeId"]]), "parentId"])
  vtm[["name"]] = as.vector(input_table3[as.vector(vtm[["nodeId"]]), "name"])
  vtm[["weight"]] = as.vector(input_table3[as.vector(vtm[["nodeId"]]), "weight"])
  
  # add the ratio between the number of protein and the sum of intensities
  tmp = unlist(lapply(O2, length))
  tmp2 = length(unique(unlist(O2)))
  tmp3 = as.vector(vtm[["weight"]]) / (tmp[as.vector(vtm[["nodeId"]])] * mean(Yn))
  vtm[["intensity"]] = tmp3
  
  return(vtm)
}

bins.description0 = bins.description
bins.description = bins.description2
Ynames = rownames(Y)


f.vtm.nonweight = function(Ynames,O2,bins.description){
  O2 = O2[lapply(O2, length) > 0]
  names(O2)[1] = "bin_R."
  
  #nodeId
  nodeId = names(O2)
  
  #parentId
  tmp = sapply(names(O2), function(x) max(gregexpr(pattern ='.',x, fixed=T)[[1]][1:(length(gregexpr(pattern ='.',x, fixed=T)[[1]]) - 1)]   ))
  parentId = substr(names(O2), 1, tmp)
  parentId[1] = "project"
  parentId[parentId == ""] = "bin_R."
  parentId[nodeId == parentId] = "bin_R."
  
  #names
  namesId = bins.description[nodeId]
  
  #weights
  weight = unlist(lapply(O2, length))
  weight = weight/max(weight)
  
  # input table
  input_table = data.frame("nodeId" = nodeId, "parentId" = parentId, "name" = namesId, "weight" = weight)
  
  rm(list=ls(pattern="tmp"))
  
  ### renormalize
  
  Yn = Y#apply(Y,2,function(x) (x/sum(x)) * 10000)
  weightn =  unlist(lapply(O2, length))
  
  
  # identify the level
  tmp = sapply(names(O2), function(x) length(gregexpr(pattern ='.',x, fixed=T)[[1]]))
  tmp2 = sapply(names(O2), function(x) min(gregexpr(pattern ='.',x, fixed=T)[[1]]))
  tmp[1] = 0
  
  #childId
  childNr = sapply(names(O2), function(x) length(grep(x,names(O2), fixed=T)))
  
  # generate new leaves in the vtm according to the weight of the category
  input_table2 = cbind(input_table, childNr)
  input_table2[,4] = weightn
  
  # give the root the child number
  input_table2["bin_R.",5] = nrow(input_table2)
  
  
  # add new leaves according to the weight
  
  f.newleaves = function(x){
    if(as.numeric(x[5]) == 1){
      tmp = round(as.numeric(x[4]))
      tmp = as.numeric(tmp)
      x2 = as.vector(t(x))
      if(tmp > 0){
        tmp_mat = data.frame("nodeId"=paste(rep(x2[1],tmp), "synth", 1:tmp,".", sep=""),"parendId"=rep(x2[1],tmp),"name"=rep(x2[3],tmp),"weight"=rep(1,tmp),"childNr"=rep(1,tmp))
      }else{
        tmp_mat = NULL
      } 
    }else{
      tmp_mat = NULL
    }
    return(tmp_mat)
  }
  
  #f.newleaves(input_table2['bin_8.2.9.',])	      
  
  new_leaves = apply(input_table2,1,f.newleaves)
  new_leaves = do.call(rbind, new_leaves)
  rownames(new_leaves) = new_leaves[["nodeId"]]
  colnames(new_leaves) = colnames(input_table2)
  
  
  
  input_table3 = rbind(input_table2, new_leaves)
  input_table3 = input_table3[,1:4]
  
  
  # reorder the table according to the hierarchy
  
  tmp0 = rownames(input_table3)
  tmp0 = gsub("bin_R.", "", tmp0)
  tmp0 = gsub("R", "", tmp0)
  tmp0 = gsub("synth", "", tmp0)
  tmp = sapply(tmp0, function(x)  strsplit(x, split=".",fixed=T)[[1]])
  tmp2 = max(unlist(lapply(tmp, length)))
  tmp3=c()
  for(i in 1:tmp2){
    tmp3 = cbind(tmp3, unlist(lapply(tmp, function(y) y[i])))
  }
  tmp3[is.na(tmp3)] = 0
  tmp4 = matrix(as.numeric(tmp3), ncol=ncol(tmp3))
  rownames(tmp4) = rownames(input_table3)
  for(i in ncol(tmp4):1){
    tmp4 = tmp4[order(tmp4[,i]),]
  }
  #tmp4 = tmp4[c(nrow(tmp4), 1:(nrow(tmp4)-1)),]
  
  
  input_table3 = input_table3[rownames(tmp4),]
  
  
  
  
  # recalculate the weights
  tmp = sapply(rownames(input_table3), function(x) sum( grepl(x,rownames(input_table3), fixed=T) & grepl("synth",rownames(input_table3), fixed=T))) 	
  input_table3[["weight_approx"]] = tmp
  input_table3[["weight_approx"]][1] = sum(grepl("synth",rownames(input_table3), fixed=T))
  input_table3 = input_table3[as.vector(input_table3[["weight_approx"]]) > 0,]
  # write the result
  
  write.table(input_table3, 'vtm.txt', sep=";", row.names=F)
  
  # run the VTM
  system(paste("java -jar Voronoi-Treemap-Library/build/libs/JVoroTreemap.jar", paste(getwd(), '/vtm.txt', sep=""), sep=" "))
  
  
  # retrive the result 
  vtm = read.table('vtm-finished.txt', header=T, sep=";")
  
  #remove the file
  system(paste("rm", paste(getwd(), "/vtm-finished.txt", sep="")))
  
  # fill the names and bin info
  vtm = vtm[is.na(vtm[["nodeId"]]) == F & vtm[["nodeId"]] != 0,]
  vtm[["nodeId"]] = as.vector(input_table3[["nodeId"]])[as.vector(vtm[["nodeId"]])]
  vtm[["parentID"]] = as.vector(input_table3[as.vector(vtm[["nodeId"]]), "parentId"])
  vtm[["name"]] = as.vector(input_table3[as.vector(vtm[["nodeId"]]), "name"])
  vtm[["weight"]] = as.vector(input_table3[as.vector(vtm[["nodeId"]]), "weight"])
  
  # add the ratio between the number of protein and the sum of intensities
  tmp = unlist(lapply(O2, function(x) sum(x %in% Ynames)))
  tmp2 = unlist(lapply(O2, length))
  tmp3 = tmp[vtm[["nodeId"]]] / tmp2[vtm[["nodeId"]]] 
  tmp3[is.na(tmp3)] = 0
  vtm[["intensity"]] = tmp3
  
  return(vtm)
}


### clean

rm(list=setdiff(ls(), c("Y", "O2","O1", "bins.description2", "f.vtm", "f.vtm.nonweight")))


### run for all


vtm = f.vtm.nonweight(rownames(Y),O2,bins.description2)

## reorder the vtm
rownames(vtm) = vtm[["nodeId"]]
tmp0 = rownames(vtm)
tmp0 = gsub("bin_", "", tmp0)
tmp0 = gsub("R", "", tmp0)
tmp0 = gsub("synth", "", tmp0)
tmp = sapply(tmp0, function(x)  strsplit(x, split=".",fixed=T)[[1]])
tmp2 = max(unlist(lapply(tmp, length)))
tmp3=c()
for(i in 1:tmp2){
  tmp3 = cbind(tmp3, unlist(lapply(tmp, function(y) y[i])))
}
tmp3[is.na(tmp3)] = 0
tmp4 = matrix(as.numeric(tmp3), ncol=ncol(tmp3))
rownames(tmp4) = rownames(vtm)
for(i in ncol(tmp4):1){
  tmp4 = tmp4[order(tmp4[,i]),]
}
#tmp4 = tmp4[c(nrow(tmp4), 1:(nrow(tmp4)-1)),]

vtm = vtm[rownames(tmp4),]

## remove synthetic nodes
vtm2 = vtm[grepl("synth", vtm[["nodeId"]]) == F,]
rownames(vtm2) = vtm2[["nodeId"]]

#identify the cell names
tmp = apply(vtm2, 1, function(x) strsplit(as.vector(t(x))[3],split=".", fixed=T)[[1]][length(strsplit(as.vector(t(x))[3],split=".", fixed=T)[[1]])]   )
vtm2[['shortname']] = tmp

#detect centroids for the labels
tmp =  apply(vtm2, 1, function(x) paste(colMeans(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T)),collapse=","))
vtm2[["center"]] = tmp

# associate the bottom with the top layer

tmp_top = as.vector(vtm2[["nodeId"]])[as.vector(vtm2[["hierarchyLevel"]]) == 2]
tmp_lowtotop = unlist(sapply(vtm2[["nodeId"]], function(x) names(which(unlist(sapply(tmp_top, function(y) grepl(y,x, fixed=T)))))[1]))

# sizes of the labels
tmp = unlist(apply(vtm2,1,function(x) max(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T)[,1]) - min(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T)[,1])))
tmp = (tmp/max(tmp))*3
tmp[tmp < 0.5] = 0.5
vtm2[["weightcex"]] = tmp



#  add average expression data to the polygons/terms

tmp = matrix(0,nrow = nrow(vtm2), ncol = 4)
for(i in 1:ncol(Y)){
  tmp[,i] = unlist(sapply(rownames(vtm2), function(x) mean(Y[O2[[x]],i],na.rm = T )))
}
colnames(tmp) = 1:4


tmp[tmp < -2] = -2
tmp[tmp >  2] = 2


vtm2[["comparisons"]] = tmp


# color the polygons

tmp2 = tmp
tmp2[] = "white"
tmp2[tmp>0] = "firebrick3"
tmp2[tmp<0] = "dodgerblue3"

tmp3 = tmp2
for( i in 1:nrow(tmp)){
  for(j in 1:ncol(tmp)){
    tmp3[i,j] = adjustcolor(tmp2[i,j], alpha = abs(tmp[i,j])/2)
  }
}

vtm2[["colors"]] = tmp3

# significance

tmp4 = tmp
tmp4[] = 1

for( i in 1:nrow(tmp)){
  for(j in 1:ncol(tmp)){
    tmp4[i,j] = wilcox.test(Y[O2[[rownames(vtm2)[i]]],j],Y[O2[[rownames(vtm2)[i]]],])$p.value
  }
}

tmp4[is.na(tmp4)] = 1

for(i in 1:4){
  tmp4[,i] = fdrtool(tmp4[,i], statistic = "pvalue")$qval
}

vtm2[["significance"]] = tmp4




# assign colors to the top layer
tmp_top = as.vector(vtm2[["nodeId"]])[as.vector(vtm2[["hierarchyLevel"]]) == 2]
vtm2_top = vtm2[tmp_top,]

# identify the bottom layer
tmp_last = which(vtm2[["nodeId"]] %in% vtm2[["parentID"]] == F)
vtm2_last = vtm2[tmp_last,]


##################################
## Plot

pdf(paste("Figure X2 (Voronoi comparison 1).pdf",sep=""), 10,10, pointsize=12)
par(xpd=T)
root = matrix(c(0,1,1,0,0,0,1,1),ncol=2)*1000
plot(root, type="n", axes=F, xlab="", ylab="", main="Time specificity")
apply(vtm2_last, 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T), col=x['colors.1'],border="light grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 3,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=1,border="dark grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 2,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=2,border="gray50"))
coverU = vtm2[,"comparisons"][,1] > 0 & vtm2[,"significance"][,1] < 0.05
coverD = vtm2[,"comparisons"][,1] < 0 & vtm2[,"significance"][,1] < 0.05
text(vtm2[coverU,6:7], vtm2[coverU,"shortname"], cex=vtm2[coverU,"weightcex"], col="firebrick4")
text(vtm2[coverD,6:7], vtm2[coverD,"shortname"], cex=vtm2[coverD,"weightcex"], col="dodgerblue4")
cover = vtm2_top[,"significance"][,1] >= 0.05
text(vtm2_top[cover,6:7], vtm2_top[cover,"shortname"], cex=vtm2_top[cover,"weightcex"], col="dark grey")
dev.off()

pdf(paste("Figure X2 (Voronoi comparison 2).pdf",sep=""), 10,10, pointsize=12)
par(xpd=T)
root = matrix(c(0,1,1,0,0,0,1,1),ncol=2)*1000
plot(root, type="n", axes=F, xlab="", ylab="", main="Time specificity")
apply(vtm2_last, 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T), col=x['colors.2'],border="light grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 3,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=1,border="dark grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 2,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=2,border="gray50"))
coverU = vtm2[,"comparisons"][,2] > 0 & vtm2[,"significance"][,2] < 0.05
coverD = vtm2[,"comparisons"][,2] < 0 & vtm2[,"significance"][,2] < 0.05
text(vtm2[coverU,6:7], vtm2[coverU,"shortname"], cex=vtm2[coverU,"weightcex"], col="firebrick4")
text(vtm2[coverD,6:7], vtm2[coverD,"shortname"], cex=vtm2[coverD,"weightcex"], col="dodgerblue4")
cover = vtm2_top[,"significance"][,2] >= 0.05
text(vtm2_top[cover,6:7], vtm2_top[cover,"shortname"], cex=vtm2_top[cover,"weightcex"], col="dark grey")
dev.off()

pdf(paste("Figure X2 (Voronoi comparison 3).pdf",sep=""), 10,10, pointsize=12)
par(xpd=T)
root = matrix(c(0,1,1,0,0,0,1,1),ncol=2)*1000
plot(root, type="n", axes=F, xlab="", ylab="", main="Time specificity")
apply(vtm2_last, 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T), col=x['colors.3'],border="light grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 3,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=1,border="dark grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 2,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=2,border="gray50"))
coverU = vtm2[,"comparisons"][,3] > 0 & vtm2[,"significance"][,3] < 0.05
coverD = vtm2[,"comparisons"][,3] < 0 & vtm2[,"significance"][,3] < 0.05
text(vtm2[coverU,6:7], vtm2[coverU,"shortname"], cex=vtm2[coverU,"weightcex"], col="firebrick4")
text(vtm2[coverD,6:7], vtm2[coverD,"shortname"], cex=vtm2[coverD,"weightcex"], col="dodgerblue4")
cover = vtm2_top[,"significance"][,3] >= 0.05
text(vtm2_top[cover,6:7], vtm2_top[cover,"shortname"], cex=vtm2_top[cover,"weightcex"], col="dark grey")
dev.off()

pdf(paste("Figure X2 (Voronoi comparison 4).pdf",sep=""), 10,10, pointsize=12)
par(xpd=T)
root = matrix(c(0,1,1,0,0,0,1,1),ncol=2)*1000
plot(root, type="n", axes=F, xlab="", ylab="", main="Time specificity")
apply(vtm2_last, 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T), col=x['colors.4'],border="light grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 3,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=1,border="dark grey" ))
apply(vtm2[as.vector(vtm2[["hierarchyLevel"]]) == 2,], 1, function(x) polygon(matrix(as.numeric(strsplit(as.vector(t(x))[9],split=",")[[1]]),ncol=2, byrow=T),lwd=2,border="gray50"))
coverU = vtm2[,"comparisons"][,4] > 0 & vtm2[,"significance"][,4] < 0.05
coverD = vtm2[,"comparisons"][,4] < 0 & vtm2[,"significance"][,4] < 0.05
text(vtm2[coverU,6:7], vtm2[coverU,"shortname"], cex=vtm2[coverU,"weightcex"], col="firebrick4")
text(vtm2[coverD,6:7], vtm2[coverD,"shortname"], cex=vtm2[coverD,"weightcex"], col="dodgerblue4")
cover = vtm2_top[,"significance"][,4] >= 0.05
text(vtm2_top[cover,6:7], vtm2_top[cover,"shortname"], cex=vtm2_top[cover,"weightcex"], col="dark grey")
dev.off()

