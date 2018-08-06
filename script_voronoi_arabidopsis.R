require(deldir)
require(nloptr)
require(DEoptim)
require(gpclib)
require(sp)
require(tripack)
require(geometry)



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
mapman[1,1] = "bin_R"
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

bins = c(unique(unlist(bins)), bins)

bins.description[1] = "root"
names(bins.description)[1] = "bin_R."



### clean
rm(mapman)
rm(list=ls(pattern="tmp"))

### add hand edited short names
bins.description2 = bins.description
load('shortnames')
names(shortnames) = gsub("bin_", "bin_R.", names(shortnames), fixed = T)
names(shortnames) = gsub("bin_R.R.", "bin_R.", names(shortnames), fixed = T)
shortnames = shortnames[names(bins.description2)]

tmp = cbind(names(bins.description2), bins.description2, shortnames)
write.table(tmp, "shortnames_edit.txt", sep = "\t")
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




