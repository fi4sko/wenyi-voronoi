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


