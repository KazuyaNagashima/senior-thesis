setwd('/home/ohganelab/Desktop/nagashima/comparison/19.10.31/ngsplot/k-means/k-means/R')

# golist = gene name list
head(go.list[[1]][[1]])

#length of genes
length(enrichList[[1]][,1])
length(enrichList[[2]][,1])
length(enrichList[[3]][,1])

#number of cluster
unique(unlist(go.list[[2]]))

#extract genes
genes <- unlist(go.list[[1]])
genes <- gsub(':.*$','', genes)
names(genes) <- NULL
cluster <- as.numeric(unlist(go.list[[2]]))
ensIds <- unlist(go.list[[1]])
ensIds <- gsub('^.*:','',ensIds)
names(ensIds) <- NULL
ensIds

#confirm data
head(genes[cluster == 1])
head(genes[cluster == 2])

#comfirm length of cluster5
enrichList[[1]][cluster == 5,]
length(genes[cluster == 5])
View(enrichList[[1]][cluster == 5,])

#comfirm length of cluster
for(i in 1:10) {
  cat("count of cluster",i,":",length(genes[cluster == i]), '\n')
}


#extract data
for (i in 1:10) {
  write.csv(genes[cluster == i],file = paste0("genename_cluster",i,".txt"),row.names = F, quote = F)
}
####最初の行にxがつくので消す####
