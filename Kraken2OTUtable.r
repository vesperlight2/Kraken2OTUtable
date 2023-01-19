#A home made script for tranforming Kraken2 results to OTU table

library(reshape2)
library(taxonomizr)

k2_res <- read.csv("Kraken2_Results.csv",sep = "\t",header = F)[,1:3]

k2_res <- k2_res[k2_res[,1]=="C",]

k2_res_reform <- cbind(k2_res[,2],str_remove(k2_res[,2],"_[0-9]*$"),str_extract(k2_res[,2],"[0-9]*$"),k2_res[,3],as.numeric(rep(1,nrow(k2_res))))

colnames(k2_res_reform) <- c("seqid","sampleid","sampleno","taxid","Num")

k2_res_reform <-as.data.frame(k2_res_reform)

k2_res_sum <- reshape2::dcast(k2_res_reform,sampleid+taxid~Num,length)

k2_res_spread <- dcast(k2_res_sum,taxid~sampleid)

k2_res_spread[is.na(k2_res_spread)] <- 0

k2_res_spread <- k2_res_spread[apply(as.data.frame(k2_res_spread[,-1]), 1, sum)>2,]

#taxadmp need to be prepared

dmp_db <- "taxadmp/"

taxanode <- read.nodes(paste0(dmp_db,"nodes.dmp"))

taxaname <- read.names(paste0(dmp_db,"names.dmp"))

taxa <- getTaxonomy2(k2_res_spread[,1],taxanode,taxaname)

taxa[is.na(taxa)] <- ""

taxa_combine <- paste0("k__",taxa[,1],";p__",taxa[,2],";c__",taxa[,3],";o__",taxa[,4],";f__",taxa[,5],";g__",taxa[,6],";s__",taxa[,7])

otu_table <- cbind(k2_res_spread,taxonomy=taxa_combine)

otu_table[,1] <- paste0("OTU_",1:nrow(otu_table))

id_list <- cbind(paste0("OTU_",1:nrow(otu_table)),k2_res_spread[,1])

id_list1 <- merge(x=id_list,y=k2_res_reform,by.x="V2",by.y="taxid",all.x=T)

id_list2 <- id_list1[duplicated(id_list1[,"V2"])==F,]

colnames(otu_table)[1] <- "#OTU ID"

write.table(otu_table,"otu_table_total.txt"),row.names = F,quote = F,sep = "\t")