#/usr/bin/env Rscript
# This is a Rscript for KEGG analys 
# Liu, Shaorui; liushaorui@mail.bnu.edu.cn; 2018-1-3

# transmit args from Rscript
ORG <- commandArgs(T)

options(bitmapType='cairo')

ID_Convert <- read.table(paste0("../ID_Convert_",ORG))

# start analysis
library(clusterProfiler)

# read DEGs file from KEGG.sh output
a <- read.table("DEGs", colClasses = "character",comment.char="")
names(a)=c("SYMBOL","ENTREZID","log2FC_colour")

# select ENGTREZID for KEGG enrichment analysis
b <- a[,2]

# enrichKEGG, org need to be an outside value
# maximize Cutoff and GSSize for maximum pathway output
kkk1 <- enrichKEGG(gene=b, organism=ORG[1],pvalueCutoff=1,
qvalueCutoff=1,minGSSize=1,maxGSSize=length(b))

# set for summary.xls output
kkk2 <- as.data.frame(kkk1)
# a xls file starting with "ID" usually causes error, change it to "PathwayID"
names(kkk2)[1] <- "PathwayID"

# convert geneID to geneSymbol
geneSymbol <- kkk2$geneID

for(i in 1:length(geneSymbol)){
  G <- geneSymbol[i]
  G1 <- strsplit(G, split="/")[[1]]
  G2 <- ID_Convert[ID_Convert[,2] %in% G1 ,1]
  G3 <- paste(G2, collapse = "/")
  geneSymbol[i] <- G3
}

kkk2 <- cbind(kkk2,geneSymbol)
write.table(kkk2, file="KEGGSummary.xls",sep='\t',row.names=F,quote=F)

# ggplot2, dotplot the result of KEGG pathway enrichment 
library(ggplot2)
# make objects in "GeneRatio"&"BgRatio" numerable, to get Rich Factor
kkk3 <- kkk2
for (i in 1:length(kkk3[,3]))
kkk3[i,3] <- unlist(strsplit(kkk3[i,3],split="/"))[1]
for (i in 1:length(kkk3[,4]))
kkk3[i,4] <- unlist(strsplit(kkk3[i,4],split="/"))[1]
# make objects in "Description", nchar[i] <= 39, to keep the width of pdf/png 
for (i in 1:length(kkk3[,2]))
(if(nchar(kkk3[i,2])<=39){kkk3[i,2] <- kkk3[i,2]}
else{kkk3[i,2] <- paste0(substr(kkk3[i,2],1,36),"...")})

# make Description factor, levels stay with the order
kkk3[,2] <- factor(kkk3[,2], levels = rev(as.vector(kkk3[,2])))

# start ggplot2, choose top 20(default by p.adjust), set x,y-axis
p <- ggplot(head(kkk3,n=20),
aes(as.numeric(GeneRatio)/as.numeric(BgRatio),Description))
# set size&color
p1 <- p + geom_point(aes(size=Count,color=-1*log10(p.adjust)))+scale_size(range=c(4,8))
# plus with text 
p2 <- p1 +geom_text(aes(label=Count),size=3)
# set gradient color
p3 <- p2 + scale_color_gradient(low = "#FFFF94",high = "red")
# set coordinate details 
p4 <- p3 + labs(color=expression(-log[10](FDR)),
size="Gene Count",x="Rich Factor",y="Pathway Description",
title ="Top 20 KEGG Pathway Enrichment")
# cancle grid lines
p5 <- p4 + theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# save as ...
ggsave("KEGG_Pathway_Enrichment.pdf",width=6.57,height=5.5,plot=p5)
p6 <- p5 + ylab("Pathway Description") + theme(axis.title.y=element_text(hjust=0.88))
ggsave("KEGG_Pathway_Enrichment.png",width=6.44,height=5.5,plot=p6)

# define a function of kegg URL convertion for each line of data kkk1$geneID 
keggURL <- function(x) {
c <- a[,2:3]
d1 <- unlist(strsplit (x,split="/"))
d2 <- data.frame(d1)
names(c) <- c("ENTREZID","Colour")
names(d2) <- "ENTREZID"
d2$Num <- 1:lengths(d2)
d3 <- merge(d2,c,by="ENTREZID",all.x=TRUE)
d4 <- d3[order(d3$Num),,drop=FALSE]
row.names(d4) <- row.names(d3)
d5 <- paste0(d4[,1],"%09%23",d4[,3])
d6 <- paste(d5,collapse="/")
return(d6)
} 
# splice matched KEGG link, write into .xls, download url to .arg file  
d <- kkk1$geneID
for (i in 1:length(d)) d[i] <- keggURL(d[i])
e1 <- kkk1$ID
e2 <- d
e3 <- gsub("/","/mmu:",e2)
e3 <- gsub("^","/mmu:",e3)
e4 <- paste("http://www.genome.jp/kegg-bin/show_pathway?",e1,e3,sep="")
e5 <- cbind(e1,e4)
e6 <- as.data.frame(e5)
names(e6) <- c("PathwayID","PathwayAddress")
write.table(e6,file="KEGG_Pathway_Link.xls",sep="\t",row.names=F,quote=F)
#for (i in 1:length(e5[,1])) download.file(url=e5[i,2],
#destfile=paste(e5[i,1],".args",sep=""))



