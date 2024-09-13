##
# this is a custom one line script to run GO Enrichment
#     usage:  Rscript runGO.new.r Genelist.txt Outputname DBname Org
#         
#         Genelist       a txt file with DEGs for enrichment     
#         Outputname     output would be 'Outputname.Org.GODBname.csv/pdf'
#         DBname         BP/MF/CC
#         Org            only Hs(human)/Mm(mouse) available
#
# Shaorui, Liu; liushaorui@mail.bnu.edu.cn
## 

args <- commandArgs(T)

if(length(args)!=4){
  
  stop(paste0("This is a custom one line script to run GO Enrichment!\n",
  "           by Shaorui, Liu; liushaorui@mail.bnu.edu.cn.\n",
  "usage: Rscript runGO.new.r Genelist.txt Outputname DBname Org\n",
  "    Genelist     a txt file with DEGs for enrichment\n",
  "    Outputname   output would be 'Outputname.Org.GODBname.csv/pdf'\n",
  "    DBname       BP/MF/CC\n",
  "    Org          only Hs(human)/Mm(mouse) available\n"))
  
}else{

library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Mm.eg.db)
library(org.Hs.eg.db)
library(ggplot2)
options(bitmapType="cairo")

input_list <- args[1]
sample_name <- args[2]
db_name <- args[3]
org_name <- args[4]

#if(is.null(db_name)){
#  db_name <- "BP"
#}
#if(is.null(org_name)){
#  org_name <- "Mm"
#}


genes <- as.vector(unlist(read.table(input_list)))

ORG <- paste0("org.",org_name,".eg.db")

ids <- mapIds(x = get(ORG),
              keys = genes,
              keytype = "SYMBOL",
              column = "ENTREZID")
ids <- na.omit(ids)
#head(ids)

enrich.go <- enrichGO(gene = ids,
                         OrgDb = get(ORG),
                         keyType = "ENTREZID",
                         ont = db_name,
                         pvalueCutoff = 0.99,
                         qvalueCutoff = 0.99)
#dotplot(enrich.go_bp, showCategory=20)
#head(enrich.go_bp@result)

result <- enrich.go@result
geneSymbol <- result$geneID
ID_Convert <- data.frame(row.names = ids, symbol= names(ids))

for(i in 1:length(geneSymbol)){
  G <- geneSymbol[i]
  G <- strsplit(G, split="/")[[1]]
  G <- ID_Convert[G,]
  G <- paste(G, collapse = "/")
  geneSymbol[i] <- G
}
result <- cbind(result,geneSymbol)
colnames(result)[1] <- "PathwayID"

# output1
write.table(result, file = paste0(sample_name,".",org_name,".GO",db_name,"_Summary.csv"), sep = ",", row.names = F, quote = F)


# make objects in "GeneRatio"&"BgRatio" numerable, to get Rich Factor
kkk3 <- result[1:50,]

for (i in 1:length(kkk3[,3]))
  kkk3[i,3] <- unlist(strsplit(kkk3[i,3],split="/"))[1]
for (i in 1:length(kkk3[,4]))
  kkk3[i,4] <- unlist(strsplit(kkk3[i,4],split="/"))[1]
# make objects in "Description", nchar[i] <= 39, to keep the width of pdf/png
for (i in 1:length(kkk3[,2]))
  (if(nchar(kkk3[i,2])<=39){kkk3[i,2] <- kkk3[i,2]}
   else{kkk3[i,2] <- paste0(substr(kkk3[i,2],1,36),"...")})

kkk3 <- kkk3[!duplicated(kkk3[,2]),]  # avoid duplicated levels !

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
                title =paste0("Top 20 GO:",db_name," Pathway Enrichment"))
# cancle grid lines
p5 <- p4 + theme_bw() + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
# save as ...
ggplot2::ggsave(paste0(sample_name,".",org_name,".GO",db_name,"_Pathway_Enrichment.pdf"),width=6.57,height=5.5,plot=p5,device = "pdf")
#p6 <- p5 + ylab("Pathway Description") + theme(axis.title.y=element_text(hjust=0.88))
ggplot2::ggsave(paste0(sample_name,".",org_name,".GO",db_name,"_Pathway_Enrichment.png"),width=6.57,height=5.5,plot=p5,device = "png")

}








