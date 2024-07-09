#!/usr/bin/env Rscript
# This is a Rscript for getting KEGG Organism code list
# Liu, Shaorui; liushaorui@mail.bnu.edu.cn; 2018-1-24

library(KEGGREST)
org <- keggList("organism")
write.table(org,"keggList.xls",sep="\t",row.names=F,quote=F)
