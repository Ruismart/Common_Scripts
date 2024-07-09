#!/bin/bash
# This is a Script for getting KEGG Organism code list
# Liu, Shaorui; liushaorui@mail.bnu.edu.cn; 2018-1-24

#Rscript /home/shaorui/Install/KEGG/keggList.R
Rscript $(which keggList.R)
cut -f 2,3 keggList.xls > keggList
rm keggList.xls 
