#!/bin/bash
# This is a Script for getting KEGG ID_Convert_$ORG
# Liu, Shaorui; liushaorui@mail.bnu.edu.cn; 2018-1-24

# transmit arguments from shell -O arg
getopts "O:" arg
 
case $arg in
  O)
   ORG=$OPTARG
    ;;
  ?)
    echo "unknow Org code, please set -O xxx"
    exit 1
    ;; 
esac

echo "The Organism code is: $ORG"

# 
wget http://rest.kegg.jp/get/br:${ORG}00001 -O ${ORG}00001.keg
export ORG
perl -alne '{if(/^C/){/PATH:$ENV{'ORG'}(\d+)/;$kegg=$1}
             else{print "$kegg\t$F[1]\t$F[2]" if /^D/ and $kegg;}
            }' ${ORG}00001.keg > ID_Convert_${ORG}_1

gawk -F "\t" '{print $3"\t"$2}' ID_Convert_${ORG}_1 |sort -k1 -u|grep ";" |sort -k1 |sed 's/;//g' > ID_Convert_${ORG}
rm ID_Convert_${ORG}_1 ${ORG}00001.keg




