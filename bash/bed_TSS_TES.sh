#!/bin/bash

#based on https://bioinformatics.stackexchange.com/questions/19493/count-reads-at-specific-gene-features

GTF=$1
GTF_OUT=$2

awk -F"\t|\"" -vOFS="\t" '$3=="gene" {
  TSS=($7=="+" ? $4-2000"\t"$4+2000 : $5-2000"\t"$5+2000)
  TES=($7=="+" ? $5-2000"\t"$5+2000 : $4-2000"\t"$4+2000)
  print $1,TSS,$10"_TSS_2kb",0,$7
  print $1,TES,$10"_TES_2kb",0,$7
}' $GTF > $GTF_OUT
