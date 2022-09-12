#!/usr/bin/env bash

BED_ratio=$1
count_TSS=$2
count_TES=$3

#function to count reads at TSS/TES
countREADS () {
cat $BED_ratio | grep $1 | awk '{print $4}' | awk -F '_' '{print $1}' | sort | uniq -c | sed 's/^ *//' > $2
}

#create file for TSS
countREADS '_TSS_' $count_TSS

#create file for TES
countREADS '_TES_' $count_TES

