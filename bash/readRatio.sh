#!/usr/bin/env bash

BED_ratio=$1
count_TSS=$2
count_TES=$3

#function to count reads at TSS/TES
countREADS () {
cat $BED_ratio | grep _TSS_ | awk '{print $4}' | awk -F '_' '{print $1}' | sort | uniq -c | sed 's/^ *//' > $1
}

#create file for TSS
countREADS $count_TSS

#create file for TES
countREADS $count_TES

