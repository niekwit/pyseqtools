#!/usr/bin/env python3

"""
Converts MACS3 peak file to bed format for peak annotation with HOMER
"""

import pandas as pd
import sys


#get parsed arguments
args = sys.argv

macs3_file = args[1]
bed_file = args[2]

#convert macs3 file to bed format
df_bed = pd.read_csv(macs3_file,
                     sep = "\t",
                     header = None)
df_bed = df_bed.drop(range(4,9), axis = 1)

#further annotate BED file (required by HOMER for peak annotation)
df_bed = pd.read_csv(bed_file, sep = "\t", header = None)
column_number = len(df_bed.columns)
if column_number == 4:
    
    df_bed["empty"] = ""
    df_bed["strand"] = "+"
    df_bed.to_csv(bed_file, sep = "\t", index = False, header = False)
elif column_number > 6:
    print(f"ERROR: too many columns found in {bed_file}")






