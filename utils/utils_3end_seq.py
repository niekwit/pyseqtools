#!/usr/bin/env python3

import glob
import os
import sys
import subprocess
import shutil



script_dir = os.path.abspath(os.path.dirname(__file__))
script_dir = os.path.dirname(script_dir)
sys.path.append(os.path.join(script_dir, "utils"))
import utils_general as utils

'''
Reads were trimmed using Cutadapt v1.9.1 (Martin, 2011) to remove 5bp 5′ anchored barcodes and then aligned against the Homo sapiens GRCh38 genome build using Bowtie2 v2.2.9 (Langmead and Salzberg, 2012). Resulting genome alignment BAM files were merged, sorted and indexed using Picard v2.1.1

A set of high-confidence cleavage sites were computed from 3′-Seq data by merging reads from all 3′-Seq samples (regardless of KO status). Reads not overlapping with any annotated Ensembl exon in a strand-specific manner or any with a mapping quality < 10 were discarded. Remaining reads mapping to forward strand exons were filtered to include those with 3′ soft-clipped regions containing an “AA” at their 3′ terminus or “AAA” in any portion of the soft-clipped region. “TT” and “TTT” were used as equivalent soft-clipped region filters at the 5′ end of reads mapping to reverse strand exons. The filtered soft-clipped reads were used to create a read-depth coverage tracks and isolate islands that were at least 10 bp deep in a strand-specific manner. This resulted in 40,215 cleavage sites. A high-confidence cleavage set was generated from cleavage sites associated with a single gene (39,256), that map to either the terminal exon (32,235) or 3′UTR (26,883) of that gene, that have a merged read count of > 300 (28,106), contain at least 1 canonical or non-canonical polyA site (33,866) and represent at least 5% of all merged-reads belonging to the total exonic portions of the gene (20,723). This resulted in 17,835 high confidence cleavage sites. The relative expression difference (RED) scores for genes harboring ≥ 2 high-confidence cleavage sites were calculated using the 3′-Seq read count in the ± 50bp interval surrounding the midpoint of each cleavage site. Cleavage sites mapping to more than one gene were not considered. The usage of cleavage sites was normalized to the usage of the 3′most cleavage site (longest isoform) and taken relative to the WT cells using a previously described method for assessing relative expression differences (RED) from 3′Seq data (Li et al., 2015). RED = log2(short#/long#)KO – log2(short#/long#)WT. Thus, a positive RED score indicates an enrichment of the short isoform in the KO relative to the WT.
'''


def trim(work_dir):
    '''
    Use cutadapt to remove 5bp 5′ anchored barcodes
    '''
    pass

def align(work_dir):
    #cutadapt trim 5bp 5' end
    #align to hg38 bowtie2
    #merge BAM files, sort, and index
    
    
    pass