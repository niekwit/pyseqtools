#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import multiprocessing
import timeit
import time


def main():

    #create top-level parser
    parser = argparse.ArgumentParser(prog = 'pyseqtools.py',
                                     description = "Analysis pipelines for a variety of NGS data, and related tools")
    subparsers = parser.add_subparsers(help ='Available pipelines/tools:',
                                       dest ='module')

    #create subparsers for each module, i.e. functionality

    ####crispr screen subparser
    parser_crispr = subparsers.add_parser('crispr',
                                          description = "Analysis pipeline for CRISPR-Cas screens",
                                          help = 'CRISPR screen analysis')

    parser_crispr.add_argument("-l", "--library",
                               required = "--csv2fasta" not in sys.argv,
                               choices = crispr_library_list,
                               help = "CRISPR library")
    parser_crispr.add_argument("-t", "--threads",
                               required = False,
                               default = 1,
                               metavar="<int>",
                               type = str,
                               help = "Number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    parser_crispr.add_argument("-r", "--rename",
                               required = False,
                               action = 'store_true',
                               help = "Rename fastq files according to rename.config")
    parser_crispr.add_argument("-m", "--mismatch",
                               required = False,
                               choices = ("0","1"),
                               metavar = "N",
                               help = "Number of mismatches (0 or 1) allowed during alignment",
                               default = 0)
    parser_crispr.add_argument("-a", "--analysis",
                               required = False,
                               default = "mageck",
                               choices = ["mageck", "bagel2"],
                               help = "Statistical analysis with MAGeCK or BAGEL2 (default is MAGeCK)")
    parser_crispr.add_argument("-f", "--fdr",
                               required = False,
                               metavar = "<FDR value>",
                               default = 0.25,help="Set FDR cut off for MAGeCK hits (default is 0.25)")
    parser_crispr.add_argument("--cnv",
                               required=False,
                               metavar = "<CCLE cell line>",
                               default = None,help="Activate CNV correction for MAGeCK/BAGEL2 with given cell line")
    parser_crispr.add_argument("--go",
                               required = False,
                               action = 'store_true',
                               default = None,
                               help = "Gene set enrichment analysis with enrichR")
    parser_crispr.add_argument("--gene-sets",
                               required = False,
                               metavar = "<GO gene set>",
                               default = ["GO_Molecular_Function_2021",
                                          "GO_Cellular_Component_2021",
                                          "GO_Biological_Process_2021"],
                               help = "Gene sets used for GO analysis (default is GO_Molecular_Function_2021, GO_Cellular_Component_2021, and GO_Biological_Process_2021). Gene sets can be found on https://maayanlab.cloud/Enrichr/#stats")
    parser_crispr.add_argument("--essential-genes",
                               required = False,
                               metavar = "<Custom essential gene list>",
                               default = os.path.join(script_dir, "core-essential-genes.csv"), #"path to Hart list"
                               help = "Essential gene list (default is Hart et al 2015 Cell)")
    parser_crispr.add_argument("--csv2fasta",
                               required = False,
                               metavar = "<CSV FILE>",
                               default = None,
                               help = "Convert CSV file with sgRNA names and sequences to fasta format. The first column should contain sgRNA names and the second sgRNA sequences (headers can be anything).")
    parser_crispr.add_argument("--skip-fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Skip FastQC/MultiQC")
    parser_crispr.add_argument("--skip-stats",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Skip MAGeCK/BAGEL2")
    parser_crispr.add_argument("--slurm",
                             required = False,
                             action = 'store_true',
                             help = "Submit jobs to Cambridge HPC using SLURM")

    ####RNA-Seq subparser
    parser_rnaseq = subparsers.add_parser('rna-seq',
                                          description = "Analysis pipeline for RNA-Seq experiments",
                                          help='RNA-Seq analysis')

    parser_rnaseq.add_argument("--md5sum",
                             action = 'store_true',
                             required = False,
                             help = "Check md5sums of fastq files")
    parser_rnaseq.add_argument("-r", "--rename",
                             required = False,
                             action = 'store_true',
                             help = "Rename fq files")
    parser_rnaseq.add_argument("-t", "--threads",
                               required = False,
                               default = 1,
                               metavar = "<int>",
                               help = "Number of CPU threads to use (default is 1). Use max to apply all available CPU threads. For Salmon 8-12 threads are optimal")
    parser_rnaseq.add_argument("-g", "--genome",
                               required = False,
                               choices = rna_seq_genomeList,
                               help = "Reference genome")
    parser_rnaseq.add_argument("--trim",
                             action = 'store_true',
                             required = False,
                             help = "Quality trimming of data using Trim_galore!")
    parser_rnaseq.add_argument("--peTags",
                               required = False,
                               default = None,
                               help = "Comma-separated paired-end file tags (e.g. _R1_001.fq.gz,_R2_001.fq.gz). Only required when --trim argument is called.")
    parser_rnaseq.add_argument("-a", "--align",
                               required = False,
                               choices = ["salmon","star"],
                               help = "Program to align fastq files")
    parser_rnaseq.add_argument("--rsemIndex",
                             required = False,
                             nargs='+',
                             help = "Generate index for STAR alignment via RSEM (--rsemIndex genome read-length out-dir")
    parser_rnaseq.add_argument("--isoformAnalysis",
                             required = False,
                             default = None,
                             choices = [None,"miso","rmats"],
                             help = "Alternative isoform analysis with RSEM/MISO or rMATS. Default is None")
    parser_rnaseq.add_argument("-f", "--scaleFactors",
                             required = False,
                             action = 'store_true',
                             help = "Calculate scale factors from yeast spike-in RNA with DESeq2")
    parser_rnaseq.add_argument("--deseq2",
                             required = False,
                             action = 'store_true',
                             help = "Calculate differential genes using DESeq2")
    parser_rnaseq.add_argument("-b", "--bigwig",
                             required = False,
                             action = 'store_true',
                             help = "Create BigWig files using bamCoverage (requires indexed BAM files")

    parser_rnaseq.add_argument("--TE",
                               required = False,
                               action = 'store_true',
                               help = "Transposable element expression analysis (Requires STAR alignment)")
    parser_rnaseq.add_argument("--go",
                               required = False,
                               action = 'store_true',
                               help = "Gene set enrichment analysis with Enrichr")
    parser_rnaseq.add_argument("--gene-sets",
                               required = False,
                               metavar = "<GO gene set>",
                               default = ["GO_Molecular_Function_2021",
                                          "GO_Cellular_Component_2021",
                                          "GO_Biological_Process_2021"],
                               help = "Gene sets used for GO analysis (default is GO_Molecular_Function_2021, GO_Cellular_Component_2021, and GO_Biological_Process_2021). Gene sets can be found on https://maayanlab.cloud/Enrichr/#stats")
    parser_rnaseq.add_argument("-p", "--pvalue",
                               required = False,
                               metavar = "<P value>",
                               default = 0.001,
                               help = "Set P value cut off (default is 0.001)")
    parser_rnaseq.add_argument("--fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Run FastQC/MultiQC")
    parser_rnaseq.add_argument("--slurm",
                             required = False,
                             action = 'store_true',
                             help = "Submit jobs to Cambridge HPC using SLURM")

    ####ChIP-Seq subparser
    parser_chip = subparsers.add_parser('chip-seq',
                                        description = "Analysis pipeline for ChIP-Seq experiments",
                                        help='ChIP-Seq analysis')

    parser_chip.add_argument("--md5sum",
                                 action = 'store_true',
                                 required = False,
                                 help = "Check md5sums of fastq files")
    parser_chip.add_argument("-t", "--threads",
                             required = False,
                             default = "1",
                             type = str,
                             help = "<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    parser_chip.add_argument("-r", "--rename",
                             required = False,
                             action = 'store_true',
                             help = "Rename fq files")
    parser_chip.add_argument("-g", "--genome",
                             required = False,
                             default = "hg38",
                             help = "Choose reference genome (default is hg38)")
    parser_chip.add_argument("--trim",
                             action = 'store_true',
                             required = False,
                             help = "Quality trimming of data using Trim_galore!")
    parser_chip.add_argument("--peTags",
                               required = False,
                               default = None,
                               help = "Comma-separated paired-end file tags (e.g. _R1_001.fq.gz,_R2_001.fq.gz). Only required when --trim argument is called.")
    parser_chip.add_argument("--align",
                             choices = ["hisat2",
                                        "bwa-mem",
                                        "bwa-aln"],
                             const = "hisat2",
                             nargs = "?",
                             required = False,
                             help = "Create BAM files using HISAT2 or BWA (mem or aln)")
    parser_chip.add_argument("-d", "--deduplication",
                             required = False,
                             action = 'store_true',
                             help = "Perform deduplication of BAM files")
    parser_chip.add_argument("--downsample",
                             required = False,
                             action = 'store_true',
                             help = "Perform downsampling of BAM files to lowest read count")
    parser_chip.add_argument("-s", "--spike-in",
                             required = False,
                             action = 'store_true',
                             help = "Perform normalisation based on Drosophila chromatin spike-in")
    parser_chip.add_argument("-b", "--bigwig",
                             required = False,
                             default = False,
                             nargs = "?",
                             const = "bigwig",
                             help = "Create BigWig or BedGraph files")
    parser_chip.add_argument("--BAMqc",
                             required = False,
                             action = 'store_true',
                             help = "Perform QC analysis of BAM files")
    parser_chip.add_argument("-p", "--peaks",
                             required = False,
                             default = False,
                             choices = ["narrow","broad"],
                             metavar = "narrow | broad",
                             help = "Call and annotate peaks with MACS3/HOMER and calculate differential peaks with DiffBind")
    parser_chip.add_argument("--metagene",
                             required = False,
                             default = False,
                             nargs = "?",
                             const = None,
                             help = "Generate metageneplots and heatmaps with plotProfile (deepTools)\nNOTE: a text file can be parsed with gene names to create a plot for just these genes")
    parser_chip.add_argument("--slurm",
                             required = False,
                             action = 'store_true',
                             help = "Submit jobs to Cambridge HPC using SLURM")
    parser_chip.add_argument("--fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Run FastQC/MultiQC")


    ####cutrun-Seq subparser
    parser_cutrun = subparsers.add_parser('cutrun',
                                          description = "Analysis pipeline for CUT & RUN experiments",
                                          help = 'CUT & RUN analysis')
    parser_cutrun.add_argument("--md5sum",
                             action = 'store_true',
                             required = False,
                             help = "Check md5sums of fastq files")
    parser_cutrun.add_argument("-t", "--threads",
                             required = False,
                             default = 1,
                             type = str,
                             help = "<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    parser_cutrun.add_argument("-g", "--genome",
                             required = False,
                             default = 'hg19',
                             help = "Choose reference genome (default is hg19)")
    parser_cutrun.add_argument("-a", "--align",
                             choices = ["bowtie",
                                        "bowtie2",
                                        "hisat2",
                                        "bwa-mem",
                                        "bwa-aln"],
                             const = "bowtie",
                             nargs = "?",
                             required = False,
                             help = "Create BAM files using Bowtie2, HISAT2 or BWA (mem or aln")
    parser_cutrun.add_argument("-d", "--deduplication",
                             required = False,
                             action = 'store_true',
                             help = "Perform deduplication of BAM files")
    parser_cutrun.add_argument("-b", "--bigwig",
                             required = False,
                             action = 'store_true',
                             help = "Create BigWig files")
    parser_cutrun.add_argument("--qc",
                             required = False,
                             action = 'store_true',
                             help = "Perform QC analysis of BAM files")
    parser_cutrun.add_argument("-p", "--peaks",
                             required = False,
                             action = 'store_true',
                             help = "Call and annotate peaks with MACS3/HOMER")
    parser_cutrun.add_argument("--metagene",
                             required = False,
                             action = 'store_true',
                             help = "Generate metageneplots and heatmaps with plotProfile (deepTools)")
    parser_cutrun.add_argument("--skip-fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Skip FastQC/MultiQC")
    parser_cutrun.add_argument("--slurm",
                             required = False,
                             action = 'store_true',
                             help = "Submit jobs to Cambridge HPC using SLURM")

    ####DamID-Seq subparser
    parser_damid = subparsers.add_parser('damid',
                                          description = "Analysis pipeline for DamID (wrapper for https://github.com/owenjm/damidseq_pipeline)",
                                          help = 'DamID analysis')
    parser_damid.add_argument("-r", "--rename",
                               required = False,
                               action = 'store_true',
                               help = "Rename fastq files according to rename.config")
    parser_damid.add_argument("-t", "--threads",
                           required = False,
                           default = 1,
                           type = str,
                           help = "<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    parser_damid.add_argument("-g", "--genome",
                            required = False,
                            default = 'hg19',
                            help = "Choose reference genome (default is hg19)")
    parser_damid.add_argument("--slurm",
                             required = False,
                             action = 'store_true',
                             help = "Submit jobs to Cambridge HPC using SLURM")
    
    ####TT-Seq subparser
    parser_ttseq = subparsers.add_parser('tt-seq',
                                        description = "Analysis pipeline for TT-Seq experiments",
                                        help='TT-Seq analysis according to https://github.com/crickbabs/DRB_TT-seq')

    parser_ttseq.add_argument("--md5sum",
                                 action = 'store_true',
                                 required = False,
                                 help = "Check md5sums of fastq files")
    parser_ttseq.add_argument("-t", "--threads",
                             required = False,
                             default = "1",
                             type = str,
                             help = "<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    parser_ttseq.add_argument("-r", "--rename",
                             required = False,
                             action = 'store_true',
                             help = "Rename fastq files")
    parser_ttseq.add_argument("-a", "--align",
                             action = 'store_true',
                             required = False,
                             help = "Alignment of fastq files using STAR")
    parser_ttseq.add_argument("-g", "--genome",
                             required = False,
                             type = str,
                             const = "hg38",
                             nargs = "?",
                             default = "hg38",
                             help = "Choose reference genome (default is hg38)")
    parser_ttseq.add_argument("--splitBAM",
                             required = False,
                             action = 'store_true',
                             help = "Generate forward and reverse strand specific BAM files")
    parser_ttseq.add_argument("--scaleFactors",
                             required = False,
                             action = 'store_true',
                             help = "Calculate size factors from yeast spike-in RNA with DESeq2")
    parser_ttseq.add_argument("-b", "--bigwig",
                             required = False,
                             action = 'store_true',
                             help = "Create BigWig files")
    parser_ttseq.add_argument("--deseq2",
                             required = False,
                             action = 'store_true',
                             help = "Calculate differential genes using DESeq2")
    parser_ttseq.add_argument("--metagene",
                             required = False,
                             choices = ["deeptools","ngsplot", None],
                             default = None,
                             help = "Generate metagene plots and heatmaps with ngs.plot or deeptools")
    parser_ttseq.add_argument("--readRatio",
                             required = False,
                             action = 'store_true',
                             help = "Calculate ratios of read numbers in areas around TSS vs TES")
    parser_ttseq.add_argument("--slurm",
                             required = False,
                             action = 'store_true',
                             help = "Submit jobs to Cambridge HPC using SLURM")
    
    ####3end-Seq subparser
    parser_3endseq = subparsers.add_parser('3end-seq',
                                          description = "Analysis pipeline for 3end-Seq experiments",
                                          help='3end-Seq analysis')
    parser_3endseq.add_argument("--md5sum",
                             action = 'store_true',
                             required = False,
                             help = "Check md5sums of fastq files")
    parser_3endseq.add_argument("-t", "--threads",
                               required = False,
                               default = 1,
                               metavar = "<int>",
                               help = "Number of CPU threads to use (default is 1). Use max to apply all available CPU threads. For Salmon 8-12 threads are optimal")
    parser_3endseq.add_argument("-g", "--genome",
                               required = False,
                               choices = rna_seq_genomeList,
                               help = "Reference genome")
    parser_3endseq.add_argument("--trim",
                             action = 'store_true',
                             required = False,
                             help = "Quality trimming of data using Trim_galore!")
    parser_3endseq.add_argument("--peTags",
                               required = False,
                               default = None,
                               help = "Comma-separated paired-end file tags (e.g. _R1_001.fq.gz,_R2_001.fq.gz)")
    parser_3endseq.add_argument("-a", "--align",
                               required = False,
                               choices = ["salmon","star"],
                               help = "Program to align fastq files")
    parser_3endseq.add_argument("--indexBAM",
                             required = False,
                             action = 'store_true',
                             help = "Index BAM files with samtools")
    parser_3endseq.add_argument("--sortBAM",
                             required = False,
                             action = 'store_true',
                             help = "Sort BAM files with samtools")
    parser_3endseq.add_argument("-f", "--scaleFactors",
                             required = False,
                             action = 'store_true',
                             help = "Calculate scale factors from yeast spike-in RNA with DESeq2")
    parser_3endseq.add_argument("--deseq2",
                             required = False,
                             action = 'store_true',
                             help = "Calculate differential genes using DESeq2")
    parser_3endseq.add_argument("-b", "--bigwig",
                             required = False,
                             action = 'store_true',
                             help = "Create BigWig files using bamCoverage (requires indexed BAM files")

    parser_3endseq.add_argument("--TE",
                               required = False,
                               action = 'store_true',
                               help = "Transposable element expression analysis (Requires STAR alignment)")
    parser_3endseq.add_argument("--go",
                               required = False,
                               action = 'store_true',
                               help = "Gene set enrichment analysis with Enrichr")
    parser_3endseq.add_argument("--gene-sets",
                               required = False,
                               metavar = "<GO gene set>",
                               default = ["GO_Molecular_Function_2021",
                                          "GO_Cellular_Component_2021",
                                          "GO_Biological_Process_2021"],
                               help = "Gene sets used for GO analysis (default is GO_Molecular_Function_2021, GO_Cellular_Component_2021, and GO_Biological_Process_2021). Gene sets can be found on https://maayanlab.cloud/Enrichr/#stats")
    parser_3endseq.add_argument("-p", "--pvalue",
                               required = False,
                               metavar = "<P value>",
                               default = str(0.001),
                               help = "Set P value cut off (default is 0.001)")
    parser_3endseq.add_argument("--fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Run FastQC/MultiQC")
    parser_3endseq.add_argument("--slurm",
                             required = False,
                             action = 'store_true',
                             help = "Submit jobs to Cambridge HPC using SLURM")
    '''
    #create subparser for gene symbol conversion
    parser_conversion = subparsers.add_parser('genesymconv',
                                          description = 'Convert human gene symbols to mouse gene symbols, or vice versa',
                                          help = "Inter-species gene symbol conversion")
    parser_conversion.add_argument("-c", "--conversion",
                                   choices = ['hm',
                                        'mh'],
                             required = True,
                             help = "Conversion type: human to mouse (hm) or mouse to human (mh)")
    parser_conversion.add_argument("-i", "--input",
                             required = True,
                             help = "Input gene list file. Each gene symbol should be on a new line")
    parser_conversion.add_argument("-o", "--output",
                             required = True,
                             help = "Output file name")
    '''

    #create dictionary with command line arguments
    args = vars(parser.parse_args())


    def crispr(args, script_dir):

        #csv to fasta conversion
        csv = args["csv2fasta"]
        if csv is not None:
            if not os.path.isfile(csv):
                sys.exit("ERROR: invalid file path given")
            else:
                crispr_utils.csv2fasta(csv,script_dir)


        ###set thread count for processing
        threads = utils.set_threads(args)
        
        ##rename files
        rename = args["rename"]
        if rename == True:
            utils.rename(work_dir)
        
        ###Check md5 checksums
        utils.checkMd5(work_dir)
       
        #determine file extension raw data
        file_extension = utils.get_extension(work_dir)

        ##Run FastQC/MultiQC
        skip_fastqc = args["skip_fastqc"]
        if not skip_fastqc:
            utils.fastqc(script_dir, work_dir,threads,file_extension)
        else:
            print("Skipping FastQC/MultiQC analysis")

        ##count reads
        #check if bowtie2 index is build for CRISPR library
        crispr_library = args["library"]
        #crispr_utils.check_index(crispr_settings, crispr_library, script_dir, work_dir)

        #check if file with just guide names exists
        crispr_utils.guide_names(crispr_settings, crispr_library)

        #count sgRNAs
        mismatch = args["mismatch"]
        crispr_utils.count(crispr_settings,
                     crispr_library,
                     mismatch,
                     threads,
                     script_dir,
                     work_dir)
        
        #plot alignment rates
        crispr_utils.plot_alignment_rate(work_dir)

        #join count files
        if not utils.file_exists(os.path.join(work_dir,
                                    "count",
                                    'counts-aggregated.tsv')):
            crispr_utils.join_counts(work_dir, crispr_settings, crispr_library)
        #normalise read count table
        if not utils.file_exists(os.path.join(work_dir,
                                    "count",
                                    "counts-aggregated-normalised.csv")):
            crispr_utils.normalise(work_dir)

        #plot sample coverage (read count / library size)
        crispr_utils.plot_coverage(work_dir, crispr_settings, crispr_library)
        
        ##run library analysis
        crispr_utils.lib_analysis(work_dir, crispr_settings, crispr_library, script_dir)
        crispr_utils.gcBias(work_dir, crispr_settings, crispr_library)

        ##run stats on counts
        analysis = args["analysis"]
        go = args["go"]
        fdr = float(args["fdr"])
        cnv = args["cnv"]

        skip_stats = args["skip_stats"]
        if not skip_stats:
            if analysis == "mageck":
                crispr_utils.mageck(work_dir,script_dir,cnv,fdr)


            elif analysis == "bagel2":
                print("Running BAGEL2")
                crispr_utils.convert4bagel(work_dir, crispr_settings, crispr_library)
                crispr_utils.remove_duplicates(work_dir)
                crispr_utils.bagel2(work_dir, script_dir, fdr, crispr_settings, crispr_library)

        #run essential gene list comparison
        essential_genes = args["essential_genes"]
        crispr_utils.essentialGenes(work_dir, analysis, essential_genes, fdr)

        if go == True:
            gene_sets = args["gene_sets"]
            crispr_utils.goPython(work_dir,
                            fdr,
                            crispr_settings,
                            crispr_library,
                            analysis,
                            gene_sets)
        

    def rna_seq(args, script_dir):
        #get parsed arguments 
        module = args["module"]
        genome = args["genome"]
        align = args["align"]
        isoformAnalysis = args["isoformAnalysis"]
        threads = args["threads"]
        deseq2 = args["deseq2"]
        TE = args["TE"]
        trim = args["trim"]
        pvalue = args["pvalue"]
        deseq2 = args["deseq2"]
        bigwig = args["bigwig"]
        scaleFactors = args["scaleFactors"]
        pe_tags = args["peTags"]
        slurm = args["slurm"] 
        rsemIndex = args["rsemIndex"]
        md5sum = args["md5sum"]
        fastqc = args["fastqc"]
        rename = args["rename"]
        
        if slurm == False:
            ###set thread count for processing
            max_threads = str(multiprocessing.cpu_count())
            threads = args["threads"]
            if threads == "max":
                threads = max_threads

            ### Check md5sums
            #utils.checkMd5(work_dir,script_dir,slurm)
    
            ###Run FastQC/MultiQC
            file_extension = utils.get_extension(work_dir)
            fastqc = args["fastqc"]
            if fastqc == True:
                utils.fastqc(script_dir, work_dir, threads, file_extension)
                
            ###Set species variable
            if "hg" in genome or genome == "gencode-v35":
                species = "human"
            elif "mm" in genome or genome == "gencode.vM1.pc_transcripts":
                species = "mouse"


            ###trim and align
            pvalue = args["pvalue"]
            align = args["align"]
            if align != None:
                if align.lower() == "salmon":
                    gtf = rna_seq_settings["gtf"][genome]
                    utils.trim(script_dir, threads, work_dir, pe_tags)
                    salmon_index = rna_seq_settings["salmon_index"][genome]
                    gtf = rna_seq_settings["salmon_gtf"][genome]
                    fasta = rna_seq_settings["FASTA"][genome]
                    rnaseq_utils.salmon(salmon_index,
                                        str(threads),
                                        work_dir,
                                        gtf,
                                        fasta,
                                        script_dir,
                                        rna_seq_settings,
                                        genome)
                    rnaseq_utils.plotMappingRate(work_dir)
                    rnaseq_utils.plotPCA(work_dir, script_dir)
                    rnaseq_utils.diff_expr(work_dir, gtf, script_dir, species, pvalue, genome)
                    rnaseq_utils.plotVolcano(work_dir)
                elif align.lower() == "star":
                    rnaseq_utils.trim(script_dir, threads, work_dir, pe_tags)
                    rnaseq_utils.STAR(work_dir, 
                                      threads, 
                                      script_dir, 
                                      rna_seq_settings, 
                                      genome, 
                                      pe_tags, 
                                      slurm)
                
            
            go = args["go"]
    
            pvalue = args["pvalue"]
            
            if deseq2 == True:
                rnaseq_utils.diff_expr(work_dir,gtf,script_dir,species,pvalue,genome, slurm)
                
            if scaleFactors == True:
                tt_seq_utils.scaleFactors(script_dir, work_dir, slurm)
            
            if bigwig == True:
                rnaseq_utils.BigWig(work_dir, threads, genome, rna_seq_settings, slurm)   
            
            if TE == True:
                rnaseq_utils.retroElements(work_dir, script_dir, rna_seq_settings, threads, genome, slurm)
    
            if go == True:
                gene_sets=args["gene_sets"]
                rnaseq_utils.geneSetEnrichment(work_dir, pvalue, gene_sets)
        else:
                          
             #Set species variable
             if genome != None:
                 if "hg" in genome or genome == "gencode-v35":
                     species = "human"
                 elif "mm" in genome or genome == "gencode.vM1.pc_transcripts":
                     species = "mouse"
             
             if md5sum == True:
                 utils.checkMd5(work_dir,script_dir,slurm)   
             
             if rename == True:
                 utils.rename(work_dirls)
             
             if fastqc == True:
                 utils.fastqcSLURM(work_dir, script_dir)
             
             if rsemIndex != None:
                 if len(rsemIndex) == 3:
                     rnaseq_utils.rsemIndex(work_dir, script_dir, rna_seq_settings, slurm, rsemIndex)
             
             if trim == True:
                 utils.trimSLURM(script_dir, work_dir,module, pe_tags)
             
             if align == "star":
                 job_id_star = rnaseq_utils.slurmSTAR(work_dir,script_dir,genome)
                 
             if isoformAnalysis != None:
                 rnaseq_utils.isoformAnalysis(work_dir, script_dir, rna_seq_settings, genome, slurm, isoformAnalysis)
                
             if deseq2 == True:
                 gtf = rna_seq_settings["gtf"][genome.split("_")[0]]
                 rnaseq_utils.diff_expr(work_dir,gtf,script_dir,species,pvalue,genome, slurm)
                 
             if scaleFactors == True:
                 tt_seq_utils.scaleFactors(script_dir, work_dir, slurm)
             
             if bigwig == True:
                 rnaseq_utils.BigWig(work_dir, threads, genome, rna_seq_settings, slurm)   
             
             if TE == True:
                 job_id_star = rnaseq_utils.slurmSTAR(work_dir,script_dir,genome,TE=True)
                 rnaseq_utils.retroElementsSLURM(work_dir,script_dir,genome,job_id_star)
             

    def chip_seq(args, script_dir, module):
        slurm = args["slurm"]
        pe_tags = args["peTags"]
        peak = args["peaks"]
        md5sum = args["md5sum"]
        threads = args["threads"]
        align = args["align"]
        genome = args["genome"]
        trim = args["trim"]
        fastqc = args["fastqc"]
        bigwig = args["bigwig"]
        metagene = args["metagene"]
        dedup = args["deduplication"]
        downscale = args["downsample"]
        BAMqc = args["BAMqc"]
        
        #set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        if threads == "max":
            threads = max_threads

        #Check md5sums
        if md5sum == True:
            utils.checkMd5(work_dir,script_dir,slurm)

        #create BAM files
        
        
        #trimming for cluster
        
        if trim == True:
            if slurm == True:
                utils.trimSLURM(script_dir, work_dir, module, pe_tags)
        
        if slurm == False:
            if align is not None:
                if align == "hisat2":
                    ##Run FastQC/MultiQC
                    skip_fastqc = args["skip_fastqc"]
                    file_extension = utils.get_extension(work_dir)
                    if not skip_fastqc:
                        utils.fastqc(script_dir, work_dir,threads, file_extension)
                    else:
                        print("Skipping FastQC/MultiQC analysis")
    
                    utils.trim(script_dir, threads, work_dir)
                    chipseq_utils.hisat2(script_dir, work_dir, threads, chip_seq_settings, genome)
                    utils.indexBam(work_dir, threads)
                elif "bwa" in align:
                    ##Run FastQC/MultiQC
                    skip_fastqc = args["skip_fastqc"]
                    file_extension = utils.get_extension(work_dir)
                    if not skip_fastqc:
                        utils.fastqc(script_dir, work_dir,threads, file_extension)
                    else:
                        print("Skipping FastQC/MultiQC analysis")
    
                    utils.trim(script_dir, threads, work_dir)
                    utils.bwa(work_dir, script_dir, args, threads, chip_seq_settings, genome)
                    utils.indexBam(work_dir, threads)
            if peak == True:
                chipseq_utils.peak(work_dir, threads, genome, chip_seq_settings)
                
            if dedup == True:
                utils.deduplicationBam(script_dir, work_dir, threads, args)
                utils.indexBam(work_dir, threads)
            
            if BAMqc == True:
                chipseq_utils.bam_bwQC(work_dir, threads)
                
        else:
            if align == "hisat2":
                chipseq_utils.hisat2SLURM(script_dir, work_dir, threads, chip_seq_settings, genome)
                
            if fastqc == True:
                utils.fastqcSLURM(work_dir, script_dir)
                
            if peak != False:
                chipseq_utils.peakSLURM(work_dir, genome, peak)
                
            if bigwig != False:
                if bigwig == "bigwig":
                    job_id_bamcoverage = utils.bamCoverageSLURM(genome,"bigwig")
                    utils.bamCompareSLURM(genome)
                    utils.pcaBwSLURM(genome,job_id_bamcoverage,"bigwig")
                elif bigwig == "bedgraph":
                    job_id_bamcoverage = utils.bamCoverageSLURM(genome,"bedgraph")
            
            if metagene != False:
                chipseq_utils.plotProfileSLURM(genome,metagene)

            if dedup == True:
                utils.deduplicationSLURM(script_dir, work_dir, genome)

            if downscale == True:
                chipseq_utils.downsample(script_dir, work_dir, threads, genome, slurm)
                
            if BAMqc == True:
                    chipseq_utils.bamQCslurm(genome)
                    


    def cutrun(args, script_dir):
        slurm = args["slurm"]
        
        #set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        threads = args["threads"]
        if threads == "max":
            threads = max_threads

        #Check md5sums
        utils.checkMd5(work_dir,script_dir,slurm)


        #create BAM files
        align = args["align"]
        genome = args["genome"]

        if align is not None:
            if align == "hisat2":
                ##Run FastQC/MultiQC
                skip_fastqc = args["skip_fastqc"]
                file_extension = utils.get_extension(work_dir)
                if not skip_fastqc:
                    utils.fastqc(script_dir, work_dir,threads, file_extension)
                else:
                    print("Skipping FastQC/MultiQC analysis")

                utils.trim(script_dir, threads, work_dir)
                chipseq_utils.hisat2(script_dir, work_dir, threads, chip_seq_settings, genome)
                utils.indexBam(work_dir, threads)
            elif "bwa" in align:
                ##Run FastQC/MultiQC
                skip_fastqc = args["skip_fastqc"]
                file_extension = utils.get_extension(work_dir)
                if not skip_fastqc:
                    utils.fastqc(script_dir, work_dir,threads, file_extension)
                else:
                    print("Skipping FastQC/MultiQC analysis")

                utils.trim(script_dir, threads, work_dir)
                utils.bwa(work_dir, script_dir, args, threads, chip_seq_settings, genome)
                utils.indexBam(work_dir, threads)
            elif align == "bowtie":
                ##Run FastQC/MultiQC
                skip_fastqc = args["skip_fastqc"]
                file_extension = utils.get_extension(work_dir)
                if not skip_fastqc:
                    utils.fastqc(script_dir, work_dir,threads, file_extension)
                else:
                    print("Skipping FastQC/MultiQC analysis")

                utils.trim(script_dir, threads, work_dir)
                cutrun_utils.bowtie(work_dir, script_dir, str(threads), cutrun_settings, genome)
            elif align == "bowtie2":
                ##Run FastQC/MultiQC
                skip_fastqc = args["skip_fastqc"]
                file_extension = utils.get_extension(work_dir)
                if not skip_fastqc:
                    utils.fastqc(script_dir, work_dir,threads, file_extension)
                else:
                    print("Skipping FastQC/MultiQC analysis")

                utils.trim(script_dir, threads, work_dir)
                cutrun_utils.bowtie2(work_dir, script_dir, str(threads), cutrun_settings, genome)
                utils.indexBam(work_dir, threads)

        dedup = args["deduplication"]
        if dedup == True:
            utils.deduplicationBam(script_dir, work_dir, threads, args)
            utils.indexBam(work_dir, threads)

        bigwig = args["bigwig"]
        if bigwig == True:
            utils.createBigWig(work_dir, threads)

        qc = args["qc"]
        if qc == True:
            chipseq_utils.bam_bwQC(work_dir, threads)

        metagene = args["metagene"]
        if metagene == True:
            chipseq_utils.plotProfile(work_dir, chip_seq_settings, genome, threads)

        peak = args["peaks"]
        if peak == True:
            chipseq_utils.peak(work_dir, threads, genome, chip_seq_settings)
    
    
    def damID(args, script_dir):
        rename = args["rename"]
        threads = args["threads"]
        genome = args["genome"]
        slurm = genome = args["slurm"]
        
        #Check md5sums
        #utils.checkMd5(work_dir,script_dir,slurm)
        
        ##rename files
        if rename == True:
            utils.rename(work_dir)
        
        #set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        if threads == "max":
            threads = max_threads
        
        #utils.trim(script_dir, threads, work_dir)
        damid_utils.damID(script_dir, work_dir, threads, genome, damid_settings)
        damid_utils.bedgraph2BigWig(script_dir, work_dir, damid_settings, genome)
    
    
    def ttSeq(args, script_dir):
              
        #get parsed arguments
        align=args["align"]
        slurm = args["slurm"]
        genome = args["genome"]
        rename = args["rename"]
        scaleFactors = args["scaleFactors"]
        md5sum = args["md5sum"]
        splitBAM = args["splitBAM"]
        metagene = args["metagene"]
                
        #set thread count for processing
        if slurm == False:
            max_threads = str(multiprocessing.cpu_count())
            threads = args["threads"]
            if threads == "max":
                threads = max_threads
            print(f"Using {threads} CPU threads for analysis")
        
        #rename files
        if rename == True:
            utils.rename(work_dir)
        
        #check md5sums
        if md5sum == True:
            utils.checkMd5(work_dir,script_dir,slurm)
        
        #quality trim fastq files and align
        if align == True:
            if slurm == True:
                job_id_trim = utils.trimSLURM(script_dir, work_dir)
                job_id_align = tt_seq_utils.STAR(work_dir, threads, script_dir, tt_seq_settings, genome, slurm, job_id_trim)
                tt_seq_utils.bamSortSLURM(work_dir, job_id_align, genome="hg38")
            else:
                utils.trim(script_dir, threads, work_dir)
                tt_seq_utils.STAR(work_dir, threads, script_dir, tt_seq_settings, genome)
                
        #split bam files into forward and reverse strand files
        if splitBAM == True:
            if slurm == True:  
                tt_seq_utils.splitBam(work_dir, genome, slurm)
            else:
                tt_seq_utils.splitBam(work_dir, genome, slurm, threads)
            
        #get scale factors from yeast spike-in
        if scaleFactors == True:
            tt_seq_utils.scaleFactors(script_dir, work_dir, slurm)
            
        #create BigWig files
        bigwig = args["bigwig"]
        if bigwig == True:
            tt_seq_utils.ttSeqBigWig(work_dir, threads, tt_seq_settings, genome, slurm)
            
        #get differential signal    
        deseq2 = args["deseq2"]
        if deseq2 == True:
            tt_seq_utils.DESeq2(script_dir, genome)
            
        readRatio = args["readRatio"]
        if readRatio == True:
            if slurm == True:
                tt_seq_utils.readRatio(work_dir, script_dir, genome, slurm)
            else:
                tt_seq_utils.readRatio(work_dir, script_dir, genome, slurm, threads)
            
        if metagene != None:
            if slurm == True:
                if metagene == "ngsplot":
                    tt_seq_utils.ngsPlotSlurm(work_dir,genome)
                elif metagene == "deeptools":
                    pass
            else:
                pass
                
    
    def geneSymConv(args, script_dir):
        conversion = args["conversion"]
        gene_list = args["input"]
        out_file = args["output"]

        subprocess.call(["Rscript",
                         os.path.join(script_dir, "R", "convert_genesymbols.R"),
                         conversion,
                         gene_list,
                         out_file])


    #check for whitespace in directory
    utils.check_whitespace(work_dir)

    #execute selected module
       
    if args["module"] == "crispr":
        crispr(args, script_dir)
    elif args["module"] == "rna-seq":
        rna_seq(args, script_dir)
    elif args["module"] == "chip-seq":
        chip_seq(args, script_dir, args["module"])
    elif args["module"] == "tt-seq":
        ttSeq(args, script_dir)
    elif args["module"] == "cutrun":
        cutrun(args, script_dir)
    elif args["module"] == "damid":
        damID(args, script_dir)
    elif args["module"] == "genesymconv":
        geneSymConv(args, script_dir)
    elif args["module"] == "subsetgtf":
        subsetGTF(args, chip_seq_settings)

if __name__ == "__main__":
    #start run timer
    start = timeit.default_timer()

    script_dir = os.path.abspath(os.path.dirname(__file__))
    work_dir = os.getcwd()

    #adds script directory to runtime for importing modules
    sys.path.append(os.path.join(script_dir,"utils"))
    import utils_crispr as crispr_utils
    import utils_general as utils
    import utils_rna_seq as rnaseq_utils
    import utils_chip_seq as chipseq_utils
    import utils_cutrun as cutrun_utils
    import utils_tt_seq as tt_seq_utils
    import utils_damid as damid_utils

    #log all command line arguments to commands.log
    utils.logCommandLineArgs(work_dir)

    ###loads available CRISPR libraries from library.yaml
    import yaml

    with open(os.path.join(script_dir,
                           "yaml",
                           "crispr-library.yaml")) as file:
        crispr_settings = yaml.full_load(file)
    crispr_library_list = list(crispr_settings.keys())

    ###loads RNA-Seq settings
    try:
        with open(os.path.join(script_dir,
                               "yaml",
                               "rna-seq.yaml")) as file:
            rna_seq_settings = yaml.full_load(file)
    except FileNotFoundError:
        sys.exit("ERROR: rna-seq.yaml not found in yaml folder. Please provide this file for further analysis.")


    rna_seq_genomeList = []
    for key, value in rna_seq_settings["FASTA"].items():
        rna_seq_genomeList.append(key)
    for key, value in rna_seq_settings["RSEM_STAR_index"].items():
        rna_seq_genomeList.append(key)
    for key, value in rna_seq_settings["STAR_index"].items():
        rna_seq_genomeList.append(key)
    rna_seq_genomeList = set(rna_seq_genomeList)

    ###loads ChIP-Seq settings
    try:
        with open(os.path.join(script_dir,
                               "yaml",
                               "chip-seq.yaml")) as file:
            chip_seq_settings = yaml.full_load(file)
    except FileNotFoundError:
        sys.exit("ERROR: chip-seq.yaml not found in yaml folder. Please provide this file for further analysis.")

    ###loads CUT & RUN settings
    try:
        with open(os.path.join(script_dir,
                               "yaml",
                               "cut-run.yaml")) as file:
            cutrun_settings = yaml.full_load(file)
    except FileNotFoundError:
        sys.exit("ERROR: cut-run.yaml not found in yaml folder. Please provide this file for further analysis.")
        
    ###loads TT-Seq settings
    try:
        with open(os.path.join(script_dir,
                               "yaml",
                               "tt-seq.yaml")) as file:
            tt_seq_settings = yaml.full_load(file)
    except FileNotFoundError:
        sys.exit("ERROR: tt-seq.yaml not found in yaml folder. Please provide this file for further analysis.")    
    
        
    ###loads DamID settings
    try:
        with open(os.path.join(script_dir,
                               "yaml",
                               "damid.yaml")) as file:
            damid_settings = yaml.full_load(file)
    except FileNotFoundError:
        sys.exit("ERROR: damid.yaml not found in yaml folder. Please provide this file for further analysis.")


    main()

    #print total run time
    stop = timeit.default_timer()
    total_time = stop - start
    ty_res = time.gmtime(total_time)
    res = time.strftime("%H:%M:%S",ty_res)
    print('Total run time: ', res)

