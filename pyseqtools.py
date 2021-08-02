#!/usr/bin/env python3

import os
import sys
import argparse
import multiprocessing
import timeit
import time
import yaml


def main():
    
    #create top-level parser
    parser = argparse.ArgumentParser(prog = 'pyseqtools.py',
                                     description = "Analysis of a variety of NGS data")
    subparsers = parser.add_subparsers(help ='Available analyses:',
                                       dest ='module')
    
    #create subparsers for each module, i.e. functionality
    
    #create subparser for crispr screen analysis commands
    parser_crispr = subparsers.add_parser('crispr', 
                                          help = 'CRISPR screens')

    parser_crispr.add_argument("-l", "--library",
                               required = "--csv2fasta" not in sys.argv,
                               choices = crispr_library_list,
                               help = "CRISPR library")
    parser_crispr.add_argument("-t", "--threads",
                               required = False, 
                               default = 1,
                               metavar="<int>",
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
                               default = os.path.join(script_dir,"core-essential-genes.csv"), #"path to Hart list"
                               help = "Essential gene list (default is Hart et al 2015 Cell)")
    parser_crispr.add_argument("--csv2fasta",
                               required = False,
                               metavar = "<CSV file>",
                               default = None,help = "Convert CSV file with sgRNA names and sequences to fasta format. The first column should contain sgRNA names and the second sgRNA sequences (headers can be anything).")
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
     
    # create the parser for RNA-Seq-analysis
    parser_rnaseq = subparsers.add_parser('rna-seq', 
                                          help='RNA-Seq')
            
    parser_rnaseq.add_argument("-t", "--threads",
                               required = False,
                               default = 1,
                               metavar = "<int>",
                               help = "Number of CPU threads to use (default is 1). Use max to apply all available CPU threads. For Salmon 8-12 threads are optimal")
    parser_rnaseq.add_argument("-r", "--reference",
                               required = False,
                               choices = rna_seq_genomeList,
                               help = "Reference genome")
    parser_rnaseq.add_argument("-s", "--species",
                               required = True,
                               choices = ["mouse","human"],
                               help = "Set species.")
    parser_rnaseq.add_argument("-a", "--align",
                               required = False,
                               choices = ["salmon","hisat2"],
                               default = "salmon",
                               help = "Choose aligner. Default is Salmon.")
    parser_rnaseq.add_argument("-p", "--pvalue",
                               required = False,
                               metavar = "<P value>",
                               default = 0.001,
                               help = "Set P value cut off (default is 0.001")
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
    parser_rnaseq.add_argument("--skip-fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Skip FastQC/MultiQC")
    
    #create subparser for ChIP-Seq analysis commands
    parser_chip = subparsers.add_parser('chip-seq', 
                                          help='ChIP-Seq')
    
    parser_chip.add_argument("-t", "--threads",
                             required = False,
                             default = 1,
                             help = "<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    parser_chip.add_argument("-r", "--rename", 
                             required = False, 
                             action = 'store_true', 
                             help = "Rename fq files")
    parser_chip.add_argument("-f", "--fastqc", 
                             required = False,
                             action = 'store_true',
                             help = "Perform FASTQC")
    parser_chip.add_argument("-g", "--genome", 
                             required = False,
                             default = 'hg19',
                             help = "Choose reference genome (default is hg19)")
    parser_chip.add_argument("-a", "--align",
                             choices = ['hisat2',
                                        'bwa'], 
                             required = False,
                             help = "Trim and align raw data to index using HISAT2 or BWA")
    parser_chip.add_argument("-d", "--deduplication", 
                             required = False, 
                             action = 'store_true',
                             help = "Perform deduplication of BAM files")
    parser_chip.add_argument("-s", "--downsample", 
                             required = False, 
                             action = 'store_true',
                             help = "Perform downsampling of BAM files")
    parser_chip.add_argument("-b", "--bigwig", 
                             required = False,
                             action = 'store_true',
                             help = "Create BigWig files")
    parser_chip.add_argument("-q", "--qc", 
                             required = False, 
                             action = 'store_true',
                             help = "Perform QC analysis of BAM files")
    parser_chip.add_argument("-p", "--peaks", 
                             required = False, 
                             action = 'store_true',
                             help = "Call and annotate peaks")
    parser_chip.add_argument("-n", "--ngsplot", 
                             required = False, 
                             action = 'store_true',
                             help = "Generate metageneplots and heatmaps with ngs.plot")
    
    #create subparser for CUT&RUN analysis commands
    parser_cutrun = subparsers.add_parser('cutrun', 
                                          help='CUT & RUN')
    
    parser_cutrun.add_argument("-t", "--threads",
                             required = False,
                             default = 1,
                             help = "<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
        
    #create dictionary with command line arguments
    args = vars(parser.parse_args())
    
        
    def crispr(args):
        ###check if software requirements are met
        utils.checkPythonPackages()
        
        #csv to fasta conversion
        csv = args["csv2fasta"]
        if csv is not None:
            if not os.path.isfile(csv):
                sys.exit("ERROR: invalid file path given")
            else:
                crispr_utils.csv2fasta(csv,script_dir)
    
        ###check if more software requirements are met
        crispr_utils.checkDeps(script_dir)
        utils.checkDeps(script_dir)
    
        ###set thread count for processing
        threads = utils.set_threads(args)
        
        ###Check md5 checksums
        utils.checkMd5(work_dir)
        
        ###run modules based on parsed arguments:
        ##rename files
        rename = args["rename"]
        if rename == True:
            utils.rename(work_dir)
    
        #determine file extension raw data
        file_extension = utils.get_extension(work_dir)
    
        ##Run FastQC/MultiQC
        skip_fastqc = args["skip_fastqc"]
        if not skip_fastqc:
            utils.fastqc(work_dir,threads,file_extension,exe_dict)
        else:
            print("Skipping FastQC/MultiQC analysis")
    
        ##count reads
        #check if bowtie2 index is build for CRISPR library
        crispr_library = args["library"]
        crispr_utils.check_index(crispr_libraries, crispr_library, script_dir, work_dir)
    
        #check if file with just guide names exists
        crispr_utils.guide_names(crispr_libraries,crispr_library)
    
        #count sgRNAs
        mismatch = args["mismatch"]
        crispr_utils.count(crispr_libraries,
                     crispr_library,
                     mismatch,
                     threads,
                     script_dir,
                     work_dir)
    
        #plot alignment rates
        crispr_utils.plot_alignment_rate(work_dir)
    
        #plot sample coverage (read count / library size)
        crispr_utils.plot_coverage(work_dir,crispr_libraries,crispr_library)
    
        #join count files
        if not utils.file_exists(os.path.join(work_dir,
                                    "count",
                                    'counts-aggregated.tsv')):
            utils.join_counts(work_dir,crispr_libraries,crispr_library)
        #normalise read count table
        if not utils.file_exists(os.path.join(work_dir,
                                    "count",
                                    "counts-aggregated-normalised.csv")):
            crispr_utils.normalise(work_dir)
    
        ##run library analysis
        crispr_utils.lib_analysis(work_dir,crispr_libraries,crispr_library,script_dir)
        crispr_utils.gcBias(work_dir,crispr_libraries,crispr_library)
    
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
                crispr_utils.remove_duplicates(work_dir)
                crispr_utils.convert4bagel(work_dir,crispr_libraries,crispr_library)
                crispr_utils.bagel2(work_dir, script_dir, exe_dict, fdr)
        
        #run essential gene list comparison
        essential_genes = args["essential_genes"]
        crispr_utils.essentialGenes(work_dir, analysis, essential_genes, fdr)
    
        if go == True:
            gene_sets = args["gene_sets"]
            crispr_utils.goPython(work_dir,
                            fdr,
                            crispr_libraries,
                            crispr_library,
                            analysis,
                            gene_sets)
    
    def rna_seq(args):
        ####set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        threads = args["threads"]
        if threads == "max":
            threads=max_threads
    
        ### Check md5sums
        utils.checkMd5(work_dir)
    
        ###Run FastQC/MultiQC
        file_extension=utils.getExtension(work_dir)
        skip_fastqc = args["skip_fastqc"]
        if not skip_fastqc:
            utils.fastqc(work_dir,threads,file_extension)
        else:
            print("Skipping FastQC/MultiQC analysis")
    
        ###Set species variable
        species = args["species"]
        
        ###trim and align
        pvalue = args["pvalue"]
        align = args["align"]
        if align.lower() == "salmon":
            utils.trim(threads,work_dir)
            salmon_index=rna_seq["salmon_index"]["gencode-v35"]
            gtf=rna_seq["salmon_gtf"]["gencode-v35"]
            fasta=rna_seq["FASTA"]["gencode-v35"]
            utils.salmon(salmon_index,str(threads),work_dir,gtf,fasta,script_dir,rna_seq)
            utils.plotMappingRate(work_dir)
            utils.plotPCA(work_dir,script_dir)
            utils.diff_expr(work_dir,gtf,script_dir,species,pvalue)
            utils.plotVolcano(work_dir)
        elif align.lower() == "hisat2":
            utils.trim(threads,work_dir)
            #hisat2()
    
    go = args["go"]
    
    try:
        pvalue = args["pvalue"]
    except KeyError:
        pass
	
    if go == True:
        gene_sets=args["gene_sets"]
        utils.geneSetEnrichment(work_dir,pvalue,gene_sets)
    
    def chip_seq(args):
        pass
    
    #execute selected module
    if args["module"] == "crispr":
        crispr(args)
    elif args["module"] == "rna-seq":
        rna_seq(args)
    elif args["module"] == "chip-seq":
        chip_seq(args)
    
    
if __name__ == "__main__":
    #start run timer
    start = timeit.default_timer()
    
    script_dir = os.path.abspath(os.path.dirname(__file__))
    work_dir = os.getcwd()
    
    #adds script directory to runtime for importing modules
    sys.path.append(script_dir)
    import utils_crispr as crispr_utils
    import utils_general as utils
    import utils_rna_seq as rnaseq_utils
    import utils_chip_seq as chipseq_utils
    import utils_cutrun as cutrun_utils
    
    
    
    ###loads available CRISPR libraries from library.yaml
    with open(os.path.join(script_dir,"yaml","crispr-library.yaml")) as file:
        crispr_libraries=yaml.full_load(file)
    crispr_library_list=list(crispr_libraries.keys())
    
    ###loads RNA-Seq settings
    try:
        with open(os.path.join(script_dir,
                               "yaml",
                               "rna-seq.yaml")) as file:
            rna_seq = yaml.full_load(file)
    except FileNotFoundError:
        print("ERROR: rna-seq.yaml not found in yaml folder. Please provide this file for further analysis.")
        sys.exit()
    
        
    rna_seq_genomeList = []
    for key, value in rna_seq["FASTA"].items():
        rna_seq_genomeList.append(key)

    ###loads ChIP-Seq settings
    
   
    
    main()
    
    
    #print total run time
    stop = timeit.default_timer()
    total_time = stop - start
    ty_res = time.gmtime(total_time)
    res = time.strftime("%H:%M:%S",ty_res)
    print('Total run time: ', res)