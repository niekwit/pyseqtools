#!/usr/bin/env python3

import os
import sys
import subprocess
import argparse
import multiprocessing
import timeit
import time
import pkg_resources
from platform import python_version

def checkPythonPackages(): #check for required python packages; installs if absent
    #check if Python >= 3.7 (required for GitPython)
    try:
        version = float(python_version().rsplit(".",1)[0])
        if version < 3.7:
            sys.exit("ERROR: please update to at least Python 3.5")
    except:
        pass

    #check packages
    required = {"shyaml", "pyyaml", "pandas", "numpy",
                "matplotlib", "seaborn", "multiqc",
                "cutadapt", "tqdm","gseapy",
                "matplotlib-venn", "pysam", "deeptools",
                "macs3", "pybedtools"}
    installed = {pkg.key for pkg in pkg_resources.working_set}
    missing = required - installed
    if missing:
        python = sys.executable
        print("Installing missing required Python3 packages")
        try:
            install_command = [python, '-m', 'pip', 'install', *missing]
            subprocess.check_call(install_command, stdout=subprocess.DEVNULL)
        except:
            sys.exit("ERROR: package installation failed")
    else:
        pass

def main():

    #create top-level parser
    parser = argparse.ArgumentParser(prog = 'pyseqtools.py',
                                     description = "Analysis pipelines for a variety of NGS data, and related tools")
    subparsers = parser.add_subparsers(help ='Available pipelines/tools:',
                                       dest ='module')

    #create subparsers for each module, i.e. functionality

    #create subparser for crispr screen analysis commands
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

    # create the parser for RNA-Seq-analysis
    parser_rnaseq = subparsers.add_parser('rna-seq',
                                          description = "Analysis pipeline for RNA-Seq experiments",
                                          help='RNA-Seq analysis')

    parser_rnaseq.add_argument("-t", "--threads",
                               required = False,
                               default = 1,
                               metavar = "<int>",
                               help = "Number of CPU threads to use (default is 1). Use max to apply all available CPU threads. For Salmon 8-12 threads are optimal")
    parser_rnaseq.add_argument("-r", "--reference",
                               required = False,
                               choices = rna_seq_genomeList,
                               help = "Reference genome")
    parser_rnaseq.add_argument("-a", "--align",
                               required = False,
                               choices = ["salmon","hisat2"],
                               default = "salmon",
                               help = "Choose aligner. Default is Salmon.")
    parser_rnaseq.add_argument("-p", "--pvalue",
                               required = False,
                               metavar = "<P value>",
                               default = 0.001,
                               help = "Set P value cut off (default is 0.001)")
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
                                        description = "Analysis pipeline for ChIP-Seq experiments",
                                        help='ChIP-Seq analysis')

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
                             default = 'hg19',
                             help = "Choose reference genome (default is hg19)")
    parser_chip.add_argument("-a", "--align",
                             choices = ["hisat2",
                                        "bwa-mem",
                                        "bwa-aln"],
                             const = "hisat2",
                             nargs = "?",
                             required = False,
                             help = "Create BAM files using HISAT2 or BWA (mem or aln")
    parser_chip.add_argument("-d", "--deduplication",
                             required = False,
                             action = 'store_true',
                             help = "Perform deduplication of BAM files")
    parser_chip.add_argument("--downsample",
                             required = False,
                             action = 'store_true',
                             help = "Perform downsampling of BAM files")
    parser_chip.add_argument("-b", "--bigwig",
                             required = False,
                             action = 'store_true',
                             help = "Create BigWig files")
    parser_chip.add_argument("--qc",
                             required = False,
                             action = 'store_true',
                             help = "Perform QC analysis of BAM files")
    parser_chip.add_argument("-p", "--peaks",
                             required = False,
                             action = 'store_true',
                             help = "Call and annotate peaks with MACS3/HOMER")
    parser_chip.add_argument("--metagene",
                             required = False,
                             action = 'store_true',
                             help = "Generate metageneplots and heatmaps with plotProfile (deepTools)")
    parser_chip.add_argument("--skip-fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Skip FastQC/MultiQC")


    #create subparser for CUT&RUN analysis commands
    parser_cutrun = subparsers.add_parser('cutrun',
                                          description = "Analysis pipeline for CUT & RUN experiments",
                                          help = 'CUT & RUN analysis')
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

    #create subparser for DamID
    parser_damid = subparsers.add_parser('damid',
                                          description = "Analysis pipeline for DamID (wrapper for https://github.com/owenjm/damidseq_pipeline)",
                                          help = 'DamID analysis')
    parser_damid.add_argument("-t", "--threads",
                           required = False,
                           default = 1,
                           type = str,
                           help = "<INT> number of CPU threads to use (default is 1). Use max to apply all available CPU threads")
    parser_damid.add_argument("-g", "--genome",
                            required = False,
                            default = 'hg19',
                            help = "Choose reference genome (default is hg19)")
    parser_damid.add_argument("--skip-fastqc",
                               required = False,
                               action = 'store_true',
                               default = False,
                               help = "Skip FastQC/MultiQC")

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

    #create subparser for subsetting GTF files
    parser_subsetgtf = subparsers.add_parser('subsetgtf',
                                             description = "Tools for subsetting GTF files according to input gene list",
                                             help = 'Subset GTF files for selected genes')

    parser_subsetgtf.add_argument("-l", "--list",
                             required = True,
                             help = "Input gene list file. Each gene symbol should be on a new line")
    parser_subsetgtf.add_argument("-i", "--input",
                             required = True,
                             help = "GTF file to be subsetted (must be specified with genome name, loaded from chip-seq.yaml)")
    parser_subsetgtf.add_argument("-o", "--output",
                             required = True,
                             help = "Output file name")



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
        crispr_utils.check_index(crispr_settings, crispr_library, script_dir, work_dir)

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

        ####set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        threads = args["threads"]
        if threads == "max":
            threads = max_threads

        ### Check md5sums
        utils.checkMd5(work_dir)

        ###Run FastQC/MultiQC
        file_extension = utils.get_extension(work_dir)
        skip_fastqc = args["skip_fastqc"]
        if not skip_fastqc:
            utils.fastqc(script_dir, work_dir, threads, file_extension)
        else:
            print("Skipping FastQC/MultiQC analysis")

        ###Set species variable
        reference = args["reference"]
        if "hg" in reference or reference == "gencode-v35":
            species = "human"
        elif "mm" in reference or reference == "gencode.vM1.pc_transcripts":
            species = "mouse"


        ###trim and align
        pvalue = args["pvalue"]
        align = args["align"]
        if align.lower() == "salmon":
            utils.trim(script_dir, threads, work_dir)
            salmon_index = rna_seq_settings["salmon_index"][reference]
            gtf = rna_seq_settings["salmon_gtf"][reference]
            fasta = rna_seq_settings["FASTA"][reference]
            rnaseq_utils.salmon(salmon_index,
                                str(threads),
                                work_dir,
                                gtf,
                                fasta,
                                script_dir,
                                rna_seq_settings,
                                reference)
            rnaseq_utils.plotMappingRate(work_dir)
            rnaseq_utils.plotPCA(work_dir, script_dir)
            rnaseq_utils.diff_expr(work_dir, gtf, script_dir, species, pvalue)
            rnaseq_utils.plotVolcano(work_dir)
        elif align.lower() == "hisat2":
            rnaseq_utils.trim(threads, work_dir)
            #hisat2()

        go = args["go"]

        pvalue = args["pvalue"]


        if go == True:
            gene_sets=args["gene_sets"]
            rnaseq_utils.geneSetEnrichment(work_dir, pvalue, gene_sets)

    def chip_seq(args, script_dir):

        #set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        threads = args["threads"]
        if threads == "max":
            threads = max_threads

        #Check md5sums
        utils.checkMd5(work_dir)


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


    def cutrun(args, script_dir):
        #set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        threads = args["threads"]
        if threads == "max":
            threads = max_threads

        #Check md5sums
        utils.checkMd5(work_dir)


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
        #Check md5sums
        utils.checkMd5(work_dir)
        
        #set thread count for processing
        max_threads = str(multiprocessing.cpu_count())
        threads = args["threads"]
        if threads == "max":
            threads = max_threads
        
        genome = args["genome"]
        
        utils.trim(script_dir, threads, work_dir)
        damid_utils.damID(script_dir, work_dir, threads, genome, damid_settings)
        damid_utils.bedgraph2BigWig(script_dir, work_dir, damid_settings, genome)
                

    def geneSymConv(args, script_dir):
        conversion = args["conversion"]
        gene_list = args["input"]
        out_file = args["output"]

        subprocess.call(["Rscript",
                         os.path.join(script_dir, "R", "convert_genesymbols.R"),
                         conversion,
                         gene_list,
                         out_file])


    def subsetGTF(args, chip_seq_settings):
        genome = args["input"]
        gene_list = args["list"]
        gtf = chip_seq_settings["gtf"][genome]
        output = os.path.join(os.path.dirname(gtf) + args["output"])
        subprocess.call(["grep", "-w", "-f", gene_list, gtf, ">", output])



    #execute selected module
    if args["module"] == "crispr":
        crispr(args, script_dir)
    elif args["module"] == "rna-seq":
        rna_seq(args, script_dir)
    elif args["module"] == "chip-seq":
        chip_seq(args, script_dir)
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
    import utils_damid as damid_utils

    #log all command line arguments to commands.log
    utils.logCommandLineArgs(work_dir)

    #check if required Python packages are available
    checkPythonPackages()

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
