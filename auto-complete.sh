#!/bin/bash

#This script enables autocompletion for the module and CRISPR library options (-l/--library) in the pyseqtools command
#Add the following line to your .bashrc file: source /path/to/pyseqtools/auto-complete.sh

#finds CRISPR libary names
SCRIPT_DIR=$(find $HOME -type d -name "pyseqtools")
lib_list=$(cat "$SCRIPT_DIR/yaml/crispr-library.yaml" | shyaml keys | tr "\n" " ")
stat_list="mageck bagel2"
module_list='crispr rna-seq chip-seq cutrun'

#enables autocompletion of `-l/--library` flag for CRISPR screen analysis



_module()
{
    local opts
    opts="crispr rna-seq chip-seq cutrun"
    case $COMP_CWORD in
        1)
            COMPREPLY=( $(compgen -W "${opts}" -- "${COMP_WORDS[COMP_CWORD]}") )
            ;;

    esac
    return 0
}

#Assign the auto-completion function _get for our command get.
complete -F _module pyseqtools.py 

function crisprLibs()
{
case $3 in
	-l) COMPREPLY=($(compgen -W "$lib_list" "${COMP_WORDS[$COMP_CWORD]}"));;
	--library) COMPREPLY=($(compgen -W "$lib_list" "${COMP_WORDS[$COMP_CWORD]}"));;
	-a) COMPREPLY=($(compgen -W "$stat_list" "${COMP_WORDS[$COMP_CWORD]}"));;
	--analysis) COMPREPLY=($(compgen -W "$stat_list" "${COMP_WORDS[$COMP_CWORD]}"));;
esac
}

complete -F crisprLibs pyseqtools.py


#enables autocomletion of pyseqtools modules

function get {
    pyseqtools.py $1
} 
