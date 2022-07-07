#!/bin/bash

#This script enables autocompletion for the module and CRISPR library options (-l/--library, -a/--analysis and positional argument) in the pyseqtools command
#Add the following line to your .bashrc file: source /path/to/pyseqtools/auto-complete.sh

#finds CRISPR libary names
SCRIPT_DIR=$(find $HOME -type d -name "pyseqtools")
lib_list=$(cat "$SCRIPT_DIR/yaml/crispr-library.yaml" | shyaml keys | tr "\n" " ")
stat_list="mageck bagel2"



function _complete()
{
case $3 in
	-l) COMPREPLY=($(compgen -W "$lib_list" "${COMP_WORDS[$COMP_CWORD]}"));;
	--library) COMPREPLY=($(compgen -W "$lib_list" "${COMP_WORDS[$COMP_CWORD]}"));;
	-a) COMPREPLY=($(compgen -W "$stat_list" "${COMP_WORDS[$COMP_CWORD]}"));;
	--analysis) COMPREPLY=($(compgen -W "$stat_list" "${COMP_WORDS[$COMP_CWORD]}"));;
esac

local opts
opts="crispr rna-seq chip-seq cutrun damid tt-seq genesymconv subsetgtf"
case $COMP_CWORD in
    1)
        COMPREPLY=( $(compgen -W "${opts}" -- "${COMP_WORDS[COMP_CWORD]}") )
        ;;

esac
return 0
}

complete -F _complete pyseqtools.py



