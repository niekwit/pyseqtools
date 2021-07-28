#!/bin/bash

#This script enables autocompletion for the CRISPR library options (-l) in crispr-pipeline.sh
#Add the following line to your .bashrc file: source /path/to/CRISPR-tools/auto-complete.sh

#finds CRISPR libary names
SCRIPT_DIR=$(find $HOME -type d -name "CRISPR-tools")
lib_list=$(cat "$SCRIPT_DIR/library.yaml" | shyaml keys | tr "\n" " ")
stat_list="mageck bagel2"

#enables autocompletion of `-l` flag
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


#write function for autocompletion of modules ($1)

#function modules()
#{
#
#}
