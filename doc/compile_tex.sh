#!/bin/bash

# This script can be used to submit all commands required to compile a .tex file which uses both
# sage and bibtex.
#
# Input parameters:
#   $1 path to tex root.
#   $2 name of tex file (without file extension).
#
# example usage:
#   sh compile_tex.sh bases bezier_bases

if [[ $# -eq 2 ]]; then
    cd $1

    pdflatex $2.tex
    bibtex   $2
    sage     $2.sagetex.sage
    pdflatex $2.tex
    pdflatex $2.tex
else
    echo "This script requires two input arguments: (path, file name)."
    exit 1
fi
