#!/bin/bash
if [ -z "$1" ]; then
    arg="report.tex"
else
    arg="$1"
fi

pdflatex $arg &&\
clear &&\
echo "$arg: $(ps2ascii report.pdf | wc -w) words"
