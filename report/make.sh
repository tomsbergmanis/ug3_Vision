#!/bin/bash

base_dir="$(dirname $(readlink -f $0))"
src_dir="$base_dir/../src"

if [ -z "$1" ]; then
    arg="report.tex"
else
    arg="$1"
fi

# print source files to HTML
mfile="$src_dir/.all.m"
vim -e +':colorscheme shine' +:TOhtml +':%s/#00ffff/#000000/g' \
    +':%s/#ffff00/#0000ff/g' +w +qa $mfile
mv "$mfile.html" "$base_dir/code.m.html"

# compile and count words
pdflatex $arg &&\
clear &&\
echo "$arg: $(ps2ascii report.pdf | wc -w) words"

pdfjoin report.pdf code.pdf -o report-full.pdf
