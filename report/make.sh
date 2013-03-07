#!/bin/bash

base_dir="$(dirname $(readlink -f $0))"
src_dir="$base_dir/../src"

if [ -z "$1" ]; then
    arg="report.tex"
else
    arg="$1"
fi

# print source files to HTML
for file in $(ls $src_dir/*.m); do
    vim -e +':colorscheme shine' +:TOhtml +':%s/#00ffff/#000000/g' \
        +':%s/#ffff00/#0000ff/g' +w +qa $file
    mvfrom="$file.html"
    mvto="$base_dir/$(basename $mvfrom)"
    mv $mvfrom $mvto
done

# compile and count words
pdflatex $arg &&\
clear &&\
echo "$arg: $(ps2ascii report.pdf | wc -w) words"
