#!/bin/bash

report_dir="$(dirname $(readlink -f $0))"
code_dir="$(readlink -f $report_dir/../src)"
src_dir="$report_dir/src"

tex_file="$src_dir/report.tex"
matlab_file="$src_dir/code.m"
pdf_file="$src_dir/$(basename $tex_file .tex).pdf"
report_file="$report_dir/$(basename $tex_file .tex)-full.pdf"

# print source files to HTML and then to PDF
if [ -f $matlab_file ]; then
    rm $matlab_file
fi
for f in $(find $code_dir -maxdepth 1 -type f -name '*.m' | sort); do
    cat $f >> $matlab_file
    echo >> $matlab_file
    echo >> $matlab_file
done
for f in $(find $code_dir -mindepth 2 -type f -name '*.m' -not -path '*/.*/*' \
           | sort); do
    cat $f >> $matlab_file
    echo >> $matlab_file
    echo >> $matlab_file
done
head -n -2 $matlab_file > $matlab_file.tmp
mv $matlab_file.tmp $matlab_file
vim -e +':colorscheme shine' +:TOhtml +':%s/#00ffff/#000000/g' \
    +':%s/#ffff00/#0000ff/g' +w +qa $matlab_file
html2ps --colour $matlab_file.html > $matlab_file.ps
ps2pdf -sPAPERSIZE=a4 -dEmbedAllFonts=true -dOptimize=true \
    -dAutoRotatePages=/None $matlab_file.ps $matlab_file.pdf

# compile report
pdflatex -output-directory $src_dir $tex_file
pdflatex -output-directory $src_dir $tex_file
pdfjoin --rotateoversize 'false' $pdf_file $matlab_file.pdf -o $report_file

# count words
echo "$report_file: $(ps2ascii $pdf_file | wc -w) words"
