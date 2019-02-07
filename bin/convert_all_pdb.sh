#!/bin/sh
for i in `ls *.ali | sed -e 's/\.ali//'`; do echo $i; /data/bioserf/bin/rewrite_modeller.pl $1 $2 $3 $4 $5 $i.ali /data/bioserf/bin/reformat.pl; done;
