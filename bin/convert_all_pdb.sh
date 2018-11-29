#!/bin/sh
for i in `ls *.ali | sed -e 's/\.ali//'`; do echo $i; /opt/bioserf/bin/rewrite_modeller.pl $1 $2 $3 $4 $5 $i.ali /opt/bioserf/bin/reformat.pl; done;
