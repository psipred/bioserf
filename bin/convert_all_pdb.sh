#!/bin/sh
for i in `ls *.ali | sed -e 's/\.ali//'`; do echo $i; /opt/bioserf/bin/rewrite_modeller.pl $TMP/$ID/ $I2 $I4 $I5 $I1 $i.ali /opt/bioserf/bin/reformat.pl; done;
