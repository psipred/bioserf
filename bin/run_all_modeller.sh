#!/bin/sh
for i in `ls *.py | sed -e 's/\.py//'`; do echo $i; /opt/modeller9.17/bin/mod9.17 $i.py; done;
