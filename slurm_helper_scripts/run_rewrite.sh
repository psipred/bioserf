#!/bin/bash
#
#SBATCH --job-name=run_modeller
#SBATCH --output=run_modeller.out
#SBATCH --error=run_modeller.err
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#

export PYTHONPATH=$HOME/Applications/modeller9.19/modlib/:$HOME/Applications/modeller9.19/lib/x86_64-intel8/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$HOME/Application/modeller9.19/lib/x86_64-intel8/
cd $HOME/working/genome3d/parse_cath_domth_output/
for i in `ls *.py | sed -e 's/\.py//'`; do echo $i; $HOME/Applications/modeller9.19/bin/mod9.19 $i.py; done;
for i in `ls *.ali | sed -e 's/\.ali//'`; do echo $i; $HOME/Applications/bioserf/bin/rewrite_modeller.pl ./ B0R5N0.mod_lookups B0R5N0.blastaligns B0R5N0.pdomaligns ../example/B0R5N0.fasta $i.ali ../bin/reformat.pl; done;`
sp_P05052_APPY_ECOLI
