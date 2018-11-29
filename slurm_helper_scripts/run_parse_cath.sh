#!/bin/bash
#
#SBATCH --job-name=parse_cath_domth
#SBATCH --output=messages/parse_cath_domth_%J.out
#SBATCH --error=messages/parse_cath_domth_%J.err
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#
#SBATCH --array=1-1000

path=$HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$1/
size=${#path}
SUBSET=$1
FILES=($HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$SUBSET/*)
out_name=${FILES[$SLURM_ARRAY_TASK_ID-1]}
out_name=${out_name:size:-6}
pdb_out_name=$out_name'.domserf.pdb.bls'
echo $out_name

# echo $HOME/working/genome3d/Genome3D.2017-09-05/ecoli_fasta/1/${FILES[$SLURM_ARRAY_TASK_ID-1]}
if [ -f $HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$SUBSET/$out_name'.fasta' ]; then
  $HOME/Applications/bioserf/bin/parse_cath_domthreader.pl $HOME/working/cath4.2/cath-domain-list-v4_2_0_annotated.txt $HOME/working/genome3d/blast_results/$out_name'.domserf.cath.bls' ${FILES[$SLURM_ARRAY_TASK_ID-1]} $HOME/working/genome3d/domthreader_output/$out_name'.pdom.presults' $HOME/working/genome3d/domthreader_output/$out_name'.pdom.align' $HOME/working/genome3d/parse_cath_domth_output/ $out_name'.blastaligns' $out_name'.ssf' $out_name'.pdomaligns' $HOME/working/cath4.1/domlib/cath_domain_tdb/
  $HOME/Applications/bioserf/bin/DomainFinder3 -i $HOME/working/genome3d/parse_cath_domth_output/$out_name'.ssf' -o $HOME/working/genome3d/parse_cath_domth_output/$out_name'.dfout'
  $HOME/Applications/bioserf/bin/make_modeller_files.pl $HOME/working/genome3d/parse_cath_domth_output/$out_name'.dfout' $HOME/working/genome3d/parse_cath_domth_output/$out_name'.blastaligns' $HOME/working/genome3d/parse_cath_domth_output/$out_name'.pdomaligns' $HOME/working/genome3d/parse_cath_domth_output/$out_name $HOME/working/genome3d/parse_cath_domth_output/$out_name'.mod_lookups' $HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$SUBSET/$out_name'.fasta' $HOME/working/cath4.1/dompdb/
fi
