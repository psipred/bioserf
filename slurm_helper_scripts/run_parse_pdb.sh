#!/bin/bash
#
#SBATCH --job-name=parse_pdb_blast
#SBATCH --output=parse_pdb_blast.out
#SBATCH --error=parse_pdb_blast.err
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#
#SBATCH --array=1-1000

FILES=($HOME/working/genome3d/Genome3D.2017-09-05/all_fasta/$1/*)
out_name=${FILES[$SLURM_ARRAY_TASK_ID-1]}
out_name=${out_name:70:-5}
pdb_out_name=$out_name'domserf.pdb.bls'
echo $out_name

#echo $HOME/working/genome3d/Genome3D.2017-09-05/ecoli_fasta/1/${FILES[$SLURM_ARRAY_TASK_ID-1]}
$HOME/Applications/bioserf/bin/parse_pdb_blast.pl $HOME/working/cath4.1/cath-domain-list-v4_2_0_annotated.txt ${FILES[$SLURM_ARRAY_TASK_ID-1]} $HOME/working/genome3d/blast_results/$pdb_out_name $HOME/working/fastadb/pdbaa $HOME/working/genome3d/pdbaa_based_models/ $HOME/working/pdb/ $HOME/Applications/bioserf/bin/reformat.pl $HOME/Applications/modeller9.19/bin/mod9.19
