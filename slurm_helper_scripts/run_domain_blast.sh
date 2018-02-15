#!/bin/bash
#
#SBATCH --job-name=cath4_blasts
#SBATCH --output=cath4_blasts.out
#SBATCH --error=cath4_blasts.err
#
#SBATCH --ntasks=1
#SBATCH --time=6:00:00
#
#SBATCH --array=1-1000

FILES=($HOME/working/genome3d/Genome3D.2017-09-05/ecoli_fasta/$1/*)
out_name=${FILES[$SLURM_ARRAY_TASK_ID-1]}
out_name=${out_name:70:-5}
cath_out_name=$out_name'domserf.cath.bls'
pdb_out_name=$out_name'domserf.pdb.bls'
# echo $SLURM_ARRAY_TASK_ID$out_name
echo srun $HOME/Applications/ncbi-blast-2.7.1+/bin/psiblast -num_threads 1 -num_alignments 1000 -outfmt 0 -num_iterations 5 -inclusion_ethresh 0.001 -db $HOME/working/fastadb/cath-domain-seqs-S100-v4_1_0.fa -query ${FILES[$SLURM_ARRAY_TASK_ID-1]} -out $HOME/working/genome3d/blast_results/$cath_out_name
echo srun $HOME/Applications/ncbi-blast-2.7.1+/bin/psiblast -num_threads 1 -num_alignments 1000 -outfmt 0 -num_iterations 5 -inclusion_ethresh 0.001 -db $HOME/working/fastadb/pdbaa -query ${FILES[$SLURM_ARRAY_TASK_ID-1]} -out $HOME/working/genome3d/blast_results/$pdb_out_name
