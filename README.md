# Prerequisites

1. Blast
2. Modeller
3. pGenTHREADER
4. HHBlits (c.f. HH-Suite)
5. Complete PDB and as fasta db
      We get pdbaa from dunbrack lab
      http://dunbrack.fccc.edu/Guoli/culledpdb_hh/pdbaa.gz
6. CATH dompain pdbs and as a fasta db
    http://download.cathdb.info/cath/releases/all-releases/v4_1_0/sequence-data/cath-domain-seqs-S100-v4_1_0.fa
    http://download.cathdb.info/cath/releases/all-releases/v4_1_0/non-redundant-data-sets/cath-dataset-nonredundant-S40-v4_1_0.pdb.tgz


# BioSerf

## Authors

S.Ward, M.Sadowski, D.Jones, D.Buchan

## BioSerf Protocol

1. Run blast over the PDBAA to find out which pdbs your sequence may hit!
TODO: THIS NEEDS CHANGED TO BLAST+

`> /scratch0/NOT_BACKED_UP/dbuchan/Applications/blast-2.2.26/bin/blastpgp -h 0.001 -d /scratch0/NOT_BACKED_UP/dbuchan/uniref/pdb_aa.fasta -o bioserf_out.bls -i example/B0R5N0.fasta`

2. Create a set of pGenTHREADER models

`> GenThreader.sh -i ~/Code/bioserf/example/B0R5N0.fasta -j B0R5N0 -s`

you need the presults (for runBioserf) and ss2 files (for HHBlits), more the ss and ss2 files to the
dir you will run HHBlits in. More the presults file to the dir you will run runBioserf in

3. Create a set of models using HHBlits (HHSuite3 btw)

`hhblits -i ../example/B0R5N0.fasta -n 3 -cpu 1 -d /scratch0/NOT_BACKED_UP/dbuchan/hhblitsdb/pdb70 -oa3m B0R5N0.hhblits.a3m`

`addss.pl B0R5N0.hhblits.a3m`

`hhsearch -realign -mact 0 -cpu 1 -E 0.001 -i B0R5N0.hhblits.a3m -d /scratch0/NOT_BACKED_UP/dbuchan/hhblitsdb/pdb70 -o B0R5N0.hhr`

note that hhmakemodel may hang if the hhr file references pdbs which are not present in the
pdb dir.

`hhmakemodel.pl -m 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 -d /scratch0/NOT_BACKED_UP/dbuchan/pdb -i B0R5N0.hhr -ts B0R5N0.hh.pdb`

4. Parse blast and output modeller mod.py files and tidy up the set of pGenTHREADER models

`> javac -cp lib/biojava-1.7.1.jar:lib/bytecode.jar:./ runBioserf.java`

`> java -cp /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/src/lib/biojava-1.7.1.jar:/cs/research/bioinf/home1/green/dbuchan/Code/bioserf/src runBioserf bioserf_out.bls ../example/B0R5N0.fasta 0.00005 ./ /scratch0/NOT_BACKED_UP/dbuchan/pdb/ /scratch0/NOT_BACKED_UP/dbuchan/uniref/pdb_aa.fasta 40 B0R5N0 ~/bin/modeller9.17/bin/mod9.17 B0R5N0.pgen.presults /cs/research/bioinf/home1/green/dbuchan/bin/modeller9.17/modlib/ /cs/research/bioinf/home1/green/dbuchan/bin/modeller9.17/lib/x86_64-intel8/ /scratch0/NOT_BACKED_UP/dbuchan/python3/bin/python /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/bin/qmodcheck_mainens /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/bin/qmodcheck /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/data/modcheckpot.dat /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/bin/tmjury3d_mq_modeller_threadsafe`

* bioserf_out.bls : the blast output from step 1
* B0R5N0.fasta : The fasta file from step 1
* 0.00005: blast evalue for models
* /pdb/: location of all the pdb files
* pdb_aa.fasta: a fasta file of all the pdb files' fasta seqs
* 40: the % overlap to build models
* B0R5N0 : A unique ID for the files to be prefixed with
* mod9.17: The location of the modeller binary  
* B0R5N0.pgen.presults: presults file from Genthreader run (you need all the pdbs too)
* modlib/ : location of modeller libraries
* lib/x86-64-intel8/ : location of modelly architecture specific libs
* python : location of the python to use
* qmodcheck_mainens : path to this Exe
* qmodcheck : path to this Exe
* modchekpot.data : path to this data file
* tmjury3d_mq_modeller_threadsafe: path to this exe

Finally this should output a file
B0R5N0.final.pdb

# DomSerf protocol

Before this you need to process the CATH domain list and CATH domall files to
build an annotated domain-list.

`> python bin/process_cath_data.py cath-domain-list-v4_2_0.txt cath-domain-boundaries-v4_2_0.txt > data/cath-domain-list-v4_2_0_annotated.txt working/cath4.1/domlib/cath_domain_tdb/ > cath-domain-list-v4_2_0_annotated.txt`

1. Run blast against the CATH db sequences and PDB

`> /scratch0/NOT_BACKED_UP/dbuchan/Applications/ncbi-blast-2.2.31+/bin/psiblast -num_threads 1 -num_alignments 1000 -outfmt 0 -num_iterations 5 -inclusion_ethresh 0.001 -db /scratch0/NOT_BACKED_UP/dbuchan/uniref/CathDomainSeqs.S100.ATOM -query ../example/B0R5N0.fasta -out B0R5N0.domserf.cath.bls`

`> /scratch0/NOT_BACKED_UP/dbuchan/Applications/ncbi-blast-2.2.31+/bin/psiblast -num_threads 1 -num_alignments 1000 -outfmt 0 -num_iterations 5 -inclusion_ethresh 0.001 -db /scratch0/NOT_BACKED_UP/dbuchan/uniref/pdb_aa.fasta -query ../example/B0R5N0.fasta -out B0R5N0.domserf.pdb.bls`

## SLURM:

1. python ~/bin/split_fasta.py all.fa all_fasta 1000

2. run_domain_blast.sh
python submitter.py /home/camp/buchand/working/genome3d/Genome3D.2017-09-05/all_fasta /home/camp/buchand/Applications/bioserf/slurm_helper_scripts/run_domain_blast.sh

2. Run parse_pdb_blast.pl

`> ../bin/parse_pdb_blast.pl /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/data/cath-domain-list.txt ../example/B0R5N0.fasta B0R5N0.domserf.pdb.bls /scratch0/NOT_BACKED_UP/dbuchan/uniref/pdb_aa.fasta . /scratch0/NOT_BACKED_UP/dbuchan/pdb/ /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/bin/reformat.pl ~/bin/modeller9.17/bin/mod9.17`

* CathDomainSummary_3.5 :  a list of cath domain IDS, CATH codes, and start stop regions
* B0R5N0.domserf.pdb.bls: PDB blast output from step 1
* pdb_aa.fasta: the PDB blast db used in step 1
* ./ : tmp dir for models and whatnot
* pdb/ : location of pdb files
* reformat.pl : location of reformat.pl
* mod9.17 : location of Modeller binary

If a good hit is found pdb files of the form NAME_start_stop.2pdb will be produced for each domain that CATH

## SLURM
sbatch run_parse_pdb.sh

## Work out which proteins need to pass through to step 3
calculate_missing.py

3. If no domain pdb files were produced (i.e. we couldn't find a PDB match which was already classified in CATH). Then we run the following.

Run domTHREADER

`> ./GenThreader.sh -i B0R5N0.fasta -j B0R5N0 -d`

## SLURM
sbatch run_domthreader.sh

4. Run runParseCathDomthreader

`> ../bin/parse_cath_domthreader.pl /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/data/cath-domain-list.txt B0R5N0.domserf.cath.bls ../example/B0R5N0.fasta ./B0R5N0.pdom.presults ./B0R5N0.pdom.align ./ B0R5N0.blastaligns B0R5N0.ssf B0R5N0.pdomaligns`

`> ../bin/DomainFinder3 -i B0R5N0.ssf -o B0R5N0.dfout`

`> ../bin/make_modeller_files.pl B0R5N0.dfout B0R5N0.blastaligns B0R5N0.pdomaligns ./B0R5N0 B0R5N0.mod_lookups ../example/B0R5N0.fasta /scratch0/NOT_BACKED_UP/dbuchan/dompdb/`

# SLURM
sbatch run_parse_cath.sh

5. Do Modelling

`export PYTHONPATH=~/bin/modeller9.17/modlib:~/bin/modeller9.17/lib/x86_64-intel8/`
`export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:~/bin/modeller9.17/lib/x86_64-intel8/`

`for i in ``ls *.py | sed -e 's/\.py//'``; do echo $i; ~/bin/modeller9.17/bin/mod9.17 $i.py; done;`

`for i in ``ls *.ali | sed -e 's/\.ali//'``; do echo $i; ../bin/rewrite_modeller.pl ./  B0R5N0.mod_lookups B0R5N0.blastaligns B0R5N0.pdomaligns ../example/B0R5N0.fasta $i.ali ../bin/reformat.pl; done;`

# SLURM
sbatch run_modeller.sh
python ~/bin/split_fasta.py ecoli.fa ecoli_fasta/all
sbatch run_rewrite.sh
