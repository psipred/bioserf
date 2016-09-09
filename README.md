# Prerequisites

1. Blast
2. Modeller
3. pGenTHREADER
4. HHBlits (c.f. HH-Suite)
5. Complete PDB and as fasta db
6. CATH db and as fasta db

# BioSerf

## Authors

S.Ward, M.Sadowski, D.Jones, D.Buchan

## Protocol

1. Run blast over the PDBAA to find out which pdbs your sequence may hit!

`> /scratch0/NOT_BACKED_UP/dbuchan/Applications/blast-2.2.26/bin/blastpgp -h 0.001 -d /scratch0/NOT_BACKED_UP/dbuchan/uniref/pdb_aa.fasta -o bioserf_out.bls -i example/B0R5N0.fasta`

2. Create a set of pGenTHREADER models

`> GenThreader.sh -i ~/Code/bioserf/example/B0R5N0.fasta -j B0R5N0 -s`

you need the presults (for runBioserf) and ss2 files (for HHBlits), more the ss and ss2 files to the
dir you will run HHBlits in. More the presults file to the dir you will run runBioserf in

3. Creat a set of models using HHBlits (HHSuite3 btw)

`hhblits -i ../example/B0R5N0.fasta -n 3 -cpu 1 -d /scratch0/NOT_BACKED_UP/dbuchan/hhblitsdb/pdb70 -oa3m B0R5N0.hhblits.a3m`

`addss.pl B0R5N0.hhblits.a3m`

`hhsearch -realign -mact 0 -cpu 1 -E 0.001 -i B0R5N0.hhblits.a3m -d /scratch0/NOT_BACKED_UP/dbuchan/hhblitsdb/pdb70 -o B0R5N0.hhr`

note that hhmakemodel will hang if the hhr file references pdbs which are not present in the
pdb dir.
`hhmakemodel.pl -m 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 50 -d /scratch0/NOT_BACKED_UP/dbuchan/pdb -i B0R5N0.hhr -ts B0R5N0.hh.pdb`

4. Parse blast and output modeller mod.py files and tidy up the set of pGenTHREADER models

`> javac -cp lib/biojava-1.7.1.jar:lib/bytecode.jar:./ runBioserf.java`

`> java -cp /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/src/lib/biojava-1.7.1.jar:/cs/research/bioinf/home1/green/dbuchan/Code/bioserf/src runBioserf bioserf_out.bls ../example/B0R5N0.fasta 0.00005 ./ /scratch0/NOT_BACKED_UP/dbuchan/pdb/ /scratch0/NOT_BACKED_UP/dbuchan/uniref/pdb_aa.fasta 40 B0R5N0 ~/bin/modeller9.17/bin/mod9.17 B0R5N0.pgen.presults /cs/research/bioinf/home1/green/dbuchan/bin/modeller9.17/modlib/ /cs/research/bioinf/home1/green/dbuchan/bin/modeller9.17/lib/x86_64-intel8/ /scratch0/NOT_BACKED_UP/dbuchan/python3/bin/python /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/bin/qmodcheck_mainens /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/bin/qmodcheck /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/data/modcheckpot.dat /cs/research/bioinf/home1/green/dbuchan/Code/bioserf/bin/tmjury3d_mq_modeller_threadsafe`

bioserf_out.bls : the blast output from step 1
B0R5N0.fasta : The fasta file from step 1
0.00005: blast evalue for models
/pdb/: location of all the pdb files
pdb_aa.fasta: a fasta file of all the pdb files' fasta seqs
40: the % overlap to build models
B0R5N0 : A unique ID for the files to be prefixed with
mod9.17: The location of the modeller binary  
B0R5N0.pgen.presults: presults file from Genthreader run (you need all the pdbs too)
modlib/ : location of modeller libraries
lib/x86-64-intel8/ : location of modelly architecture specific libs
python : location of the python to use
qmodcheck_mainens : path to this Exe
qmodcheck : path to this Exe
modchekpot.data : path to this data file
tmjury3d_mq_modeller_threadsafe: path to this exe

Finally this should output a file
B0R5N0.final.pdb

# DomSerf protocol
runParsePdbBlast(strParsePdbBlast, strCathSummary, fPfilt.getCanonicalPath(), fBlastOutPdb.getCanonicalPath(), strPdbDB, strTmp, strPdb, strReformat, strModellerBin);
            boolean pdbFound = gatherModels(strTmp);
            boolean domsFound = false;
            if(pdbFound == false)
            {
                runParseCathDomthreader(strParseCathDom, strCathSummary, fBlastOutCath.getCanonicalPath(), fPfilt.getCanonicalPath(), fPResults.getCanonicalPath(), fAlign.getCanonicalPath(), tempNameRoot);

                if(fSSF.exists() && fSSF.length() > 0)
                {
                    runDomainFinder(strDomainFinder, fSSF.getCanonicalPath(), fSSFProcessed.getCanonicalPath());
                    if(fSSFProcessed.exists())
                    {
                        FileInputStream fileinputstream = new FileInputStream(fSSFProcessed.getCanonicalPath());
                        int numberBytes = fileinputstream.available();
                        byte byteArray[] = new byte[numberBytes];
                        int size_of_read = fileinputstream.read(byteArray);
                        localJob.UpdateStatus(StatusClass.DomainTemplates, StatusCode.Unknown, "Found Domain templates" , byteArray);
                    }

                    //make_modeller_files.pl e9aa184a-754a-4caf-bb02-23f333aa1c69.ssf_out b4f69346-fb49-4852-88f0-d85858c190dc.blastaligns b4f69346-fb49-4852-88f0-d85858c190dc.pdomaligns /webdata/tmp/NewPredServer/b4f69346-fb49-4852-88f0-d85858c190dc b4f69346-fb49-4852-88f0-d85858c190dc.mod_lookups b4f69346-fb49-4852-88f0-d85858c190dc.pfilt /webdata/data/dompdb/
                    runMakeModeller(strMakeModeller, fSSFProcessed.getCanonicalPath(), fCathAligns.getCanonicalPath(), fPDomAligns.getCanonicalPath(), tempNameRoot, fModLookups.getCanonicalPath(), fPfilt.getCanonicalPath(), strDomPDB );
                    runModeller(rootUUID, strTmp, strModellerBin);
                    //rewrite_modeller.pl /webdata/tmp/NewPredServer b4f69346-fb49-4852-88f0-d85858c190dc.mod_lookups b4f69346-fb49-4852-88f0-d85858c190dc.blastaligns b4f69346-fb49-4852-88f0-d85858c190dc.pdomaligns b4f69346-fb49-4852-88f0-d85858c190dc.pfilt b4f69346-fb49-4852-88f0-d85858c190dc_1.ali
                    runRewrite(strRewriteModeller, strTmp, fModLookups.getCanonicalPath(), fCathAligns.getCanonicalPath(), fPDomAligns.getCanonicalPath(), fPfilt.getCanonicalPath(), rootUUID);
                    domsFound = gatherModels(strTmp);
                }
                else
                {
                    //TODO: return a fault if the sff file disappeared after ParseCathDom ran.
                    //I simply have no idea how this could happen but I suppose 1 in a million times it might
                }

            }
            else
            {   //93
                if(fPDBTemplates.exists())
                {
                    FileInputStream fileinputstream = new FileInputStream(fPDBTemplates.getCanonicalPath());
                    int numberBytes = fileinputstream.available();
                    byte byteArray[] = new byte[numberBytes];
                    int size_of_read = fileinputstream.read(byteArray);
                    localJob.UpdateStatus(StatusClass.PDBTemplates, StatusCode.Unknown, "Found PDB templates" , byteArray);
                }
                localJob.UpdateStatus(StatusClass.domserfCathMatch, StatusCode.JobRunning, "A full length PDB chain with CATH domains was matched", strError);
            }
