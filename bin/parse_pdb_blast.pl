#!/usr/bin/perl -w

use strict;
use FileHandle;
use English;
use Data::Dumper;
use Digest::MD5;

#readCathDomainSummary
#read_fasta
#read_pdb_blast
#get hit seq
#print modeller files based on coords from CathDomainSummary
#Possibly print out something about not finding a result

my $hCathSummary = {};
my $CathDomainSummary = $ARGV[0];
print("READING CATH DOMAIN SUMMARY DATA\n");
readCathDomainSummary();

my $query_id = '';
my $query_file_id = '';
my $query_seq = '';
my $query = $ARGV[1];
print("READING FASTA SEQ\n");
read_fasta();
print($query_id);
my $query_length = length $query_seq;

my $blast_data = $ARGV[2];
print("READING BLAST RESULTS\n");
my $hBlastData = readBlast();
#print Dumper $hBlastData;
my $pdb_fasta = $ARGV[3];
print("GETTING PDB LENGTHS\n");
getPdbLengths();

my $aData = transformData();

my $hBestHit = {};
print("FINDING BEST HITS\n");
#print Dumper $aData;
$hBestHit = findBestHit();
# print Dumper $hCathSummary;
#print Dumper $hBestHit;

my $modellerDir = $ARGV[4];
$modellerDir =~ s!/*$!/!;
my $pdbDir = $ARGV[5];
my $reformatBin = $ARGV[6];
my $modellerBin = $ARGV[7];


if(ref($hBestHit) eq 'HASH')
{
	print("FOUND BEST HIT AND MODELLING\n");
	# runModeller();
	chopDomains();
}

sub testOverlap
{
	my($dom_start, $subjct_stop)  = @ARG;

	if($dom_start < $subjct_stop)
	{
		return(1);
	}
	else
	{
		return(0);
	}
}

sub chopDomains
{
	$hBestHit->{DOMAINID} =~ /^(.{5})/;
	my $pdb_code = $1;
	$pdb_code =~ /(.{4})(.)/;
	my $chain = $2;
	$chain = uc $chain;
	my $pdb = $1;
	$pdb_code = $pdb.$chain;
	my $pdb_first_residue = getPdbFirstRes($pdb,$chain);
	my $alignment_offset = $hBestHit->{SSTART};
	my $query_offset = $hBestHit->{SSTART}-$hBestHit->{QSTART};
	#print Dumper $hCathSummary;
	#print($pdb_code."\n");
	#print Dumper $hCathSummary->{$pdb_code};
	# print Dumper $hBestHit;

	my $fhTemplates = new FileHandle($modellerDir.$query_file_id."_pdb_templates.txt","w");
	print $fhTemplates "PDB CHAIN,CATH DOMAIN,Q START,Q STOP,S START,S STOP,EVALUE\n";
	foreach my $dom_id (keys %{$hCathSummary->{$pdb_code}})
	{
		#print "hi\n";
		#First test if domain overlaps with Query
		my $overlap_test = testOverlap($hCathSummary->{$pdb_code}{$dom_id}{START},$hBestHit->{SSTOP});
		if($overlap_test == 0)
		{
			next;
		}

		#slice out alignment
		my ($qstr, $sstr) = getAlignStrings($hBestHit->{ALIGNMENT});

		#the domain start and stops as pdb chain coordinates
		my $cath_dom_start = $hCathSummary->{$pdb_code}{$dom_id}{START};
		my $cath_dom_stop  = $hCathSummary->{$pdb_code}{$dom_id}{STOP};

		#the sequences as arrays
		my $aS_chars = [];
		@$aS_chars = split //, $sstr;
		my $aQ_chars = [];
		@$aQ_chars = split //, $sstr;

		#the counter starts
		my $align_position = 0;
		my $s_position = $hBestHit->{SSTART}-1;
		my $q_position = $hBestHit->{QSTART}-1;

		my $query_domain_start = 0;
		my $query_domain_stop = 0;
		my $sbjct_domain_start = 0;
		my $sbjct_domain_stop = 0;

		my $align_segment_start = 0;
		my $align_segment_stop = 0;
		foreach my $s_char (@$aS_chars)
		{
			if(@$aS_chars[$align_position] !~ /-/)
			{
				$s_position++;
			}
			if(@$aQ_chars[$align_position] !~ /-/)
			{
				$q_position++;
			}

			#test if the s_position is within in the domain boundaries
			if($s_position >= $cath_dom_start && $s_position <= $cath_dom_stop)
			{
				if($query_domain_start == 0)
				{
					$query_domain_start = $q_position;
					$align_segment_start = $align_position;
					$sbjct_domain_start = $s_position;
				}
				$query_domain_stop = $q_position;
				$align_segment_stop = $align_position;
				$sbjct_domain_stop = $s_position;
			}
			$align_position++;
		}
		#print $query_domain_start." ".$query_domain_stop."\n";
		#print $align_segment_start." ".$align_segment_stop."\n";


		#Need to correctly handle '-' characters
		#what we need to do is run along the $qstr from the start domain position until we'be got as many
		#residues from the domain length.
		my $q_segment = substr $qstr, $align_segment_start,($align_segment_stop-$align_segment_start)+1;
		my $s_segment = substr $sstr, $align_segment_start,($align_segment_stop-$align_segment_start)+1;

		#print $q_segment."\n";
		#print $s_segment."\n";

		my $ctx = Digest::MD5->new;
    $ctx->add($query_seq);
		my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my @abbr = qw( 01 02 03 04 05 06 07 08 09 10 11 12 );
    $year+=1900;
    $mon = $abbr[$mon];
    $mday=~ s/^(\d)$/0$1/;

		my $genome3d_remarks = '';
		$genome3d_remarks.="REMARK GENOME3D NAME ".$query_id."_".$query_domain_start."_".$query_domain_stop."\n";
		$genome3d_remarks.="REMARK GENOME3D UNIPROT_ID ".$query_id."\n";
		$genome3d_remarks.="REMARK GENOME3D UNIPROT_MD5 .".$ctx->hexdigest()."\n";
		$genome3d_remarks.="REMARK GENOME3D TIMESTAMP ".$year."-".$mon."-".$mday."\n";
    $genome3d_remarks.="REMARK GENOME3D TOTAL TEMPLATES 1\n";
    $genome3d_remarks.="REMARK GENOME3D SELECTION psiblast\n";
    $genome3d_remarks.="REMARK GENOME3D SELECTION pDomTHREADER\n";
    $genome3d_remarks.="REMARK GENOME3D SELECTION DomainFinder3\n";
    $genome3d_remarks.="REMARK GENOME3D MODELLING modeller\n";
    $genome3d_remarks.="REMARK GENOME3D ALIGNMENT_SOURCE t_coffee\n";
    $genome3d_remarks.="REMARK GENOME3D START ".$query_domain_start."\n";
    $genome3d_remarks.="REMARK GENOME3D STOP ".$query_domain_stop."\n";

		$genome3d_remarks.="REMARK GENOME3D >".$query_id.":".$query_domain_start."-".$query_domain_stop."\n";
		($qstr, $sstr) = convertA2M($q_segment, $s_segment);
		my $aSegments = [];
		@$aSegments = ( $q_segment =~ m/.{1,60}/sg );
		foreach my $seg (@$aSegments)
		{
			$genome3d_remarks.="REMARK GENOME3D ".$seg."\n";
		}

		$genome3d_remarks.="REMARK GENOME3D >".$hBestHit->{DOMAINID}.":".$sbjct_domain_start."-".$sbjct_domain_stop."\n";
		@$aSegments = ( $s_segment =~ m/.{1,60}/g );
		foreach my $seg (@$aSegments)
		{
			$genome3d_remarks.="REMARK GENOME3D ".$seg."\n";
		}

		print $fhTemplates $pdb_code.",".$dom_id.",".$query_domain_start.",".$query_domain_stop.",".$sbjct_domain_start.",".$sbjct_domain_stop.",".$hBestHit->{EVAL}."\n";

		#if(-e "/cluster/project2/domserf/modeller_out/".$query_id.".B99990001.pdb")
		if(-e $modellerDir.$query_file_id.".B99990001.pdb")
		{
			my $fhOut = new FileHandle($modellerDir.$query_file_id."_".$query_domain_start."_".$query_domain_stop.".pdb","w");
			print $fhOut $genome3d_remarks;

			my $fhIn = new FileHandle($modellerDir.$query_file_id.".B99990001.pdb", "r");
			my $res_count = 0;
			my $atom_count = 0;
			my $res_id = '';
			my $res_number = 0;
			while(my $line = $fhIn->getline)
			{
				if($line =~ /^ATOM\s+(\d+)\s+.+?\s+(.{3})\s+(\d+)\s+/)
				{
					my $tmpatom_count = $1;
					my $tmpres_id = $2;
					my $tmpres_number = $3;
					if($tmpres_number >=  $query_domain_start && $tmpres_number <= $query_domain_stop)
					{
						$atom_count = $tmpatom_count;
						$res_id = $tmpres_id;
						$res_number = $tmpres_number;
						print $fhOut $line;
					}
				}
				elsif($line =~ /^TER/)
				{
					print $fhOut "TER    ".$atom_count."      ".$res_id."   ".$res_number."\n";
				}
				else
				{
					print $fhOut $line;
				}
			}

		}
	}
	$fhTemplates->close;

}

sub convertA2M
{
	my($query,$subject) = @ARG;
	my $time = time;

	my $fhOut = new FileHandle($modellerDir.$query_file_id."_".$time.".fsa","w");
	print $fhOut ">Query\n";
	print $fhOut $query."\n";
	print $fhOut ">Subjct\n";
	print $fhOut $subject."\n";
	my $cmd = $reformatBin." fas a2m ".$modellerDir.$query_file_id."_".$time.".fsa ".$modellerDir.$query_file_id."_".$time.".a2m";
	print($cmd);
	`$cmd`;
	sleep 1;

	my $rmcmd = "rm ".$modellerDir.$query_file_id."_".$time.".fsa";
	`$rmcmd`;
	$fhOut ->close;

	my $fhIn = new FileHandle($modellerDir.$query_file_id."_".$time.".a2m","r");
	my $sub_found = 0;
	while(my $line = $fhIn->getline)
	{
		chomp $line;
		if($line =~ /^>Query/)
		{
			$query = '';
		}
		elsif($line =~ /^>Subjct/)
		{
			$sub_found = 1;
			$subject = '';
		}
		else
		{
			if($sub_found == 1)
			{
				$subject.=$line;
			}
			else
			{
				$query.=$line;
			}
		}
	}
	$fhIn->close;
	$rmcmd = "rm ".$modellerDir.$local_query_id."_".$time.".a2m";
	`$rmcmd`;

	return($query,$subject);
}

sub getPdbFirstRes
{
	my($pdb,$chain) = @ARG;
	my $pdb_first_residue = '';
	print $modellerDir.$query_file_id.".B99990001.pdb\n";
	if(-e $modellerDir.$query_file_id.".B99990001.pdb")
	{
		#open the pdb file find out the offset for the chains so we can map the
		if(-e $pdbDir."pdb".$pdb.".ent")
		{
			my $fhIn = new FileHandle($pdbDir."pdb".$pdb.".ent","r");
			while(my $line = $fhIn->getline)
			{
				if($line =~ /^ATOM.{17}$chain\s+(\d+)\s+/)
				{
					$pdb_first_residue =$1;
					last;
				}
			}
		}
		#CATH domain locations properly
	}
	return($pdb_first_residue);
}

sub runModeller
{
	$hBestHit->{DOMAINID} =~ /(.{5})/;
	my $pdb_code = $1;

	#print Dumper $hCathSummary->{$pdb_code};

	my ($query_str, $sub_str) = getAlignStrings($hBestHit->{ALIGNMENT});
	print $query_str."\n";
	print $sub_str."\n";

	$pdb_code =~ /(.{4})(.)/;
	my $chain = $2;
	$chain = uc $chain;
	my $pdb = $1;
	$query_id =~ s/\|/_/g;
	print "Printing: ".$query_file_id.".ali\n";
	my $fhOut = new FileHandle($modellerDir.$query_file_id.".ali","w");
	print $fhOut ">P1;$pdb\n";
	print $fhOut "structure:".$pdb.":FIRST:".$chain.":LAST:".$chain."::::\n";
	# print $fhOut $sub_str."*\n";
  print $fhOut "*\n";

	print $fhOut ">P1;$query_id\n";
	print $fhOut "sequence:".$query_id."::::::::\n";
	print $fhOut $query_str."*\n";
	$fhOut->close;

	print "Printing: ".$query_file_id.".py\n";
	$fhOut = new FileHandle($modellerDir.$query_id.".py","w");
	print $fhOut "from modeller import *\nfrom modeller.automodel import *\n\nenv = environ()\n";
	print $fhOut "env.io.atom_files_directory = '".$pdbDir."'\n";

	print $fhOut "a = automodel(env, alnfile='".$modellerDir.$query_id.".ali',\n";
    print $fhOut "  knowns=('".$pdb."'),\n";
    print $fhOut "  sequence='".$query_id."')\n";

	print $fhOut "a.starting_model = 1\na.ending_model = 1\na.auto_align()\na.make()\nok_models = filter(lambda x: x['failure'] is None, a.outputs)\n";
	print $fhOut "key = 'molpdf'\nok_models.sort(lambda a,b: cmp(a[key], b[key]))\nm = ok_models[0]\nprint \"Top model:%s\" % m['name']\n";
	chdir $modellerDir;
	my $cmd = $modellerBin." ".$modellerDir.$query_file_id.".py";
	print STDERR $cmd."\n";
	`$cmd`;
}

sub getAlignStrings
{
	my ($alignment) = @ARG;
	my $q_str = '';
	my $s_str = '';
	my $aLines = [];
	@$aLines = split /\n/, $alignment;

	foreach my $line (@$aLines)
	{
		if($line =~ /Query\s+\d+\s+(.+)\s+\s+/)
		{
			$q_str.=$1;
		}
		if($line =~ /Sbjct\s+\d+\s+(.+)\s+\s+/)
		{
			$s_str.=$1;
		}
	}

	return($q_str,$s_str);
}


sub findBestHit
{
	foreach my $hit (@$aData)
	{
		if(($hit->{EVAL} == 0.00 || $hit->{EVAL} < 0.00005) && $hit->{QUERY_OVERLAP} >= 80)
		{
			return $hit;
			last;
		}
	}
}


sub transformData
{
	my $aArray = [];
	foreach my $dom_id (keys %$hBlastData)
	{
		if(length $hBlastData->{$dom_id}{ALIGNMENT} > 0)
		{
			@$aArray[$hBlastData->{$dom_id}{COUNT}] = $hBlastData->{$dom_id};
		}
	}
	return $aArray;
}

sub getPdbLengths
{
	my $fhIn = new FileHandle($pdb_fasta,"r");

	my $current_code = '';
	while(my $line = $fhIn->getline)
	{
		chomp $line;
		if($line =~ />(\S{6})\s/)
		{
			$current_code = $1;
			if(exists $hBlastData->{$current_code})
			{
				$hBlastData->{$current_code}{SEQ}='';
			}
		}
		else
		{
			if(exists $hBlastData->{$current_code})
			{
				$hBlastData->{$current_code}{SEQ}.=$line;
				$hBlastData->{$current_code}{LENGTH}.= length $hBlastData->{$current_code}{SEQ};

			}
		}
	}
}

sub readBlast
{
    my $fhIn = new FileHandle($blast_data,"r");

    my $passed_count = 0;
    my $found_align = 0;
    my $hData = {};
    my $current_id = '';
    while ( my $line = $fhIn->getline )
    {
        chomp $line;
        #if($line =~ /^\s\s(\d.{3}[A-Z]\d)\s+.+\s+(\S+?)\s+(\S+?)\s*$/)
        if ( $line =~ /^\s\s(\d.{3}[A-Z]\d*)\s+.+\s+(.+?)\s+(\S+?)\s*$/ )
        {
            my $id    = $1;
            my $score = $2;
            my $eval  = $3;
            #if($id !~ /3eo1G0/){next;}
            #print $id." S: ".$score." E:".$eval."\n";
            if ( $eval == 0.0 || $eval < 0.00005 ) {
                $current_id = $id;
                if ( !exists $hData->{$current_id} ) {
                    $hData->{$current_id}{SCORE}     = $score;
                    $hData->{$current_id}{EVAL}      = $eval;
                    $hData->{$current_id}{QSTART}    = 0;
                    $hData->{$current_id}{QSTOP}     = 0;
                    $hData->{$current_id}{SSTART}    = 0;
                    $hData->{$current_id}{SSTOP}     = 0;
                    $hData->{$current_id}{ALIGNMENT} = '';
                    $hData->{$current_id}{DOMAINID}  = lc $id;
                    $hData->{$current_id}{COUNT}     = $passed_count;
                    $passed_count++;
                }
            }
        }
        if ( $line =~ /^>\s(.{5,6})\s/ )
        {
            $current_id = $1;
            if ( $found_align < $passed_count
                && exists $hData->{$current_id} )
            {
                $hData->{$current_id}{ALIGNMENT_HEADER} = $line;
            }
            $found_align++;
        }
        if (   $line =~ /^(Query|Sbjct)\s/
            && $found_align <= $passed_count
            && exists $hData->{$current_id} )
        {
            $hData->{$current_id}{ALIGNMENT} .= $line . "\n";
            if ( $line =~ /^Query\s+(\d+)\s.+\s(\d+)/ )
            {
                my $tmp_start = $1;
                my $tmp_stop  = $2;
                if ( $hData->{$current_id}{QSTART} == 0 )
                {
                    $hData->{$current_id}{QSTART} = $tmp_start;
                }
                if ( $tmp_stop > $hData->{$current_id}{QSTOP} ) {
                    $hData->{$current_id}{QSTOP} = $tmp_stop;
                }
            }
            if ( $line =~ /^Sbjct\s+(\d+)\s.+\s(\d+)/
                && exists $hData->{$current_id} )
            {
                my $tmp_start = $1;
                my $tmp_stop  = $2;
                if ( $hData->{$current_id}{SSTART} == 0 ) {
                    $hData->{$current_id}{SSTART} = $tmp_start;
                }
                if ( $tmp_stop > $hData->{$current_id}{SSTOP} ) {
                    $hData->{$current_id}{SSTOP} = $tmp_stop;
                }
            }
        }
        if (   $line =~ /Positives\s=\s(\d+)\/\d+/
            && $found_align <= $passed_count
            && exists $hData->{$current_id} )
        {
            $hData->{$current_id}{POSITIVES} = $1;
            $hData->{$current_id}{QUERY_OVERLAP} =
              ( $1 / $query_length ) * 100;
        }
    }
    return ($hData);

}



sub readCathDomainSummary
{
	my $fh = new FileHandle($CathDomainSummary,"r");
	while(my $line = $fh->getline)
	{
    if($line =~ /^#/){next;}
		my $aEntries = [];
		@$aEntries = split /\s+/, $line;
		@$aEntries[0] =~ /^(.{5})/;
		my $pdb_code = $1;
		$hCathSummary->{$pdb_code}{@$aEntries[0]}{CATHCODE}=@$aEntries[1].".".@$aEntries[2].".".@$aEntries[3].".".@$aEntries[4];
		$hCathSummary->{$pdb_code}{@$aEntries[0]}{START}=@$aEntries[9];
		$hCathSummary->{$pdb_code}{@$aEntries[0]}{STOP}=@$aEntries[10];
	}
}

sub read_fasta
{
    my $fhIn = new FileHandle($query,"r");
    while(my $line = $fhIn->getline)
    {
        if($line =~ /^>\s*(.+?)\s+/)
        {
          $query_id = $1;
          $query_file_id=$query_id;
          $query_file_id =~ tr/\|/_/;
          next
        }
        chomp $line;
        $query_seq.=$line;
    }
}
