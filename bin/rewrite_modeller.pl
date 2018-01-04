#!/usr/bin/perl -w

use strict;
use FileHandle;
use DirHandle;
use English;
use Data::Dumper;
use Digest::MD5;

my $tmp = $ARGV[0];

my $alidir = $tmp;
my $modeller_coords = $ARGV[1];
my $fastadir = $tmp;
my $modeldir = $tmp;
my $rewritedir = $tmp;

my $blast_aligns = $ARGV[2];
my $pdom_aligns = $ARGV[3];
my $fasta_file = $ARGV[4];
my $ali_file = $ARGV[5];
my $reformat_bin = $ARGV[6];


#my $hLookup = read_lookup();
my $hLookup = make_lookup();
my $hCoords = read_coords();

print Dumper $hLookup;
print Dumper $hCoords;

read_ali();

sub make_lookup
{
	#my $dhIn = new DirHandle($fasta_files);
	my $hData = {};
	my $id = '';;

	# $fasta_file =~ /^(.+)\.pfilt/;
	# $id = $1;

	  if($fasta_file =~ /\//)
	  {
	    if($fasta_file =~ /^.+\/(.+)\.pfilt/){$id = $1;}
	  	elsif($fasta_file =~ /^.+\/(.+)\.fasta/){$id = $1;}
			elsif($fasta_file =~ /^.+\/(.+)\.fa/){$id = $1;}
			elsif($fasta_file =~ /^.+\/(.+)\.input/){$id = $1;}
	  	else
	  	{
	  		print "COULD NOT GET ID FROM FASTA FILE\n";
	  		exit;
	  	}
	  }
	  else
	  {
	  	if ($fasta_file =~ /^(.+)\.fasta/ ){$id = $1;}
	  	elsif ($fasta_file =~ /^(.+)\.pfilt/){$id = $1;}
			elsif ($fasta_file =~ /^(.+)\.fa/){$id = $1;}
			elsif ($fasta_file =~ /^(.+)\.input/){$id = $1;}
	  	else
	  	{
	  		print "COULD NOT GET ID FROM FASTA FILE\n";
	  		exit;
	  	}
	  }
	my $fhIn = new FileHandle($fasta_file,"r");

	LOOP: while(my $line = $fhIn->getline)
	{
			chomp $line;
			# print $line;
			if($line =~ /^>(.+?)\s/)
			{
				$hData->{$1} = $id;
				#print $fhLookup $1." : ".$id."\n";
				$fhIn->close;
				last LOOP;
			}
		}
	#}
	return($hData);
}


sub read_ali
{
		$ali_file =~ /(.+)\.ali/;
		my $number = $1;
		print $alidir."/".$ali_file."\n";
		my $fhIn = new FileHandle($alidir."/".$ali_file,"r");
		print "Processing : ".$ali_file."\n";
		my $modelid = '';
		my $q_seq = '';
		my $s_seq = '';
		my $q_start = 0;
		my $q_stop = 0;
		my $s_stop = 0;
		my $s_start = 0;
		my $targetid = '';
		#also need to get the start and stop regions of the query seq
		#get the query data
		my $line = $fhIn->getline;
		if($line =~ />P1;(.+)/)
		{
			$modelid = $1;
      my $lookupid = $modelid;
      $lookupid =~ s/\|/_/g;
      if(exists $hCoords->{$modelid})
			{
				$q_start = $hCoords->{$modelid}{START};
				$q_stop  = $hCoords->{$modelid}{STOP};
			}
			else
			{
				print STDERR $modelid." NO COORD LOOKUP\n";
				next;
			}
		}
		$line = $fhIn->getline;
		#if($line =~ /sequence:::::::::/)
		#{
		#	$q_start = $1;
		#	$q_stop = $2;
		#}
		$line = $fhIn->getline;
		$q_seq = $line;
		$q_seq =~ s/\*$//;
		chomp $q_seq;

		#get the target data
		 $line = $fhIn->getline;
		if($line =~ />P1;(.+)/)
		{
			$targetid = $1;
		}
		$line = $fhIn->getline;
		if($line =~ /StructureX:.+:(\d+):.:(\d+):.::::/)
		{
			$s_start = $1;
			$s_stop = $2;
		}
		#print $s_start." ".$s_stop."\n";
		$line = $fhIn->getline;
		$s_seq = $line;
		$s_seq =~ s/\*$//;
		chomp $s_seq;

		#if($modelid !~ /A0SXL3/){next;}
		# print $modelid."\n";
		my $modelfile = $modelid.".B99990001.pdb";
		my $model = '';
		print $modeldir.$modelfile."\n";
		if(-e $modeldir.$modelfile)
		{
			my $fhMod = new FileHandle($modeldir.$modelfile,"r");

			while(my $line2 = $fhMod->getline)
			{
				$model.=$line2;
			}

		}
		else
		{
			print STDERR $modelid." NO MODEL\n";
			next;
		}

		reprintModel($model, $modelid, $targetid,$s_seq,$q_seq,$q_start,$q_stop,$s_start,$s_stop);
		#exit;

}

sub reprintModel
{
	my($model, $modelid, $targetid,$s_seq,$q_seq,$q_start,$q_stop,$s_start,$s_stop) = @ARG;
	$modelid =~ s/_.+$//;

	my $full_seq = get_full_seq($modelid);
	#print $full_seq;
	my $ctx = Digest::MD5->new;
    $ctx->add($full_seq);
	my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
    my @abbr = qw( 01 02 03 04 05 06 07 08 09 10 11 12 );
    $year+=1900;
    $mon = $abbr[$mon];
    $mday=~ s/^(\d)$/0$1/;

	my $genome3d_remarks = '';
	$genome3d_remarks.="REMARK GENOME3D NAME ".$modelid."_".$q_start."_".$q_stop."\n";
	$genome3d_remarks.="REMARK GENOME3D UNIPROT_ID ".$modelid."\n";
	$genome3d_remarks.="REMARK GENOME3D UNIPROT_MD5 .".$ctx->hexdigest()."\n";
	$genome3d_remarks.="REMARK GENOME3D TIMESTAMP ".$year."-".$mon."-".$mday."\n";
    $genome3d_remarks.="REMARK GENOME3D TOTAL TEMPLATES 1\n";
    $genome3d_remarks.="REMARK GENOME3D SELECTION psiblast\n";
    $genome3d_remarks.="REMARK GENOME3D SELECTION pDomTHREADER\n";
    $genome3d_remarks.="REMARK GENOME3D SELECTION DomainFinder3\n";
    $genome3d_remarks.="REMARK GENOME3D MODELLING modeller\n";
    $genome3d_remarks.="REMARK GENOME3D ALIGNMENT_SOURCE psiblast\n";
	$genome3d_remarks.="REMARK GENOME3D ALIGNMENT_SOURCE pDomTHREADER\n";
	$genome3d_remarks.="REMARK GENOME3D START ".$q_start."\n";
    $genome3d_remarks.="REMARK GENOME3D STOP ".$q_stop."\n";

	$genome3d_remarks.="REMARK GENOME3D >".$modelid.":".$q_start."-".$q_stop."\n";
	my ($sstr, $qstr) = convertA2M($s_seq, $q_seq, $modelid);
	my $aSegments = [];
	@$aSegments = ( $qstr =~ m/.{1,60}/sg );
	foreach my $seg (@$aSegments)
	{
		$genome3d_remarks.="REMARK GENOME3D ".$seg."\n";
	}

	$genome3d_remarks.="REMARK GENOME3D >".$targetid.":".$s_start."-".$s_stop."\n";
	@$aSegments = ( $sstr =~ m/.{1,60}/g );
	foreach my $seg (@$aSegments)
	{
		$genome3d_remarks.="REMARK GENOME3D ".$seg."\n";
	}
	#print $genome3d_remarks;
	if($q_start == 0 || $q_stop == 0){print STDERR $modelid." BAD BOUNDARIES\n";return;}
	my $fhOut = new FileHandle($rewritedir.$modelid."_".$q_start."_".$q_stop.".pdb","w");
	print "Printing: ".$modelid."_".$q_start."_".$q_stop."\n";
	print $fhOut $genome3d_remarks;
	my $aLines = [];
	@$aLines = split /\n/, $model;
	my $res_count = 0;
	my $atom_count = 0;
	my $res_id = '';
	my $res_number = 0;

	foreach my $line (@$aLines)
	{
		if($line =~ /^ATOM\s+(\d+)\s+.+?\s+(.{3})\s+(\d+)\s+/)
		{
			my $tmpatom_count = $1;
			my $tmpres_id = $2;
			my $tmpres_number = $3;
			#if($tmpres_number >=  $q_start && $tmpres_number <= $q_stop)
			#{
				$atom_count = $tmpatom_count;
				$res_id = $tmpres_id;
				$res_number = $tmpres_number;
				print $fhOut $line."\n";
			#}
		}
		elsif($line =~ /^TER/)
		{
			print $fhOut "TER    ".$atom_count."      ".$res_id."   ".$res_number."\n";
		}
		else
		{
			print $fhOut $line."\n";
		}
	}

}

sub convertA2M
{
	my($query,$subject,$modelid) = @ARG;
	my $time = time;
	my $fhOut = new FileHandle($alidir.$modelid."_".$time.".fsa","w");
	print $fhOut ">Query\n";
	print $fhOut $query."\n";
	print $fhOut ">Subjct\n";
	print $fhOut $subject."\n";

	my $this_modelid = $modelid;
	my $cmd = $reformat_bin." fas a2m '".$alidir."/".$this_modelid."_".$time.".fsa' '".$alidir."/".$this_modelid."_".$time.".a2m'";
  print $cmd."\n";
	#exit;
	`$cmd`;
	sleep 1;
#exit;
	my $rmcmd = "rm '".$alidir."/".$modelid."_".$time.".fsa'";
	`$rmcmd`;
	$fhOut ->close;

	my $fhIn = new FileHandle($alidir."/".$modelid."_".$time.".a2m","r");
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
	$rmcmd = "rm '".$alidir."/".$modelid."_".$time.".a2m'";
	`$rmcmd`;

	return($query,$subject);
}




sub get_full_seq
{
	#my($modelid) = @ARG;

	my $seq = '';
	#if(exists $hLookup->{$modelid})
	#{
		my $fhIn = new FileHandle($fasta_file);
		$fhIn->getline;
		while(my $line = $fhIn->getline)
		{
			chomp $line;
			$seq.=$line;
		}
	#}
	return $seq;
}

sub read_coords
{
	my $fhLookup = new FileHandle($modeller_coords, "r");
	my $hData = {};
	while(my $line = $fhLookup->getline)
	{
		chomp $line;
		$line =~ /(.+)\s(\d+)\s(\d+)/;
    my $id = $1;
    $id =~ s/\|/_/g;
		$hData->{$id}{START} = $2;
		$hData->{$id}{STOP} = $3;
	}
	return($hData);
}
