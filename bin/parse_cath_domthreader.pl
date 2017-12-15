#!/bin/perl -w

use strict;
use FileHandle;
use DirHandle;
use English;
use Data::Dumper;

my $CathDomainSummary = $ARGV[0];
my $tmpPath = $ARGV[5];
my $dom_lib = $ARGV[9];
my $line_counts = `wc -l $dom_lib/*`;
print($line_counts);
exit();
# print $tmpRoot.".blastaligns\n";
# my $fhBlastAlignOut = new FileHandle($tmpRoot.".blastaligns","w");
# my $fhSSF = new FileHandle($tmpRoot.".ssf","w");
# my $fhAligns = new FileHandle($tmpRoot.".pdomaligns","w");
my $fhBlastAlignOut = new FileHandle($tmpPath.$ARGV[6],"w");
my $fhSSF = new FileHandle($tmpPath.$ARGV[7],"w");
my $fhAligns = new FileHandle($tmpPath.$ARGV[8],"w");
my $hBlastData ={};
my $length = 0;
my $hPDomData = {};

my $hCathSummary = {};
print("Reading CATH summary\n");
readCathDomainSummary();
# print Dumper $hCathSummary;
$hBlastData ={};
$length = 0;
my $ID = '';
$hPDomData = {};

#if($i != 15){next;}
my $file = $ARGV[1];
my $fasta = $ARGV[2];
my $prot_id = '';
if($fasta =~ /\//)
{
  if($fasta =~ /^.+\/(.+)\.pfilt/){$prot_id = $1;}
	elsif($fasta =~ /^.+\/(.+)\.fasta/){$prot_id = $1;}
  elsif($fasta =~ /^.+\/(.+)\.fa/){$prot_id = $1;}
  elsif($fasta =~ /^.+\/(.+)\.input/){$prot_id = $1;}
	else
	{
		print "COULD NOT ID FASTA FILE\n";
		exit;
	}
}
else
{
	if ($fasta =~ /^(.+)\.fasta/ ){$prot_id = $1;}
	elsif ($fasta =~ /^(.+)\.pfilt/){$prot_id = $1;}
	elsif ($fasta =~ /^(.+)\.fa/){$prot_id = $1;}
  elsif($fasta =~ /^(.+)\.input/){$prot_id = $1;}
	else
	{
		print "COULD NOT ID FASTA FILE\n";
		exit;
	}
}
#exit;
if(-e $file)
{
	print $fhBlastAlignOut $file."\n";
	$length = read_fasta($fasta);
	$hBlastData = read_blast_data($file);
	print $fhBlastAlignOut "-------\n";
}
my $pgen_file = $ARGV[3];
my $align_file = $ARGV[4];
#exit;

if(-e $pgen_file && -e $align_file)
{
  # print($pgen_file);
	$hPDomData = read_pdom_data($pgen_file);
}
print Dumper $hBlastData;
print Dumper $hPDomData;
remove_low_overlaps();

print_ssf();
print_alignments();


sub print_alignments
{
	foreach my $id (keys %$hBlastData)
	{
		print $fhBlastAlignOut $ID." ".$id."\n";
		print $fhBlastAlignOut $hBlastData->{$id}{ALIGNMENT_HEADER}."\n";
		print $fhBlastAlignOut $hBlastData->{$id}{ALIGNMENT}."\n";
	}
	foreach my $id (keys %$hPDomData)
	{
		print $fhAligns $ID." ".$id."\n";
		print $fhAligns $hPDomData->{$id}{ALIGNMENT_HEADER}."\n";
		print $fhAligns $hPDomData->{$id}{ALIGNMENT}."\n";
	}
}

sub print_ssf
{
	foreach my $id (keys %$hBlastData)
	{
		print $fhSSF $ID." ".$id." 0 0 0 0 0 0 0 0 ".$hBlastData->{$id}{EVAL}." 0.00 0.00 1 ".$hBlastData->{$id}{START}.":".$hBlastData->{$id}{STOP}."\n";
	}
	foreach my $id (keys %$hPDomData)
	{
		print $fhSSF $ID." ".$id." 0 0 0 0 0 0 0 0 ".$hPDomData->{$id}{PVAL}." 0.00 0.00 1 ".$hPDomData->{$id}{START}.":".$hPDomData->{$id}{STOP}."\n";
	}
}

sub remove_low_overlaps
{
	foreach my $id (keys %$hBlastData)
	{
		my $domId=$hBlastData->{$id}{DOMAINID};
		if(! exists $hCathSummary->{$domId})
		{
      print($domId."\n");
			delete $hBlastData->{$id};
			next;
		}
		my $align_length = $hBlastData->{$id}{STOP}-$hBlastData->{$id}{START};
		#print $domId."\n";
		my $ratio = $align_length/$hCathSummary->{$domId}{LENGTH};
    print $ratio."\n";
		if($ratio < 0.4)
		{
				delete $hBlastData->{$id};
		}
	}

	foreach my $id (keys %$hPDomData)
	{
		my $domId=$hPDomData->{$id}{DOMAINID};
		if(! exists $hCathSummary->{$domId})
		{
			delete $hPDomData->{$id};
			next;
		}
		my $align_length =$hPDomData->{$id}{STOP}-$hPDomData->{$id}{START};
		my $ratio = $align_length/$hCathSummary->{$domId}{LENGTH};
    print $ratio."\n";
		if($ratio < 0.4)
		{
			delete $hPDomData->{$id};
		}

	}
}

$fhBlastAlignOut->close;

sub read_fasta
{
	my ($ffile) = @ARG;
	$ffile =~ s/bls/fasta/;
	my $fhIn = new FileHandle($ffile,"r");
	my $seq = '';
	while(my $line = $fhIn->getline)
	{
		if($line =~ /^>(.+)/)
		{$ID = $1;
		 $ID = substr($ID, 0, index($ID, ' '));
     next;}
		chomp $line;
		$seq.=$line;
	}
	my $length = length $seq;

  # before we return in the input files starts with an
  # we'll use the first
  # portion as the name instead
  if($ffile =~ /(.{8})-.{4}-.{4}-.{4}-.{12}/)
  {
    $ID=$1;
  }
  #print($ID);
	return($length);
}

sub read_blast_data
{
	my ($file) = @ARG;
	my $fhIn = new FileHandle($file,"r");

	my $found_alignments = 0;
	my $passed_count = 0;
	my $found_align = 0;
	my $hData = {};
	my $current_id = '';
	my $current_length = 0;
	my $align_count = 0;
	my $pdb_id ='';
	my $get_align = 0;
	while(my $line = $fhIn->getline)
	{
		chomp $line;

		if($line =~ /^>\sdomain\|(.{7})\|.+/)
		{
			$current_id = "CDOM|".$1.":";
			$pdb_id = $1;
			$found_alignments = 1;
			$align_count=-1;
			$get_align = 0;
		}
		if($line =~ /^>\scath\|4_1_0\|(.{7})\/.+/)
		{
			$current_id = "CDOM|".$1.":";
			$pdb_id = $1;
			$found_alignments = 1;
			$align_count=-1;
			$get_align = 0;
		}

		if($line =~ /^Length=(\d+)/)
		{
			$current_length = $1;
		}
		if($line =~ /^\sScore\s=\s+(.+)\s+bits\s.+,\s+Expect\s=\s+(.+),\sMethod/)
		{
			my $score = $1;
			my $eval = $2;
			if($eval == 0.00 || $eval < 0.00005)
			{

				$align_count++;

				$hData->{$current_id.$align_count}{START} = 0;
				$hData->{$current_id.$align_count}{STOP} = 0;
				$hData->{$current_id.$align_count}{ALIGNMENT} = '';
				$hData->{$current_id.$align_count}{SCORE} = $score;
				$hData->{$current_id.$align_count}{EVAL} = $eval;
				$hData->{$current_id.$align_count}{DOMAINID} = $pdb_id;
				$hData->{$current_id.$align_count}{DOMAINLENGTH} = $current_length;
				$hData->{$current_id.$align_count}{ALIGNMENT_HEADER} = $current_id.$align_count;

				$get_align = 1;
			}
			else
			{
				$get_align = 0;
			}
		}


		if($line =~ /^(Query|Sbjct)\s/ && $get_align == 1)
		{
			$hData->{$current_id.$align_count}{ALIGNMENT} .= $line."\n";
			if($line =~ /^Query\s+(\d+)\s+.+\s+(\d+)/)
			{
				#print $line;
				my $tmp_start = $1;
				my $tmp_stop = $2;
				if($hData->{$current_id.$align_count}{START} == 0)
				{
					$hData->{$current_id.$align_count}{START} =$tmp_start;
				}
				if($tmp_stop > $hData->{$current_id.$align_count}{STOP})
				{
					$hData->{$current_id.$align_count}{STOP} = $tmp_stop;
				}

			}
		}

	}
	#print Dumper $hData;
	#exit;
	return($hData);
}

sub read_pdom_data
{
	my ($pfile) = @ARG;

	$pfile =~ s/bls/pdom.presults/;
  my $fhIn = new FileHandle($pfile,"r");
	my $hit_count = 0;
	my $hData = {};
	#print $pfile."\n";
	while(my $line = $fhIn->getline)
	{
		#if($line =~ /^CERT|^HIGH/)
    if($line =~ /^CERT|^HIGH|^MEDIUM/)
    {
			$hit_count++;
			chomp $line;
			my $aEntries = [];
			@$aEntries = split /\s+/, $line;
			my $current_id = "PDOM|".@$aEntries[11].":".$hit_count;
			$hData->{$current_id}{SCORE} = @$aEntries[1];
			$hData->{$current_id}{PVAL} = @$aEntries[2];
			$hData->{$current_id}{START} = @$aEntries[9];
			$hData->{$current_id}{STOP} = @$aEntries[10];
			$hData->{$current_id}{DOMAINID} = @$aEntries[11];

		}
	}

	$pfile =~ s/presults/align/;
	$fhIn = new FileHandle($pfile,"r");
	my $align_count = 0;
	my $current_id = '';
	my $id = '';
	while(my $line = $fhIn->getline)
	{
		if($line =~ />>>\sAlignment\swith\s(.+):/)
		{
			$align_count++;
			$id = $1;
			$current_id = "PDOM|".$id.":".$align_count;
			if(exists $hData->{$current_id})
			{
				$hData->{$current_id}{ALIGNMENT_HEADER} = ">".$current_id;
				$hData->{$current_id}{ALIGNMENT} = '';
			}
		}
    # print $prot_id."\n";
		if(exists $hData->{$current_id})
		{
			if($line =~ /^$prot_id|^Query/)
			{
				$hData->{$current_id}{ALIGNMENT}.=$line;
			}
			if($line =~ /^$id/)
			{
				$hData->{$current_id}{ALIGNMENT}.=$line;
			}
		}
	}
  #print Dumper $hData;
	return($hData);
}

sub readCathDomainSummary
{
	my $fh = new FileHandle($CathDomainSummary,"r");
	while(my $line = $fh->getline)
	{
		if($line =~ /^#/){next;}
		my $aEntries = [];
		@$aEntries = split /\s+/, $line;
		@$aEntries[0] =~ /^(.{7})/;
		my $pdb_code = $1;
		$hCathSummary->{$pdb_code}{CATHCODE}=@$aEntries[1].".".@$aEntries[2].".".@$aEntries[3].".".@$aEntries[4];
		$hCathSummary->{$pdb_code}{START}=@$aEntries[12];
		$hCathSummary->{$pdb_code}{STOP}=@$aEntries[13];
    $hCathSummary->{$pdb_code}{LENGTH}=@$aEntries[13]-@$aEntries[12];
	}
}
