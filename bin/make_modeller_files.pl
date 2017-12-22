#!/usr/bin/perl -w
#
# Read in the domainFindered domains and output a modeller .py and aln file for
# each domain
#

use strict;
use FileHandle;
use DirHandle;
use English;
use Data::Dumper;

my $ssf_file = $ARGV[0];
my $blast_aligns = $ARGV[1];
my $pdom_aligns = $ARGV[2];
my $tmpDir = $ARGV[3];
my $fasta_files = $tmpDir;
#my $lookup = "/cluster/project2/domserf/lookup.txt";
my $dompdb = $ARGV[6];
my $ali_dir = $tmpDir;
my $py_dir = $tmpDir;
my $mod_lookups = $ARGV[4];
my $fasta_file = $ARGV[5];
my $fhLook = new FileHandle($mod_lookups,"w");
my $id = '';
my $hSsf = read_ssf();
my $hLookup = make_lookup();
#my $hLookup = read_lookup();
#print Dumper $hLookup;
#print Dumper $hSsf;
#exit:
my $hDomCount = {};
#aligns are converted to PIR format and output
my $output_count = 0;
read_blast_aligns();
read_pdom_aligns();
#print $output_count." FILES OUTPUT\n";


sub make_lookup
{
	#my $dhIn = new DirHandle($fasta_files);
	my $hData = {};

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

sub read_pdom_aligns
{
	my $fhIn = new FileHandle($pdom_aligns, "r");
	my $hData = {};
	my $align_name = 'START';
	my $get_align = 0;
	my $pdb_id = '';
	my $q_seq = '';
	my $s_seq = '';

	while(my $line = $fhIn->getline)
	{
		chomp $line;
		# print("PDB:".$pdb_id."\n");
    if($line =~ /^(.+\sPDOM\|.+)/)
    {
        my $name = $1;
        $align_name = $name;
        # print($align_name."\n");
        $hData->{$align_name}{'q_seq'} = '';
        $hData->{$align_name}{'s_seq'} = '';
        if($align_name =~ /PDOM\|(.+):\d+/)
        {
          $pdb_id=$1;
          # print($pdb_id."\n");
        }
    }
    if($line =~ /^Query\s+(.+)/)
  	{
  		$hData->{$align_name}{'q_seq'}.= $1;
  	}
  	if($line =~ /^$pdb_id\s+(.+)/)
  	{
  		$hData->{$align_name}{'s_seq'}.= $1;
  	}
	}
  foreach my $thisid (keys $hData)
  {
    if(exists $hSsf->{$thisid})
  	{
      $output_count++;
      print_ali($thisid,0,0,0,0,$hData->{$thisid}{'q_seq'},$hData->{$thisid}{'s_seq'},1);
    	print_py($thisid,1);
    }
  }
}

sub read_blast_aligns
{
	my $fhIn = new FileHandle($blast_aligns, "r");
	my $hData = {};
	my $align_name = 'START';
	my $get_align = 0;
	my $pdb_id = '';
	my $q_seq = '';
	my $s_seq = '';

	while(my $line = $fhIn->getline)
	{
		chomp $line;
		# print("PDB:".$pdb_id."\n");
    if($line =~ /^(.+\sCDOM\|.+)/)
    {
        my $name = $1;
        $align_name = $name;
        #print($align_name."\n");
        $hData->{$align_name}{'q_seq'} = '';
        $hData->{$align_name}{'q_start'} = 0;
        $hData->{$align_name}{'q_stop'} = 0;
        $hData->{$align_name}{'s_seq'} = '';
        $hData->{$align_name}{'s_start'} = 0;
        $hData->{$align_name}{'s_stop'} = 0;
        if($align_name =~ /CDOM\|(.+):\d+/)
        {
          $pdb_id=$1;
          # print($pdb_id."\n");
        }
    }
    if($line =~ /^Query\s+(\d+)\s+(.+?)\s+(\d+)/)
  	{
      my $start = $1;
    	my $stop = $3;
  		$hData->{$align_name}{'q_seq'}.= $2;
      if($hData->{$align_name}{'q_start'} == 0)
    	{
    		$hData->{$align_name}{'q_start'} = $start;
    	}
    	if($stop > $hData->{$align_name}{'q_stop'})
    	{
    		$hData->{$align_name}{'q_stop'} = $stop;
    	}
  	}
  	if($line =~ /^Sbjct\s+(\d+)\s+(.+?)\s+(\d+)/)
  	{
      my $start = $1;
    	my $stop = $3;
  		$hData->{$align_name}{'s_seq'}.= $2;
      if($hData->{$align_name}{'s_start'} == 0)
    	{
    		$hData->{$align_name}{'s_start'} = $start;
    	}
    	if($stop > $hData->{$align_name}{'s_stop'})
    	{
    		$hData->{$align_name}{'s_stop'} = $stop;
    	}
  	}
	}
  foreach my $thisid (keys $hData)
  {
    if(exists $hSsf->{$thisid})
  	{
      $output_count++;
      my $start = read_domain_start($thisid);
      print_ali($thisid,$hData->{$thisid}{'q_start'},
                $hData->{$thisid}{'q_stop'},
                ($hData->{$thisid}{'s_start'}+$start)-1,
                ($hData->{$thisid}{'s_stop'}+$start)-1,
                $hData->{$thisid}{'q_seq'},$hData->{$thisid}{'s_seq'},0);
    	print_py($thisid,0);
    }
  }
}

sub read_blast_aligns_alt
{
	my $fhIn = new FileHandle($blast_aligns, "r");
	my $hData = {};
	my $get_align = 0;
	my $q_start = 0;
	my $q_stop = 0;
	my $s_start = 0;
	my $s_stop = 0;
	my $q_seq = '';
	my $s_seq = '';
	my $align_name = 'START';

	while(my $line = $fhIn->getline)
	{
		chomp $line;
		if($line =~ /^(.+\sCDOM\|.+)/)
		{
			my $name = $1;
			if($align_name !~ /START/ && $get_align==1)
			{
				print $align_name."\n";
				$output_count++;
				my $start = read_domain_start($align_name);
				print_ali($align_name,$q_start,$q_stop,($s_start+$start)-1,($s_stop+$start)-1,$q_seq,$s_seq,0);
				print_py($align_name,0);
				$get_align = 0;
				$q_start = 0;
				$q_stop = 0;
				$s_start = 0;
				$s_stop = 0;
				$q_seq = '';
				$s_seq = '';
				$get_align = 0;
				#exit;
			}
			$align_name = $name;

			if(exists $hSsf->{$align_name})
			{
				$get_align = 1;
				#if($line !~ /^Q86XP0\sCATHDOM\|1cjyB02:0/){$get_align=0}
			}
			else
			{
				$get_align = 0;
			}

		}
		if($line =~ /^Query\s+(\d+)\s+(.+?)\s+(\d+)/ && $get_align ==1)
		{
			my $start = $1;
			my $stop = $3;
			$q_seq.=$2;
			if($q_start == 0)
			{
				$q_start = $start;
			}
			if($stop > $q_stop)
			{
				$q_stop = $stop;
			}
		}
		if($line =~ /^Sbjct\s+(\d+)\s+(.+?)\s+(\d+)/ && $get_align ==1)
		{
			my $start = $1;
			my $stop = $3;
			$s_seq.=$2;
			if($s_start == 0)
			{
				$s_start = $start;
			}
			if($stop > $s_stop)
			{
				$s_stop = $stop;
			}
		}
	}
	if($align_name !~ /START/ && $get_align==1)
			{
				print $align_name."\n";
				$output_count++;
				my $start = read_domain_start($align_name);
				print_ali($align_name,$q_start,$q_stop,($s_start+$start)-1,($s_stop+$start)-1,$q_seq,$s_seq,0);
				print_py($align_name,0);
				$get_align = 0;
				$q_start = 0;
				$q_stop = 0;
				$s_start = 0;
				$s_stop = 0;
				$q_seq = '';
				$s_seq = '';
				$get_align = 0;
				#exit;
			}
}

sub read_domain_start
{
	my($align_name) = @ARG;
	$align_name =~ /^(.+?)\sCDOM\|(.{4})(.{1})(.{2}):\d+/;

	my $seq_id = $1;
	my $pdb_id = $2.$3.$4;

	my $fhIn = new FileHandle($dompdb.$pdb_id,"r");
	my $line = $fhIn->getline;
	$line =~ /^ATOM.{7}\s\s.{4}.{3}.{2}(.{4})/;
	my $start = $1;
	$start =~ s/\s+//;
	#print $start."\n";
	$fhIn->close;
	return($start);
}

sub print_py
{
	my($align_name,$pdom_ctrl) = @ARG;
	my $dom_type = 'CDOM';
	if($pdom_ctrl == 0)
	{
		$dom_type = 'CDOM';
	}
	elsif($pdom_ctrl == 1)
	{
		$dom_type = 'PDOM';
	}
	$align_name =~ /^(.+?)\s$dom_type\|(.{4})(.{1})(.{2}):\d+/;
	my $seq_id = $1;
	my $pdb_id = $2.$3.$4;
	print "Outputting ".$py_dir."_".$output_count.".py\n";
	my $fhOut = new FileHandle($py_dir."_".$output_count.".py","w");
	print $fhOut "from modeller import *\n";
	print $fhOut "from modeller.automodel import *\n\n";
	print $fhOut "env = environ()\n";
	print $fhOut "env.io.atom_files_directory = \'".$dompdb."\'\n";
	print $fhOut "a = automodel(env, alnfile='".$ali_dir."_".$output_count.".ali',\n";
	print $fhOut "\t\tknowns=('".$pdb_id."'),\n";
	print $fhOut "\t\tsequence=\'".$seq_id."_".$hDomCount->{$seq_id}."\')\n";
	print $fhOut "a.starting_model = 1\n";
	print $fhOut "a.ending_model = 1\n";
  print $fhOut "a.auto_align()\n";
	print $fhOut "a.make()\n";
	print $fhOut "ok_models = filter(lambda x: x[\'failure\'] is None, a.outputs)\n";
	print $fhOut "key = \'molpdf\'\n";
	print $fhOut "ok_models.sort(lambda a,b: cmp(a[key], b[key]))\n";
	print $fhOut "m = ok_models[0]\n";
	print $fhOut "print \"Top model:\%s\" \% m['name']\n";
	$fhOut->close;
}

sub print_ali
{
	my($align_name,$q_start,$q_stop,$s_start,$s_stop,$q_seq,$s_seq,$pdom_ctrl) = @ARG;

	if($pdom_ctrl != 0)
  {
      ($q_start,$q_stop) = get_coords($q_seq,$s_seq);
  }

	($q_seq, $s_seq) = remove_unaligned_ends($q_seq, $s_seq);

	$align_name =~ /^(.+?)\s(.+)/;
	my $uniprot_id = $1;
	my $struct_id = $2;
	my $dom_type = 'CDOM';
	if($pdom_ctrl == 0)
	{
		$dom_type = 'CDOM';
	}
	elsif($pdom_ctrl == 1)
	{
		$dom_type = 'PDOM';
	}
	$struct_id =~ /$dom_type\|(.{4})(.{1})(.{2}):\d+/;
	my $pdb_id = $1.$2.$3;
	my $chain = $2;
	print "Outputting ".$ali_dir."_".$output_count.".ali\n";
	my $fhOut = new FileHandle($ali_dir."_".$output_count.".ali","w");
	if(exists $hDomCount->{$uniprot_id})
	{
		$hDomCount->{$uniprot_id}++;
	}
	else
	{
		$hDomCount->{$uniprot_id} = 1;
	}

	print $fhOut ">P1;".$uniprot_id."_".$hDomCount->{$uniprot_id}."\n";
	if($pdom_ctrl == 0)
	{
		print $fhOut "sequence:::::::::\n";
		print $fhLook $uniprot_id."_".$hDomCount->{$uniprot_id}." ".$q_start." ".$q_stop."\n";
	}
	else
	{
		print $fhOut "sequence:::::::::\n";

#		($q_start,$q_stop) = get_coords($q_seq,$s_seq);
		print $fhLook $uniprot_id."_".$hDomCount->{$uniprot_id}." ".$q_start." ".$q_stop."\n";

	}
	print $fhOut $q_seq."*\n";
	print $fhOut ">P1;".$pdb_id."\n";
	if($pdom_ctrl == 0)
	{
		print $fhOut "StructureX:".$pdb_id."::".$chain."::".$chain."::::\n";
	}
	else
	{
		print $fhOut "StructureX:".$pdb_id."::".$chain."::".$chain."::::\n";
	}
	print $fhOut $s_seq."*\n";
	$fhOut->close;
}

sub get_coords
{
	my($q_seq,$s_seq) = @ARG;
	my $start = 0;
	my $stop = 0;

	my $aSseq = [];
	@$aSseq = split //, $s_seq;
	my $aQseq = [];
	@$aQseq = split //, $q_seq;

	my $res_count = 0;
	my $align_pos = 0;
	foreach my $res (@$aQseq)
	{
		if(@$aQseq[$align_pos] !~ /-/)
		{
			$res_count++;
		}
		if($start == 0)
		{
			if(@$aQseq[$align_pos] !~ /-/ && @$aSseq[$align_pos] !~ /-/)
			{
				$start = $res_count;
			}
		}
		if(@$aQseq[$align_pos] !~ /-/ && @$aSseq[$align_pos] !~ /-/)
		{
			$stop = $res_count;
		}
		$align_pos++;
	}

	return($start,$stop);
}

sub read_ssf
{
	my $fhIn = new FileHandle($ssf_file, "r");
	my $hData = {};
	while(my $line = $fhIn->getline)
	{
		chomp $line;
		my $aEntries = [];
		@$aEntries = split /\s+/, $line;
		my $align_name = @$aEntries[0]." ".@$aEntries[1];
		$hData->{$align_name}{UNIPROTID} = @$aEntries[0];
		$hData->{$align_name}{ALIGNMENTID} = @$aEntries[1];
		$hData->{$align_name}{EVAL} = @$aEntries[10];
		@$aEntries[14] =~ /(\d+):(\d+)/;
		$hData->{$align_name}{START} = $1;
		$hData->{$align_name}{STOP} = $2;
	}
	return($hData);
}


sub remove_unaligned_ends
{
	my($q_seq, $s_seq) = @ARG;

	$q_seq =~ s/\s//g;
	$s_seq =~ s/\s//g;

	#print $q_seq."\n";
	#print $s_seq."\n";
	($q_seq, $s_seq) = remove_leading($q_seq, $s_seq);

	$q_seq = reverse $q_seq;
	$s_seq = reverse $s_seq;
	#print $q_seq."\n";
	#print $s_seq."\n";

	($q_seq, $s_seq) = remove_leading($q_seq, $s_seq);
	$q_seq = reverse $q_seq;
	$s_seq = reverse $s_seq;
	#print $q_seq."\n";
	#print $s_seq."\n";

	return($q_seq, $s_seq);
}

sub remove_leading
{
	my($q_seq, $s_seq) = @ARG;
  print("Q".$q_seq."\n");
  print("S".$a_seq."\n");

	my $align_length = length $s_seq;
	my $s_res = [];
	@$s_res = split //, $s_seq;
	my $q_res = [];
	@$q_res = split //, $q_seq;

	my $res_count = 0;

	foreach my $res (@$s_res)
	{
		if($res !~ /-/ && @$q_res[$res_count] !~ /-/)
		{

			last;
		}
		$res_count++;
	}

	my $new_q = '';
	my $new_s = '';
	if($res_count > 0 && $res_count != ($align_length-1))
	{
		 $new_q  = substr $q_seq, $res_count;
		 $new_s  = substr $s_seq, $res_count;
		 return($new_q, $new_s);
	}
	else
	{
		return($q_seq, $s_seq);
	}

}
