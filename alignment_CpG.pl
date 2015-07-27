###############################################

=head1 description
Author : Xiaoyu Zhuo
Purpose: Regenerate phylogenetic tree and all ancestor sequences for an alignment file using prank, and calculate CpG/nonCpG substituion ratio.

=head2 USAGE

perl alignment_CpG.pl -f <fasta file> ［-h help]

options: 
-f the fasta alignment file

-h
help

=cut

use FindBin qw( $RealBin);
use lib $RealBin;
use strict;
use warnings;
use Func;
use Bio::SimpleAlign;
use Bio::SeqIO;
use Bio::AlignIO;
use Bio::TreeIO;
use IO::String;
use Data::Dumper;
use Getopt::Std;
my %opts=();
getopts("hf:", \%opts);
my $usage = "perl CpG_cal.pl -f <fasta alignment file>"; 
die "$usage" if $opts{h};

my $fasta = $opts{f} or die "$usage";
my $prankpath =  "~/bioinfo/prank/bin/prank";
#my $prankpath = "~/bioinfo/prank-msa/src/prank";

my $prank_hash_ref = {"d" => $fasta,
			"showanc" => 1,
			"showevents" => 1,
			"keep" => 1,
			"o" => "$fasta.CpG.out",
		};
CpG_with_prank($prank_hash_ref);
exit;


sub CpG_with_prank{
	my $prank_hash_ref = shift or die $!; #prank parameters
	my $sub_root_id = shift;
	#$sub_root_id ||= 0;
	run_prank($prank_hash_ref);
	my $fasta;
	my $events;
	my $tree_file;
	my $out;
	my $domain = $prank_hash_ref->{'o'};
	if ($prank_hash_ref->{'keep'}){
		$fasta = $domain.".anc.fas";
		$events = $domain.".events";
		$tree_file = "$domain.anc.dnd";
		open $out, ">$domain".".CpG.out" or die "$!";
	}
	else{
		$fasta = $domain.".best.anc.fas";
		$events = $domain.".best.events";
		$tree_file = $domain.".best.anc.dnd";
		open $out, ">$domain".".best.CpG.out" or die "$!";
	}
	my $str = Bio::AlignIO->new(-file => $fasta,
					-format => 'fasta',
				);
	my $aln = $str->next_aln();
	my $tree_string; #newick format tree string
	my $tree;

	#get $tree_string and $tree from the dnd file.
	open my $Tree_fh, "<$tree_file" or die "$!";
	$tree_string = <$Tree_fh>;
	my $io = IO::String->new($tree_string);
	my $treeio = Bio::TreeIO->new(-fh => $io,
					-format => 'newick',
				);
	$tree = $treeio->next_tree;
	close ($Tree_fh);
	#close tree filehandle.

	#set up the internal node to calculate CpG ratio within all descendent of it.
	my $sub_root;
	my @sub_tree_nodes;
	my @sub_ids;
	if ($sub_root_id){    #if sub_root_id is provided, do the following things:
		$sub_root = $tree->find_node($sub_root_id);
		@sub_tree_nodes = $sub_root->get_all_Descendents;
		foreach (@sub_tree_nodes) {
			push @sub_ids, $_->id;
		}
	}

	my $sub_mutCpA = 0;
	my $sub_mutTpG = 0;
	my $sub_numCpG = 0;
	my $sub_mutC = 0;
	my $sub_mutG = 0;
	my $sub_numC = 0;
	my $sub_numG = 0;

	my $all_mutCpA = 0;
	my $all_mutTpG = 0;
	my $all_numCpG = 0;
	my $all_mutC = 0;
	my $all_mutG = 0;
	my $all_numC = 0;
	my $all_numG = 0;
	
	open my $Fh, "<$events" or die "$!";
	while (<$Fh>) {
		chomp;
		next unless $_;
		#read tree from dnd file instead.
		#if(/;$/){
		#	$tree_string = $_;
		#	my $io = IO::String->new($tree_string);
		#	my $treeio = Bio::TreeIO->new(-fh => $io,
		#					-format => 'newick',
		#					);
		#	$tree = $treeio->next_tree;
		#}
		my @line = split();
		if($line[0] eq "branch"){ #if this is branch name line
			my $current_node_id = $line[1]; 
			my $node = $tree->find_node($current_node_id);
			next unless $node->ancestor;
			my $ancestor_id = $node->ancestor->id; #find the node leads to current branch
			print "$ancestor_id\t$current_node_id\n";
			my $anc_obj = $aln->get_seq_by_id($ancestor_id);
			my $cur_obj = $aln->get_seq_by_id($current_node_id);
			my ($numC, $mutC, $numG, $mutG, $numCpG, $mutTpG, $mutCpA) = Func::CpG_rate($anc_obj,$cur_obj);
			$all_numC += $numC;
			$all_mutC += $mutC;
			$all_numG += $numG;
			$all_mutG += $mutG;
			$all_numCpG += $numCpG;
			$all_mutTpG += $mutTpG;
			$all_mutCpA += $mutCpA;
			if($sub_root_id && $current_node_id ~~ @sub_ids){ #if sub_root_id is provided, do the following things:
				$sub_numC += $numC;
				$sub_mutC += $mutC;
				$sub_numG += $numG;
				$sub_mutG += $mutG;
				$sub_numCpG += $numCpG;
				$sub_mutTpG += $mutTpG;
				$sub_mutCpA += $mutCpA;
			}

			my $branchratio;
			eval { $branchratio = ($mutTpG/$numCpG)/($mutC/$numC);};
			$branchratio = $@?"NA":sprintf("%.2f",$branchratio); #leave only 2 digits.
			$tree_string =~ s/$current_node_id/$current_node_id $branchratio/;
			#$node->description($branchratio);
			print $out "$current_node_id\t$numC\t$mutC\t$numG\t$mutG\t$numCpG\t$mutTpG\t$mutCpA\n";
		}
	}
	#my $treeout = Bio::TreeIO->new(-file => ">$domain.CpG.new",
	#				-format => 'newick',
	#			);
	#$treeout->write_tree($tree);
	#print Dumper($tree);
	print $out "$tree_string\n";
	print $out "$all_numC\t$all_mutC\t$all_numG\t$all_mutG\t$all_numCpG\t$all_mutTpG\t$all_mutCpA\n";
	my $ratio;
	eval { $ratio = ($all_mutTpG/$all_numCpG)/($all_mutC/$all_numC);};
	$ratio = "NA" if $@;
	printf $out "%.2f", $ratio;
	print $out "\n";

	if($sub_root_id){  #if sub_root_id is provided, do the following things:
		my $subratio;
		eval {$subratio = ($sub_mutTpG/$sub_numCpG)/($sub_mutC/$sub_numC);};
		$subratio = $@?"NA":sprintf("%.2f",$subratio);
		print $out "$sub_root_id\n$sub_numC\t$sub_mutC\t$sub_numG\t$sub_mutG\t$sub_numCpG\t$sub_mutTpG\t$sub_mutCpA\n";
		print $out "$sub_root_id subratio is $subratio\n";
	}
	close ($Fh);
	close ($out);
}


#how to build the prank_hash: "d" => "some.fas" for -d=some.fas, "F" =>0 for -F in command line.
sub run_prank{ # argument is the prank argument hash ref.
	my $prank = $prankpath;
	my %inhash = %{shift()};
	my @params = ($prank);
	foreach my $key (keys %inhash){
		#	if ($inhash{$key} == 1){
		#	$para = "-".$key;
		#}
		#else{
		#		$para = "-".$key."=".$inhash{$key};
		#}
		next unless $inhash{$key};
		my $param;
		{
			no warnings; #disable warnings for this block (string in numeric ===)
			$param = ($inhash{$key} == 1)?("-".$key):("-".$key."=".$inhash{$key});
		}
		push @params, $param;
	}
	my $params = join ' ', @params;
	print "$params\n";
	system($params);
}



