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
getopts("hkm:f:", \%opts);
my $usage = "perl CpG_cal.pl -f <fasta alignment file> [-m <c|a>] [-k]"; 
die "$usage" if $opts{h};
my $fasta = $opts{f} or die "$usage";

$opts{m} ||= "c"; #set default model as consensus
die "$usage" unless $opts{m} eq "c" || $opts{m} eq "a";
my $model = $pots{m};

if ($model eq 'c'){
	my $str = Bio::AlignIO->new(-file => $fasta,
					-format => 'fasta',
				);
	my $msa = $str->next_aln();
	my $consensus = $msa->consensus_string();
	my $consensus_obj = Bio::LocatableSeq->new(-seq => $consensus,
						-id => "consensus",
						-alphabet => "dna",
					);
	my $copies = $msa->num_sequences();
	my $length = $msa->length();
	my $numCpG = 0; #number of all CpG sites in repSeq
	my $mutTpG = 0; #number of mutated CpG to TpG
	my $mutCpA = 0; #number of mutated CpG to CpA
	my $numC = 0; #number of all C sites and G sites in repSeq (excluding CpG sites)
	my $mutC = 0; #number of all C to T transitions in repSeq (excluding CpG sites)
	my $numG = 0; #number of all G sites in repseq
	my $mutG = 0; #number of all G to A transitions in repseq
	foreach my $genoSeq_obj ($msa->each_seq){
		my ($C, $T, $G, $A, $CpG, $TpG, $CpA) = Func::CpG_rate($consensus_obj, $genoSeq_obj);
		$numC += $C;
		$mutC += $T;
		$numG += $G;
		$mutG += $A;
		$numCpG += $CpG;
		$mutTpG += $TpG;
		$mutCpA += $CpA;
	}
	print "$fasta\t$copies\t$numC\t$mutC\t$numG\t$mutG\t$numCpG\t$mutTpG\t$mutCpA\t$length\n";
	

}
if($model eq 'a'){
	my $keepfile;
	$opts{k}?$keepfile=1:$keepfile=0;
	my $prank_hash_ref = {"d" => $fasta,
				"showanc" => 1,
				"showevents" => 1,
				"keep" => 1,
				"uselogs" => 1,
				"o" => "$fasta.CpG.out",
			};
	Func::CpG_with_prank($prank_hash_ref, $keepfile);
}
exit;
