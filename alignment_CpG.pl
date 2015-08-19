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
getopts("hkf:", \%opts);
my $usage = "perl CpG_cal.pl -f <fasta alignment file> [-k]"; 
die "$usage" if $opts{h};
my $keepfile;
$opts{k}?$keepfile=1:$keepfile=0;
my $fasta = $opts{f} or die "$usage";

my $prank_hash_ref = {"d" => $fasta,
			"showanc" => 1,
			"showevents" => 1,
			"keep" => 1,
			"uselogs" => 1,
			"o" => "$fasta.CpG.out",
		};
Func::CpG_with_prank($prank_hash_ref, $keepfile);
exit;
