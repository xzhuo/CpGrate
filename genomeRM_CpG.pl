##################################################

=head1 description

Author  : Xiaoyu Zhuo
Purpose : Read a repeatmasker alignment output and regenerate the consensus sequence and then calculate CpG/nonCpG substitution ratio. Note: insertions are ignoared during regeneration. Directly modified from trueCONSENSUS
modified from John Pace's script.

=head2 USAGE

perl genomeRM_CpG.pl -a <RM align file> [-m <c|a>] [-p <con|CG|aln|all>] [-c <class name>] [-r <repeat name>] [-s <INT> -l <INT>] [-b <INT>] [-d <replibrary name>] [-u <INT>] [-h help]

options:
-a the repeatmasker alignment file.

-m the model, "c" for consensus based CpG calculation, "a" for alignment based calculation. Default value is "c".

-p 
con: output only regenerated consensus, default. Mutually exclusive with "aln".
CG: output only CpG/nonCpG calculate output.
aln: output whole alignment instead of consensus, Mutually exclusive with "con".
all: output both regenerated consensus and CpG/nonCpG calculation.

-c
target class. with this option, only selected class is calculated, SINE|LINE|LTR, etc.

-r
target class. Only selected repeat is calculated. LTR5B, etc.

-s
maxstart. int. Only repeats with repstart <= maxstart are used for consensus regeneration always come with -l. Recommend 50.

-l
maxleft. int. Only repeats with repleft <= maxleft are used for consensus regeneration. Work together with -s to defined complete TE. Recommend 50.

-b
minlength. int. only repeats with sequence length >= minlength are used for consensus regeneration. Recommend 100.

-d
repeat library used for repeatmasker. Highly recommended to verify sequence alignment information.

-u
number of processes in multithreads.

ps: if -s, -l and -b are not defined, any fragment matches particular TE would be included in the alignment for calculation.

-h
help

-i
Include insertions in the alignment or not. If it is false, the MSA length would be exactly the consensus length. if it is true, then insertions of individual copy would be inclued in the MSA. Default option is false. it is only compatible with "-p aln".

=cut
#use forks;
#use forks::shared;
use Parallel::ForkManager;
use FindBin qw($RealBin);
use lib $RealBin;
use warnings;
use strict;
use Func;
use Getopt::Std;
use Bio::SeqIO;
use Bio::SimpleAlign;
use Bio::LocatableSeq;
use Bio::AlignIO;
use Bio::Align::Utilities qw(cat);
use Bio::DB::Fasta;
use File::Temp qw/ tempfile/;
use Fcntl qw(:flock SEEK_END O_CREAT O_RDWR);
use IO::Handle;
use Statistics::Basic qw(:all);
use List::Util qw(max);
use BerkeleyDB;
use MLDBM qw(BerkeleyDB::Hash);
# use Set::IntervalTree;
use Data::Dumper;
use File::Basename;

STDERR->autoflush(1);
STDOUT->autoflush(1);

my %opts=();
getopts("hia:m:p:c:r:s:l:b:d:u:", \%opts);
my $usage = "perl ParseRMalign_regenerate_consensus.pl -a <RM align file> [-m <c|a>] [-p <con|CG|aln|all>] [-c <class name>] [-r <repeat name>] [-s <INT> -l <INT>] [-b <INT>] [-d <replibrary name>] [-u <INT>] [-i]";
die "$usage" if $opts{h};
my $fileName = $opts{a} or die "$usage";
# my $RMfile = $opts{o} or die "$usage";
$opts{p} ||= "con"; #set default parameter for p: only output consensus.
die "$usage" unless $opts{p} eq "con" || $opts{p} eq "CG" || $opts{p} eq "all" || $opts{p} eq "aln";
$opts{m} ||= "c"; #set default model as consensus
die "$usage" unless $opts{m} eq "c" || $opts{m} eq "a";
die "$usage" if $opts{i} and ($opts{p} eq "con" || $opts{p} eq "CG" || $opts{p} eq "all");
$opts{u} ||= 1;
my $cpus = $opts{u};
my $purpose = $opts{p};
my $model = $opts{m};
my $targetClass = $opts{c} if defined $opts{c};
my $targetRepeat= $opts{r} if defined $opts{r};
my $maxStart;
my $maxLeft;
if (defined $opts{s} && defined $opts{l}){
	$maxStart = $opts{s};
	$maxLeft = $opts{l};
}
elsif (defined $opts{s} || defined $opts{l}){
	die "$usage";
}
my $minlength = $opts{b};
my @alignArray = ();

my $db;
my $replib;
my $tempfasta;
if (defined $opts{d}){
	print "reading and loading repeat library\n";
	print "$opts{d}\n";
	#convert embl file to fasta file
	my $replib = $opts{d};
	$tempfasta = substr($replib,0,-4)."fa";
	print "$tempfasta\n";
	my ($fh,$tmpembl) = tempfile();
	#somehow the tempfile does not get deleted.
	open my $oriFh, "<$replib";
	open $fh, ">$tmpembl";
	while (<$oriFh>){
		if(substr($_,0,2) eq "ID"){
			print $fh "$_";
			last;
		}
	}
	while (<$oriFh>){
		print $fh "$_";
	}
	close($oriFh);
	my $embl = Bio::SeqIO->new(-file => $tmpembl,
					-format => 'EMBL',
				);
	my $fasta = Bio::SeqIO->new(-file=> ">$tempfasta",
					-format => 'Fasta',
				);
	while (my $seq = $embl->next_seq()){
		my $emblseqstr = $seq->seq();
		$emblseqstr =~ s/[^ATGC]/-/g; # non-GATC NTs are not correctly revcom in repeatmasker alignment.
		$seq->seq($emblseqstr);
		$fasta->write_seq($seq);
	}
	close($fh);
	unlink $tmpembl; #delete the file manually
	#build db::fasta
	$db = Bio::DB::Fasta->new($tempfasta);
	print "tempfasta: $tempfasta\n"; #for debug
	print "repeat library loaded!\n";
}

# print "reading and processing ".$RMfile."\n";
# open the RM output file
# my %trees = ();
# open my $RMout, "<$RMfile";
# while(<$RMout>){
#	# get the columns:blabla
#	if(exists $trees{$chr}){
# 		$trees{$chr}->insert($te, $chr_start, $chr_end);
# 	}
# 	else{
# 		$trees{$chr} = Set::IntervalTree->new;
# 		$trees{$chr}->insert($te, $chr_start, $chr_end);
# 	}
# }
# close($RMout);

# print "RMout entry in IntervalTree ready for query!\n";

print "Reading and processing ".$fileName."\n";

#Open the align file
open my $alignfile, "<$fileName";
my $outName;
my $outFileName = basename($fileName);
my $outFasta;
if (defined $targetClass || defined $targetRepeat || defined $maxStart || defined $maxLeft || defined $minlength){
	my $target = "";
	$target = $target.$targetClass if defined $targetClass;
	$target = $target.$targetRepeat if defined $targetRepeat;
	$target = $target.".$maxStart-$maxLeft" if defined $maxStart && defined $maxLeft;
	$target = $target.".$minlength" if defined $minlength;
	$outName = "$outFileName.$target.CpG.out";
	$outFasta = "$outFileName.$target.new_consensus.fa";
#	$outName = "fAlb15.$target.CpG.out";
#	$outFasta = "fAlb15.$target.new_consensus.fa";
}
else{
	$outName = "$outFileName.CpG.out";
	$outFasta = "$outFileName.new_consensus.fa";
#	$outName = "fAlb15.CpG.out";
#	$outFasta = "fAlb15.new_consensus.fa";
}

#Put the lines from the file into an array, but only the ones that aren't just returns or the tv/ts line
#Also split the lines appropriately so they can be handled later
while(<$alignfile>){
	if(($_ ne "\n") && (!($_ =~ /^   /))){ #only check for the lines that are not just returns orthe tv/ts line
		my @tempArray = (); #temporary array to hold the token from the line that is split
		@tempArray = split(); #split the line and assign it to @tempArray
		my $temp = ''; #temporary variable used for concatenation of tempArray rows
		for(my $i=0;$i<=$#tempArray;$i++){ #loop through tempArray
			$temp=$temp.$tempArray[$i]."\t"; #concatenate tempArray rows
		}
		push(@alignArray, $temp); #push $temp into alignArray
	}
}
#close the alignfile
close($alignfile);
print "\nThe align file has been loaded\n";

#open the outfile
my $outFile;
#my @outFile  : shared;
my $newFasta;
#my @newFasta : shared;

#parse alnHash to regenerate consensus for all repeats:


if ($purpose eq "CG" || $purpose eq "all"){
	open $outFile, ">$outName";
	print $outFile "te\tcopies\tnumC\tmutC\tnumG\tmutG\tnumCpG\tmutTpG\tmutCpA\tlength\n";
#	foreach (@outFile){
#		print $outFile "$_\n";
#	}
#	close ($outFile);
}
if ($purpose eq "con" || $purpose eq "aln" || $purpose eq "all"){
	open $newFasta, ">$outFasta";
#	foreach (@newFasta){
#		print $newFasta ">$_->[0]\n$_->[1]\n";
#	}
#	close ($newFasta);
}


#a HASH holding each repeat with its MSA
for my $dbfile (qw(mldbm)){
	unlink $dbfile;
	die "unable to delete dbfile!\n\n" if -e $dbfile;
}

my %repHash = (); #This is simple perl hash access.

#Loop through the array, parsing the lines and add to an output file
for(my $i = 0; $i<= $#alignArray;$i++){
#	print "$i\n" if (($i % 10000) == 0); #print line number of every 10000 lines

	#Declare the variables
	my $genoName = '';
	my $genoStart = 0;
	my $genoEnd = 0;
	my $repName = '';
	my $repClass = '';
	my $repFamily = '';
	my $strand = '+';
	my $repStart = 0;
	my $repEnd = 0;
	my $repLeft = 0;
	my $chrSeq = '';
	my $repSeq = '';
	my $beginningRow = 0; #beginningRow is used to keep track of the start of the repeat seq block
	my $endingRow = 0;  #endingRow is used to keep track of the end of the repeat seq block
	my $id = 0;

	#Check to see if the first two characters of the array row are digits (indicate this is header line of each record)
	if($alignArray[$i] =~ /^\d+\s/){
		$beginningRow = $i; #Set beginningRow

		#Find the subsequent row that begins with 'Gap_init'
		for(my $j = ($i+1); $j<= $#alignArray;$j++){
			if($alignArray[$j] =~ /^Gap_init/){
				$endingRow = $j;  #Set the endingRow
				last;  #End the for loop
			}
		}
		#At this point, we have the beginning and ending rows
		#Check to see if the repeat is in + or - (C) orientation
		#Split the array row and see if the 9th row is a C
		my @arraySplit = ();
		my $headline = $alignArray[$i];
		$headline =~ s/#/\t/;
		@arraySplit = split(/\t/, $headline);
		if($arraySplit[8] eq 'C'){	  #This means the repeat is in - (C) orientation
			#Set the variables
			$genoName = $arraySplit[4];
			$genoStart = $arraySplit[5];
			$genoEnd = $arraySplit[6];
			$strand = $arraySplit[8];
			$repName = $arraySplit[9];
			$repClass = $arraySplit[10];
			($repLeft) = $arraySplit[11] =~ /\((.+)\)/;
			$repEnd = $arraySplit[12];
			$repStart = $arraySplit[13];
			$id = $arraySplit[-1];
			#repClass has to be split further (modified from John Pace's)
			($repClass, $repFamily) = split(/\//, $repClass);
		}

		else { #The repeat is in + orientation
			#Set the variables
			$genoName = $arraySplit[4];
			$genoStart = $arraySplit[5];
			$genoEnd = $arraySplit[6];
			$strand = '+';
			$id = $arraySplit[-1];
			if ($arraySplit[8] eq '+'){
				$repName = $arraySplit[9]; #note: the mm10 and hg19 are different at these lines!
				$repClass = $arraySplit[10];
				$repStart = $arraySplit[11];
				$repEnd = $arraySplit[12];
				($repLeft) = $arraySplit[13] =~ /\((.+)\)/;
			}
			else{
				$repName = $arraySplit[8]; #note: the mm10 and hg19 are different at these lines!
				$repClass = $arraySplit[9];
				$repStart = $arraySplit[10];
				$repEnd = $arraySplit[11];
				($repLeft) = $arraySplit[12] =~ /\((.+)\)/;
			}
			#repClass has to be split further
			($repClass, $repFamily) = split(/\//, $repClass);
		}

		next if $repLeft < 0; #some problem with repeatMasker program? skip this suspicious record if $repLeft <0.
		#example in hg19: 5114  22.26 8.60 4.26  chr3  165270866 165271015 (32751415) + MLT1-int        LTR/ERVL-MaLR       1634   1734 (-363)  410

		next if ((defined $targetClass) && ($repClass ne $targetClass)) || ((defined $targetRepeat) && ($repName ne $targetRepeat)); #only do this with the $targetClass or $targetRepeat.

		#All the variables are set now, except chrSeq and repSeq.
		for (my $k = ($beginningRow +1);$k<=($endingRow-1);$k++){
			#Split alignArray[$k] to see if it is a sequence line.
			#First check for + strand
			@arraySplit = ();  #Clear arraySplit
			@arraySplit = split(/\t/, $alignArray[$k]);
			if (($arraySplit[0] ne 'Matrix') &&($arraySplit[0] ne 'Transitions') && ($arraySplit[0] ne 'Gap_init') && ($arraySplit[0] ne 'Kimura') && ($strand eq '+')){
				#Check to see if this is a chr sequence or a rep sequence
				if ( $arraySplit[0] eq substr($genoName, 0, length($arraySplit[0])) || $arraySplit[0] =~ /^\Q$genoName\E/){	#This is a chr seq
					$chrSeq = $chrSeq.$arraySplit[2];
				}
				else{  #This is the rep seq
					$repSeq = $repSeq.$arraySplit[2];
				}
			}

			#Next check for the - (C) strand
			if (($arraySplit[0] ne 'Matrix') && ($arraySplit[0] ne 'Transitions') && ($arraySplit[0] ne 'Gap_init') && ($arraySplit[0] ne 'Kimura') && ($strand eq 'C')){
				#Check to see if this is a chr sequence or a rep sequence
				if ($arraySplit[0] eq substr($genoName, 0, length($arraySplit[0])) || $arraySplit[0] =~ /^\Q$genoName\E/){	#This is a chr seq
					$chrSeq = $chrSeq.$arraySplit[2];
				}
				else{  #This is a rep seq
					$repSeq = $repSeq.$arraySplit[3];
				}
			}
		}
		$chrSeq = Func::revcom($chrSeq) if $strand eq 'C';
		$repSeq = Func::revcom($repSeq) if $strand eq 'C';
		#delete leading gaps and trailing gaps in repSeq, and delete correlated nt in chrSeq
		$repSeq =~ /^(-*)[^-].*[^-](-*)$/;
		my $leadinggap_length = length($1);
		my $trailinggap_length = length($2);
		if ($trailinggap_length > 0){
			$repSeq = substr($repSeq,$leadinggap_length,-$trailinggap_length);
			$chrSeq = substr($chrSeq,$leadinggap_length,-$trailinggap_length);
		}
		else{
			$repSeq = substr($repSeq,$leadinggap_length);
			$chrSeq = substr($chrSeq,$leadinggap_length);
		}
		#delete insertions and create new mutiple alignment
		$chrSeq = uc($chrSeq);
		$repSeq = uc($repSeq);
		my $gaps_ref = Func::find($repSeq, "-"); #gap in repSeq means insertions in chrSeq, we have to find them out. save them to another hash if $opts{i}, else delete them.
		my @descending_gaps = sort { $b<=>$a } @$gaps_ref;
		my @indel_matrix = ();  # an array of hashes with each item like: {"pos"=>3, "seq"=> 'TA'}
		foreach my $gap(@descending_gaps){
			my $indel_seq = substr($chrSeq, $gap-1, 1);
			if (scalar(@indel_matrix)){
				for my $indel (@indel_matrix){
					$indel->{'pos'} -= 1;
				}
				if ($indel_matrix[-1]{'pos'} == $gap){
					$indel_matrix[-1]{'seq'} = $indel_seq.$indel_matrix[-1]{'seq'};
				}
				else{
					push @indel_matrix, {'pos' => $gap, 'seq' => $indel_seq};
				}
			}
			else{
				push @indel_matrix, {'pos' => $gap, 'seq' => $indel_seq};
			}
			substr($chrSeq, $gap-1, 1) = "";
			substr($repSeq, $gap-1, 1) = "";
		}
		@indel_matrix = sort {$a->{'pos'}<=>$b->{'pos'}} @indel_matrix;  # sort indel matrix by position ascending.
		$chrSeq =~ s/[^ATGC]/-/g;
		$repSeq =~ s/[^ATGC]/-/g; # non-GATC NTs are not correctly revcom in repeatmasker alignment.

		next if $chrSeq eq ""; #remove suspicious sequences (I don't know why they are in the alignment in the first place. but the fact that they have different seq length causes a problem for MSA building).
		#next if Func::pairwise_identity($chrSeq,$repSeq) > 40;

		my $seq_hash_ref = { "chrSeq" => $chrSeq,
					"repSeq" => $repSeq,
					"chr" => $genoName,
					"chrStart" => $genoStart,
					"chrEnd" => $genoEnd,
					"repStart" => $repStart,
					"repEnd" => $repEnd,
					"repLeft" => $repLeft,
					"strand" => $strand,
					"id" => $id,
					"indel_matrix" => \@indel_matrix,
				};
		if(exists $repHash{$repName}){
			#with perl hash:
			push @{$repHash{$repName}}, $seq_hash_ref;
		}
		else{
			$repHash{$repName} = [$seq_hash_ref];
		}
	}
}
undef @alignArray;
print "everything te in seq_hash_ref now,building repdb:\n";
#use mldb:
tie my %repDB, 'MLDBM', -Filename => 'mldbm', 
				-Flags => DB_CREATE 
				or die $!;

foreach my $key(keys %repHash){
	$repDB{$key} = $repHash{$key};
}
undef %repHash;
print "now repDB is done!\n";
my $pm = Parallel::ForkManager->new($cpus);
while (my ($te, $array_ref) = each(%repDB)){
	my $pid = $pm->start and next;
	print "te name is $te\n";
	print "working on $te:\n";
	@$array_ref = sort {$a->{"chr"} cmp $b->{"chr"} or $a->{"chrStart"} <=> $b->{"chrStart"} or $a->{"chrEnd"} <=> $b->{"chrEnd"}} @$array_ref; 
	my $curr_ref;
	my $repfull;
	my $msa;
	my %msa_indels;
	if (defined $opts{d}){
		$repfull = $db->get_Seq_by_id($te);
		my $repfullseq = $repfull->seq();
	}
	foreach my $hash_ref(@$array_ref){
		
		#next correct repStart and repEnd if they are wrong:
		if (defined $opts{d}){ 
			my $seqaln = $hash_ref->{"repSeq"};
			my $tempStart = $hash_ref->{"repStart"};
			my $tempEnd = $hash_ref->{"repEnd"};
			my $tempLeft = $hash_ref->{"repLeft"};
			my $seqpos = $repfull->subseq($tempStart,$tempEnd);
			unless ($seqaln eq $seqpos){
				my $tmpCheck = 0;
				for my $diff(-1,1,2,-2,-3,3,-4,4,-5,5,-6,6,-7,7,-8,8,-9,9,-10,10){
					$seqpos = $repfull->subseq($tempStart + $diff,$tempEnd + $diff);
					#print "$diff\n$seqpos\n";
					if ($seqaln eq $seqpos){
						$hash_ref->{"repStart"} = $tempStart + $diff;
						$hash_ref->{"repEnd"} = $tempEnd + $diff;
						$hash_ref->{"repLeft"} = $tempLeft - $diff;
						$tmpCheck = 1;	#swtich it!
						last;
					}
				}
				die "wrong alignment somehow!\n$hash_ref->{'chr'}\t$hash_ref->{'chrStart'}\t$hash_ref->{'chrEnd'}\n$seqaln\n$seqpos\n\n" unless $tmpCheck;
			}
		}

		unless (defined $curr_ref){
			$curr_ref = $hash_ref;
		}

		elsif ($hash_ref->{"id"} == $curr_ref->{"id"} and $hash_ref->{"chr"} eq $curr_ref->{"chr"} and $hash_ref->{"strand"} eq $curr_ref->{"strand"} and $hash_ref->{"chrStart"} < $curr_ref->{"chrEnd"} + 10000 and $curr_ref->{"strand"} eq "+"?abs($hash_ref->{"repStart"}-$curr_ref->{"repEnd"})<20:abs($hash_ref->{"repEnd"}-$curr_ref->{"repStart"})<20){
			#merge record
			if ($curr_ref->{"strand"} eq "+"){
				if($hash_ref->{"repStart"}>$curr_ref->{"repEnd"}){
					$curr_ref->{"chrSeq"} = $curr_ref->{"chrSeq"}."-"x($hash_ref->{"repStart"}-$curr_ref->{"repEnd"}-1).$hash_ref->{"chrSeq"};
					$curr_ref->{"repSeq"} = $curr_ref->{"repSeq"}."-"x($hash_ref->{"repStart"}-$curr_ref->{"repEnd"}-1).$hash_ref->{"repSeq"};
					for my $indel (@{$hash_ref->{"indel_matrix"}}){
						push @{$curr_ref->{"indel_matrix"}}, {"pos" => $indel->{"pos"}+$hash_ref->{"repStart"}-$curr_ref->{"repStart"}, "seq" => $indel->{"seq"}};
					}

					$curr_ref->{"repEnd"} = $hash_ref->{"repEnd"};
					$curr_ref->{"repLeft"} = $hash_ref->{"repLeft"};
					$curr_ref->{"chrEnd"} = $hash_ref->{"chrEnd"};
				}
				else{
					#deal with overlap
					my $hash_overlap = $hash_ref->{"repEnd"}-$hash_ref->{"repStart"}+1;
					my $curr_overlap = $curr_ref->{"repEnd"}-$hash_ref->{"repStart"}+1;
					my $min_overlap = $hash_overlap<$curr_overlap?$hash_overlap:$curr_overlap;
					my $curr_chrSeq = substr($curr_ref->{"chrSeq"},0,0-$curr_overlap);
					my $curr_repSeq = substr($curr_ref->{"repSeq"},0,0-$curr_overlap);
					my $hash_chrSeq;
					my $hash_repSeq;
					if ($min_overlap == $curr_overlap){
						$hash_chrSeq = substr($hash_ref->{"chrSeq"},$min_overlap);
						$hash_repSeq = substr($hash_ref->{"repSeq"},$min_overlap);
					}
					else{
						$hash_chrSeq = substr($curr_ref->{"chrSeq"},$min_overlap-$curr_overlap);
						$hash_repSeq = substr($curr_ref->{"repSeq"},$min_overlap-$curr_overlap);
					}
					my $tempseq1 = substr($curr_ref->{"chrSeq"},0-$curr_overlap,$min_overlap);
					my $tempseq2 = substr($hash_ref->{"chrSeq"},0,$min_overlap);
					my $tempOverlap_length1 = 0;
					my $tempOverlap_length2 = 0;
					for my $individual_indel1 (@{$curr_ref->{"indel_matrix"}}){
						if ($individual_indel1->{"pos"} > $hash_ref->{"repStart"}-$curr_ref->{"repStart"} and $individual_indel1->{"pos"} < $hash_ref->{"repStart"}+$min_overlap-$curr_ref->{"repStart"}+1){
							$tempOverlap_length1 += length($individual_indel1->{"seq"});
						}
					}
					for my $individual_indel2 (@{$hash_ref->{"indel_matrix"}}){
						if ($individual_indel2->{"pos"} > 0 and $individual_indel2->{"pos"} < $min_overlap){
							$tempOverlap_length2 += length($individual_indel2->{"seq"});
						}
					}
					if ($tempOverlap_length1 <= $tempOverlap_length2){  # insertions in curr_ref overlap region is shorter, so curr_ref overlap sequence is better.
						if ($min_overlap == $curr_overlap){
							# first 2 section from curr_ref, last 1 section from hash_ref: keep all from curr_ref, push last part of hash_ref.
							for my $indel (@{$hash_ref->{"indel_matrix"}}){
								if ($indel->{"pos"} >= $min_overlap){
									push @{$curr_ref->{"indel_matrix"}}, {"pos" => $indel->{"pos"}+$hash_ref->{"repStart"}-$curr_ref->{"repStart"}, "seq" => $indel->{"seq"}};
								}
							}
						}
						#else all from curr_ref
					}
					else{  # insertions in hash_ref overlap region is shorter, so hasg_ref overlap sequence is better.
						if ($min_overlap == $curr_overlap){
							#first 1 section from curr_ref, last 2 sections from hash_ref: remove last part of curr_ref, push all hash_ref
							for my $indel (@{$curr_ref->{"indel_matrix"}}){
								if ($indel->{"pos"} > $hash_ref->{"repStart"}-$curr_ref->{"repStart"}){
									$indel->{"seq"} = "";
									$indel->{"pos"} = 0;  # is 0 ok here?
								}
							}
							for my $indel (@{$hash_ref->{"indel_matrix"}}){
								push @{$curr_ref->{"indel_matrix"}}, {"pos" => $indel->{"pos"}+$hash_ref->{"repStart"}-$curr_ref->{"repStart"}, "seq" => $indel->{"seq"}}; 
							}

						}
						else{
							#first and 3rd from curr_ref, middle one from hash_ref
							for my $indel (@{$curr_ref->{"indel_matrix"}}){
								if ($indel->{"pos"} > $hash_ref->{"repStart"}-$curr_ref->{"repStart"} and $indel->{"pos"} < $hash_ref->{"repEnd"}-$curr_ref->{"repStart"}+1){
									$indel->{"seq"} = "";
									$indel->{"pos"} = 0;  # is 0 ok here?
								}
							}
							for my $indel (@{$hash_ref->{"indel_matrix"}}){
								push @{$curr_ref->{"indel_matrix"}}, {"pos" => $indel->{"pos"}+$hash_ref->{"repStart"}-$curr_ref->{"repStart"}, "seq" => $indel->{"seq"}}; 
							}
						}
					}
					my $consensus = substr($hash_ref->{"repSeq"},0,$min_overlap);
					my $overseq = "";
					#print "$curr_chrSeq\n$tempseq1\n$tempseq2\n$consensus\n$hash_chrSeq\n\n";
					for(my $i=0;$i<length($consensus);$i++){
						my $nt = substr($tempseq1,$i,1) eq substr($tempseq2,$i,1)?substr($tempseq1,$i,1):substr($consensus,$i,1);
						$overseq=$overseq.$nt;
					}
					$curr_ref->{"chrSeq"} = $curr_chrSeq.$overseq.$hash_chrSeq;
					$curr_ref->{"repSeq"} = $curr_repSeq.$overseq.$hash_repSeq;
					if ($curr_ref->{"repEnd"} <= $hash_ref->{"repEnd"}){
						$curr_ref->{"repEnd"} = $hash_ref->{"repEnd"};
						$curr_ref->{"repLeft"} = $hash_ref->{"repLeft"};
					}
					$curr_ref->{"chrEnd"} = $hash_ref->{"chrEnd"} if $hash_ref->{"chrEnd"} > $curr_ref->{"chrEnd"};
				}
			}
			else{  # - strand
				if($hash_ref->{"repEnd"}<$curr_ref->{"repStart"}){
					$curr_ref->{"chrSeq"} = $hash_ref->{"chrSeq"}."-"x($curr_ref->{"repStart"}-$hash_ref->{"repEnd"}-1).$curr_ref->{"chrSeq"};
					$curr_ref->{"repSeq"} = $hash_ref->{"repSeq"}."-"x($curr_ref->{"repStart"}-$hash_ref->{"repEnd"}-1).$curr_ref->{"repSeq"};
					for my $indel (@{$curr_ref->{"indel_matrix"}}){
						$indel->{"pos"} = $indel->{"pos"} + $curr_ref->{"repStart"} - $hash_ref->{"repStart"};
					}
					push @{$curr_ref->{"indel_matrix"}}, @{$hash_ref->{"indel_matrix"}};
					$curr_ref->{"chrEnd"} = $hash_ref->{"chrEnd"};
					$curr_ref->{"repStart"} = $hash_ref->{"repStart"};
				}
				else{
					#deal with overlap
					my $hash_overlap = $hash_ref->{"repEnd"}-$hash_ref->{"repStart"}+1;
					my $curr_overlap = $hash_ref->{"repEnd"}-$curr_ref->{"repStart"}+1;
					my $min_overlap = $hash_overlap<$curr_overlap?$hash_overlap:$curr_overlap;
					my $curr_chrSeq = substr($curr_ref->{"chrSeq"},$curr_overlap);
					my $curr_repSeq = substr($curr_ref->{"repSeq"},$curr_overlap);
					my $hash_chrSeq;
					my $hash_repSeq;
					if($min_overlap == $curr_overlap){
						$hash_chrSeq = substr($hash_ref->{"chrSeq"},0,0-$curr_overlap);
						$hash_repSeq = substr($hash_ref->{"repSeq"},0,0-$curr_overlap);
					}
					else{
						$hash_chrSeq = substr($curr_ref->{"chrSeq"},0,$curr_overlap-$min_overlap);
						$hash_repSeq = substr($curr_ref->{"repSeq"},0,$curr_overlap-$min_overlap);
					}
					my $tempseq1 = substr($curr_ref->{"chrSeq"},$curr_overlap-$min_overlap,$min_overlap);
					my $tempseq2 = substr($hash_ref->{"chrSeq"},0-$min_overlap);
					my $tempOverlap_length1 = 0;
					my $tempOverlap_length2 = 0;
					for my $individual_indel1 (@{$curr_ref->{"indel_matrix"}}){
						if ($individual_indel1->{"pos"} < $curr_overlap and $individual_indel1->{"pos"} > $curr_overlap-$min_overlap){
							$tempOverlap_length1 += length($individual_indel1->{"seq"});
						}
					}
					for my $individual_indel2 (@{$hash_ref->{"indel_matrix"}}){
						if ($individual_indel2->{"pos"} > $hash_overlap-$min_overlap){
							$tempOverlap_length2 += length($individual_indel2->{"seq"});
						}
					}
					if ($tempOverlap_length1 <= $tempOverlap_length2){
						if ($min_overlap == $curr_overlap){
							# last 2 section from curr_ref, first 1 section from hash_ref: keep all from curr_ref, push first part of hash_ref.
							for my $indel (@{$curr_ref->{"indel_matrix"}}){
								$indel->{"pos"} += $curr_ref->{"repStart"}-$hash_ref->{"repStart"};
							}
							for my $indel (@{$hash_ref->{"indel_matrix"}}){
								if($indel->{"pos"} <= $curr_ref->{"repStart"}-$hash_ref->{"repStart"} and $indel->{"pos"} > 0){
									push @{$curr_ref->{"indel_matrix"}}, $indel;
								}
							}
						}
						#else all from curr_ref
					}
					else{
						if ($min_overlap == $curr_overlap){
							#last 1 section from curr_ref, first 2 sections from hash_ref: remove first part of curr_ref, push all hash_ref
							for my $indel (@{$curr_ref->{"indel_matrix"}}){
								if ($indel->{"pos"} < $curr_overlap){
									$indel->{"pos"} = 0;
									$indel->{"seq"} = "";
								}
								else{
									$indel->{"pos"} += $curr_ref->{"repStart"}-$hash_ref->{"repStart"};
								}
							}
							for my $indel (@{$hash_ref->{"indel_matrix"}}){
								push @{$curr_ref->{"indel_matrix"}}, $indel;
							}
						}
						else{
							#first and 3rd from curr_ref, middle one from hash_ref
							for my $indel (@{$curr_ref->{"indel_matrix"}}){
								if ($indel->{"pos"} < $curr_overlap and $indel->{"pos"} > $hash_ref->{"repStart"}-$curr_ref->{"repStart"}){
									$indel->{"pos"} = 0;
									$indel->{"seq"} = "";
								}
							}
							for my $indel (@{$hash_ref->{"indel_matrix"}}){
								push @{$curr_ref->{"indel_matrix"}}, {"pos"=>$indel->{"pos"}+$hash_ref->{"repStart"}-$curr_ref->{"repStart"}, "seq"=>$indel->{"seq"}};
							}
						}
					}
					my $consensus = substr($hash_ref->{"repSeq"},0-$min_overlap);
					my $overseq = "";
					for(my $i=0;$i<length($consensus);$i++){
						my $nt = substr($tempseq1,$i,1) eq substr($tempseq2,$i,1)?substr($tempseq1,$i,1):substr($consensus,$i,1);
						$overseq=$overseq.$nt;
					}
					$curr_ref->{"chrSeq"} = $hash_chrSeq.$overseq.$curr_chrSeq;
					$curr_ref->{"repSeq"} = $hash_repSeq.$overseq.$curr_repSeq;
					if($curr_ref->{"repStart"} >= $hash_ref->{"repStart"}){
						$curr_ref->{"repStart"} = $hash_ref->{"repStart"};
					}
					$curr_ref->{"chrEnd"} = $hash_ref->{"chrEnd"} if $hash_ref->{"chrEnd"} > $curr_ref->{"chrEnd"};
				}
			}
		}
		else{
			unless (((defined $maxStart) && ($curr_ref->{"repStart"} > $maxStart)) || ((defined $maxLeft) && ($curr_ref->{"repLeft"} > $maxLeft)) || ((defined $minlength) && ($curr_ref->{"chrEnd"} - $curr_ref->{"chrStart"} + 1 < $minlength))){ #filter out maxstart, maxend and minlength parameter
				#print "$curr_ref->{'chr'}\t$curr_ref->{'chrStart'}\t$curr_ref->{'chrEnd'}\t$curr_ref->{'strand'}\nstart is $curr_ref->{'repStart'}\nend is $curr_ref->{'repLeft'}\n\n";
				#assign new bio::seq to alignemnt
				my $chrSeq = "-"x($curr_ref->{"repStart"}-1).$curr_ref->{"chrSeq"}."-"x($curr_ref->{"repLeft"});
#				print "$headline\n$chrSeq\n";
				my $tempId = $curr_ref->{"chr"}."_".$curr_ref->{"chrStart"}."-".$curr_ref->{"chrEnd"}.$curr_ref->{"strand"};
				my $tempSeq = Bio::LocatableSeq->new(-seq => $chrSeq,
									-id => $tempId,
									-alphabet => "dna",
								);
				#add the chrSeq with insertions deleted to the MSA:
				if (defined $msa){
					$msa->add_seq($tempSeq);
					for my $indels ($curr_ref->{"indel_matrix"}){
						for my $single_insertion (@$indels){
							$msa_indels{$single_insertion->{'pos'} + $curr_ref->{"repStart"}-1}{$tempId} = $single_insertion->{'seq'};
						}
					}
				}
				else{
					$msa = Bio::SimpleAlign->new(-seqs => [$tempSeq]);
					for my $indels ($curr_ref->{"indel_matrix"}){
						for my $single_insertion (@$indels){
							$msa_indels{$single_insertion->{'pos'} + $curr_ref->{"repStart"}-1}{$tempId} = $single_insertion->{'seq'};
						}
					}
				}
			}
			else{
				#print "not included:\n$curr_ref->{'chr'}\t$curr_ref->{'chrStart'}\t$curr_ref->{'chrEnd'}\t$curr_ref->{'strand'}\nstart is $curr_ref->{'repStart'}\nend is $curr_ref->{'repLeft'}\n\n";
			}
			$curr_ref = $hash_ref;
		}
	}
	unless (((defined $maxStart) && ($curr_ref->{"repStart"} > $maxStart)) || ((defined $maxLeft) && ($curr_ref->{"repLeft"} > $maxLeft)) || ((defined $minlength) && ($curr_ref->{"chrEnd"} - $curr_ref->{"chrStart"} + 1 < $minlength))){ #filter out maxstart, maxend and minlength parameter
		#print "the last record!!!\n";
		#print "$curr_ref->{'chr'}\t$curr_ref->{'chrStart'}\t$curr_ref->{'chrEnd'}\t$curr_ref->{'strand'}\nstart is $curr_ref->{'repStart'}\nend is $curr_ref->{'repLeft'}\n\n";
		
		#assign last bio::seq to alignment
		my $chrSeq = "-"x($curr_ref->{"repStart"}-1).$curr_ref->{"chrSeq"}."-"x($curr_ref->{"repLeft"});
#		print "$headline\n$chrSeq\n";
		my $tempId = $curr_ref->{"chr"}."_".$curr_ref->{"chrStart"}."-".$curr_ref->{"chrEnd"}.$curr_ref->{"strand"};
		my $tempSeq = Bio::LocatableSeq->new(-seq => $chrSeq,
							-id => $tempId,
							-alphabet => "dna",
						);
		#add the chrSeq with insertions deleted to the MSA:
		if (defined $msa){
			$msa->add_seq($tempSeq);
			for my $indels ($curr_ref->{"indel_matrix"}){
				for my $single_insertion (@$indels){
					$msa_indels{$single_insertion->{'pos'} + $curr_ref->{"repStart"}-1}{$tempId} = $single_insertion->{'seq'};
				}
			}
		}
		else{
			$msa = Bio::SimpleAlign->new(-seqs => [$tempSeq]);
			for my $indels ($curr_ref->{"indel_matrix"}){
				for my $single_insertion (@$indels){
					$msa_indels{$single_insertion->{'pos'} + $curr_ref->{"repStart"}-1}{$tempId} = $single_insertion->{'seq'};
				}
			}
		}
	}
	else {
		#print "the last record not included!!!\n";
		#print "not included:\n$curr_ref->{'chr'}\t$curr_ref->{'chrStart'}\t$curr_ref->{'chrEnd'}\t$curr_ref->{'strand'}\nstart is $curr_ref->{'repStart'}\nend is $curr_ref->{'repLeft'}\n\n";
	}
	# sort $msa_indels key descending
	# foreach key, insert $msa_indels->{"seq"} to associated $seq_obj->id. if use the longest length, and shorter seq supplimented with gap.
	if ($opts{i}){
		for my $pos (sort {$b <=> $a} keys %msa_indels){
			# maxium length
			my @all_indel_lengths = map {length($msa_indels{$pos}{$_})} keys %{$msa_indels{$pos}};
			my $maxium_length = max(@all_indel_lengths);
			if ($maxium_length > 0){
				my $temp_insertionMSA = Bio::SimpleAlign->new();
				foreach my $seq_obj($msa->each_seq){
					my $temp_displayId = $seq_obj->display_id();
					my $indel_seq;
					if (exists $msa_indels{$pos}{$temp_displayId}) {
						$indel_seq = $msa_indels{$pos}{$temp_displayId};
					}
					else{
						$indel_seq = "";
					}
					my $insertion_seq = '-'x($maxium_length-length($indel_seq)).$indel_seq;
					my $tempSeq = Bio::LocatableSeq->new(-seq => $insertion_seq,
										-id => $temp_displayId,
										-alphabet => "dna",
									);
					$temp_insertionMSA->add_seq($tempSeq);
				}
				$msa = cat($msa->slice(1,$pos-1,1), $temp_insertionMSA, $msa->slice($pos,$msa->length(),1))
			}
		}
	}
	print "$te alignment done!\n";
	$pm->finish unless defined $msa;
	unless ($msa->is_flush){ # quit and print error message if not all seqs in alignment have same length.
		print "$te alignment is not flush!!\n";
		my @all_length = ();
		foreach my $seq_obj($msa->each_seq){
			push @all_length, $seq_obj->length();
		}
		my $median = median(@all_length);
		foreach my $seqobj($msa->each_seq){
			unless ($seqobj->length() == $median){
				print STDERR "removed seq: $te\t".$seqobj->display_id()."\t";
				print STDERR $seqobj->length()."\n";
				print STDERR $seqobj->seq()."\n";
				$msa->remove_seq($seqobj);
			}
		}
		#die "$te alignment is not flush!!\n";
	}
	my $consensus = $msa->consensus_string();
	my $consensus_obj = Bio::LocatableSeq->new(-seq => $consensus,
						-id => "consensus",
						-alphabet => "dna",
					);
	my $copies = $msa->num_sequences();
	my $length = $msa->length();
	flock $newFasta, LOCK_EX or die "can't lock!!" if $purpose eq "con" || $purpose eq "all";
	print $newFasta ">$te\n$consensus\n" if $purpose eq "con" || $purpose eq "all";
	flock $newFasta, LOCK_UN or die "can't unlock!!" if $purpose eq "con" || $purpose eq "all";

#	push (@newFasta, [$te,$consensus]);
	if ($purpose eq "aln"){
		foreach my $seq_obj($msa->each_seq){
			my $temp_id = $seq_obj->display_id();
			my $temp_seq = $seq_obj->seq();
			flock $newFasta, LOCK_EX or die "can't lock!!";
			print $newFasta ">$temp_id\n$temp_seq\n";
			flock $newFasta, LOCK_UN or die "can't unlock!!";
#			push (@newFasta, [$temp_id,$temp_seq]);
		}
	}
	$pm->finish if $copies == 1;
	if ($purpose eq "CG" || $purpose eq "all"){
		my $numCpG = 0; #number of all CpG sites in repSeq
		my $mutTpG = 0; #number of mutated CpG to TpG
		my $mutCpA = 0; #number of mutated CpG to CpA
		my $numC = 0; #number of all C sites and G sites in repSeq (excluding CpG sites)
		my $mutC = 0; #number of all C to T transitions in repSeq (excluding CpG sites)
		my $numG = 0; #number of all G sites in repseq
		my $mutG = 0; #number of all G to A transitions in repseq
		if ($model eq "c"){
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
		}
		if ($model eq "a"){
			##write fasta to $fasta
			my ($filehandle, $fasta) = tempfile(DIR => '/home/xiaoyu/temp');
			my $out = Bio::AlignIO->new(-file => ">$fasta",
							-format =>'fasta',
							-displayname_flat => 1, #to avoid an AlignIO bug.
						);
			$out->write_aln($msa);
			my $prank_hash_ref = {"d" => $fasta,
						"showanc" => 1,
						"showevents" => 1,
						"keep" => 1,
						"uselogs" => 1,
						"quiet" => 1,
						"o" => "$fasta.CpG.out",
			};
			($numC,$mutC,$numG,$mutG,$numCpG,$mutTpG,$mutCpA) = Func::CpG_with_prank($prank_hash_ref,0);
			unlink($fasta);
		}
		flock $outFile, LOCK_EX or die "can't lock!!";
		print $outFile "$te\t$copies\t$numC\t$mutC\t$numG\t$mutG\t$numCpG\t$mutTpG\t$mutCpA\t$length\n";
		flock $outFile, LOCK_UN or die "can't unlock!!";
#		push (@outFile, "$te\t$copies\t$numC\t$mutC\t$numG\t$mutG\t$numCpG\t$mutTpG\t$mutCpA"); 
	}
	print "$te analysis done!!!\n";
	$pm->finish;
}
for my $dbfile (qw(mldbm)){
	unlink $dbfile;
	print "hmm, unable to delete dbfile!\n\n" if -e $dbfile;
}
	
$pm->wait_all_children;
print "all threads done!!!\n\n";

#if (defined $opts{d}){
#	unlink($tempfasta,$tempfasta.".index");
#}

if ($purpose eq "CG" || $purpose eq "all"){
#	open $outFile, ">$outName";
#	print $outFile "te\tcopies\tnumC\tmutC\tnumG\tmutG\tnumCpG\tmutTpG\tmutCpA\n";
#	foreach (@outFile){
#		print $outFile "$_\n";
#	}
	close ($outFile);
}
if ($purpose eq "con" || $purpose eq "aln" || $purpose eq "all"){
#	open $newFasta, ">$outFasta";
#	foreach (@newFasta){
#		print $newFasta ">$_->[0]\n$_->[1]\n";
#	}
	close ($newFasta);
}
#close($outFile) if ($purpose eq "CG" || $purpose eq "all");
#close($newFasta) if ($purpose eq "con" || $purpose eq "aln" || $purpose eq "all");

exit;

