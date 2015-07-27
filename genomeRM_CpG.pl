##################################################

=head1 description

Author  : Xiaoyu Zhuo
Purpose : Read a repeatmasker alignment output and regenerate the consensus sequence and then calculate CpG/nonCpG substitution ratio. Note: insertions are ignoared during regeneration. Directly modified from trueCONSENSUS
modified from John Pace's script.

=head2 USAGE

perl genomeRM_CpG.pl -a <RM align file>  [-p <con|CG|aln|all>] [-c <class name>] [-r <repeat name>] [-s <INT> -l <INT>] [-h help]

options:
-a the repeatmasker alignment file.

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

ps: if -s and -l are not defined, any fragment matches particular TE would be included in the alignment for calculation.

-h
help

=cut

use FindBin qw($RealBin);
use lib $RealBin;
use warnings;
use strict;
use Func;
use Getopt::Std;
use Bio::SimpleAlign;
use Bio::LocatableSeq;

my %opts=();
getopts("ha:p:c:r:s:l:", \%opts);
my $usage = "perl ParseRMalign_regenerate_consensus.pl -a <RM align file>  [-p <con|CG|aln|all>] [-c <class name>] [-r <repeat name>] [-s <INT> -l <INT>]";
die "$usage" if $opts{h};
my $fileName = $opts{a} or die "$usage";
$opts{p} ||= "con"; #set default parameter for p: only output consensus.
die "$usage" unless $opts{p} eq "con" || $opts{p} eq "CG" || $opts{p} eq "all" || $opts{p} eq "aln";
my $purpose = $opts{p};
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

my @alignArray = ();
print "Reading and processing ".$fileName."\n";

#Open the align file
open my $alignfile, "<$fileName";
my $outName;
my $outFasta;
if (defined $targetClass || defined $targetRepeat || defined $maxStart || defined $maxLeft){
	my $target = "";
	$target = $target.$targetClass if defined $targetClass;
	$target = $target.$targetRepeat if defined $targetRepeat;
	$target = $target.".$maxStart-$maxLeft" if defined $maxStart && defined $maxLeft;
	$outName = "$fileName.$target.CpG.out";
	$outFasta = "$fileName.$target.new_consensus.fa";
}
else{
	$outName = "$fileName.CpG.out";
	$outFasta = "$fileName.new_consensus.fa";
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

#a HASH holding each repeat with its MSA
my %repHash = ();

#Loop through the array, parsing the lines and add to an output file
for(my $i = 0; $i<= $#alignArray;$i++){
	print "$i\n" if (($i % 10000) == 0); #print line number of every 10000 lines

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
			#repClass has to be split further (modified from John Pace's)
			($repClass, $repFamily) = split(/\//, $repClass);
		}

		else { #The repeat is in + orientation
			#Set the variables
			$genoName = $arraySplit[4];
			$genoStart = $arraySplit[5];
			$genoEnd = $arraySplit[6];
			$strand = '+';
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

		next if ((defined $targetClass) && ($repClass ne $targetClass)) || ((defined $targetRepeat) && ($repName ne $targetRepeat)) || ((defined $maxStart) && ($repStart > $maxStart)) || ((defined $maxLeft) && ($repLeft > $maxLeft)); #only do this with the $targetClass or $targetRepeat and repeat satisfy $maxStart and $maxLeft.

		#All the variables are set now, except chrSeq and repSeq.
		for (my $k = ($beginningRow +1);$k<=($endingRow-1);$k++){
			#Split alignArray[$k] to see if it is a sequence line.
			#First check for + strand
			@arraySplit = ();  #Clear arraySplit
			@arraySplit = split(/\t/, $alignArray[$k]);
			if (($arraySplit[0] ne 'Matrix') &&($arraySplit[0] ne 'Transitions') && ($arraySplit[0] ne 'Gap_init') && ($arraySplit[0] ne 'Kimura') && ($strand eq '+')){
				#Check to see if this is a chr sequence or a rep sequence
				if ( $arraySplit[0] eq substr($genoName, 0,13) || $arraySplit[0] =~ /^\Q$genoName\E/){	#This is a chr seq
					$chrSeq = $chrSeq.$arraySplit[2];
				}
				else{  #This is the rep seq
					$repSeq = $repSeq.$arraySplit[2];
				}
			}

			#Next check for the - (C) strand
			if (($arraySplit[0] ne 'Matrix') && ($arraySplit[0] ne 'Transitions') && ($arraySplit[0] ne 'Gap_init') && ($arraySplit[0] ne 'Kimura') && ($strand eq 'C')){
				#Check to see if this is a chr sequence or a rep sequence
				if ($arraySplit[0] eq substr($genoName, 0,13) || $arraySplit[0] =~ /^\Q$genoName\E/){	#This is a chr seq
					$chrSeq = $chrSeq.$arraySplit[2];
				}
				else{  #This is a rep seq
					$repSeq = $repSeq.$arraySplit[3];
				}
			}
		}
		$chrSeq = Func::revcom($chrSeq) if $strand eq 'C';
		$repSeq = Func::revcom($repSeq) if $strand eq 'C';

		#delete insertions and create new mutiple alignment
		$chrSeq = uc($chrSeq);
		$repSeq = uc($repSeq); 
		my @gaps = Func::find($repSeq, "-"); #gap in repSeq means insertions in chrSeq, we have to delete them!
		my @descending_gaps = sort { $b<=>$a } @gaps;
		foreach my $gap(@descending_gaps){
			substr($chrSeq, $gap-1, 1) = "";
		}
		next if $chrSeq eq ""; #remove suspicious sequences (I don't know why they are in the alignment in the first place. but the fact that they have different seq length causes a problem for MSA building).
		$chrSeq = "-"x($repStart-1).$chrSeq."-"x($repLeft);
#		print "$headline\n$chrSeq\n";
		my $tempSeq = Bio::LocatableSeq->new(-seq => $chrSeq,
							-id => "$genoName".":$genoStart"."-$genoEnd$strand",
							#	-start => $repStart,
							#	-end => $repEnd,
							-alphabet => "dna",
						);
		#add the chrSeq with insertions deleted to the MSA:
		if (exists $repHash{$repName}){
			$repHash{$repName}->add_seq($tempSeq);
		}
		else{
			my $tempAln = Bio::SimpleAlign->new(-seqs => [$tempSeq]);
			$repHash{$repName} = $tempAln;
		}
	}
}
undef @alignArray;

#open the outfile
my $outFile;
my $newFasta;
if ($purpose eq "CG" || $purpose eq "all"){
	open $outFile, ">$outName";
	print $outFile "te\tcopies\tnumC\tmutC\tnumG\tmutG\tnumCpG\tmutTpG\tmutCpA\n";
}
if ($purpose eq "con" || $purpose eq "aln" || $purpose eq "all"){
	open $newFasta, ">$outFasta";
}
#parse repHash to regenerate consensus for all repeats:
while (my ($te, $msa) = each(%repHash)){
	print "working on $te:\n";
	unless ($msa->is_flush){ # quit and print error message if not all seqs in alignment have same length.
		foreach my $seq_obj($msa->each_seq){
			print $seq_obj->display_id()."\t";
			print $seq_obj->length()."\n";
			print $seq_obj->seq()."\n";
		}
		die "$te alignment is not flush!!\n";
	}
	my $consensus_obj = $msa->consensus_meta();
	my $consensus = $consensus_obj->seq();
	my $copies = $msa->num_sequences();
	print $newFasta ">$te\n$consensus\n" if $purpose eq "con" || $purpose eq "all";
	if ($purpose eq "aln"){
		foreach my $seq_obj($msa->each_seq){
			my $temp_id = $seq_obj->display_id();
			my $temp_seq = $seq_obj->seq();
			print $newFasta ">$temp_id\n$temp_seq\n";
		}
	}
	if ($purpose eq "CG" || $purpose eq "all"){
		my $numCpG = 0; #number of all CpG sites in repSeq
		my $mutTpG = 0; #number of mutated CpG to TpG
		my $mutCpA = 0; #number of mutated CpG to CpA
		my $numC = 0; #number of all C sites and G sites in repSeq (excluding CpG sites)
		my $mutC = 0; #number of all C to T transitions in repSeq (excluding CpG sites)
		my $numG = 0; #number of all G sites in repseq
		my $mutG = 0; #number of all G to A transitions in repseq
		foreach my $genoSeq_obj ($msa->each_seq){
			my ($C, $T, $G, $A, $CpG, $TpG, $CpA) = Func::CpG_rate($consensus_obj, $genoSeq_obj);
			$numC = $numC + $C;
			$mutC = $mutC + $T;
			$numG = $numG + $G;
			$mutG = $mutG + $A;
			$numCpG = $numCpG + $CpG;
			$mutTpG = $mutTpG + $TpG;
			$mutCpA = $mutCpA + $CpA;
		}
		print $outFile "$te\t$copies\t$numC\t$mutC\t$numG\t$mutG\t$numCpG\t$mutTpG\t$mutCpA\n";
	}
}

close($outFile) if ($purpose eq "CG" || $purpose eq "all");
close($newFasta) if ($purpose eq "con" || $purpose eq "all");

exit;

