package Func;
use strict;
use warnings;

sub which_mutation{ #return 0 for same residue, 1 for transition, 2 for transversion.
	;
	my $res1 = shift;
	my $res2 = shift;
	my %res = ( #standard IUPAC single letter symbol
		"G" => "purine",
		"A" => "purine",
		"T" => "pyrimidine",
		"C" => "pyrimidine",
		"U" => "pyrimidine",
		"R" => "purine",
		"Y" => "pyrimidine",
	);
	my $result;
	if($res1 eq $res2){
		$result = 0;
	}
	else{
		$result = ($res{$res1} eq $res{$res2})?1:2;
	}
	return $result;
}

sub revcom{
	my $seq = $_[0];
	my $revcom = reverse($seq);
	$revcom =~ tr/ACGTRYMKWSBDHVacgtrymkwsbdhv/TGCAYRKMSWVHDBtgcayrkmswvhdb/;
	return $revcom;
}

sub find{
	my @pos = ();
	for(my $i=0;$i >= 0;){
		my $pos = index($_[0], $_[1], $i); 
		if($pos >= 0){
			$i = $pos+1;
			push(@pos, $i);
		}
		else{
			$i=$pos;
		}
	}
	return \@pos;
}

sub CpG_rate{
	my $anc_obj = shift;
	my $cur_obj = shift;
	my $pairwise = Bio::SimpleAlign->new(-seqs => [$anc_obj, $cur_obj],
						);
	my $pairwise_aln = $pairwise->remove_gaps("-",1); #flag 1 means only remove all gap columns.
	$anc_obj = $pairwise_aln->get_seq_by_pos(1);
	$cur_obj = $pairwise_aln->get_seq_by_pos(2);
	my $anc_seq = $anc_obj->seq();
	my $cur_seq = $cur_obj->seq();
	my $numCpG = 0; #number of all CpG sites in repSeq
	my $mutTpG = 0; #number of mutated CpG to TpG
	my $mutCpA = 0; #number of mutated CpG to CpA
	my $numC = 0; #number of all C sites and G sites in repSeq (excluding CpG sites)
	my $mutC = 0; #number of all C to T transitions in repSeq (excluding CpG sites)
	my $numG = 0; #number of all G sites in repseq
	my $mutG = 0; #number of all G to A transitions in repseq
	my $CpG_ref = &find($anc_seq, "CG");
	my $CpA_ref = &find($anc_seq, "CA");
	my $TpG_ref = &find($anc_seq, "TG");
	my $C_ref = &find($anc_seq, "C");
	my $G_ref = &find($anc_seq, "G");
	foreach my $CpG(@$CpG_ref){
		my $res = substr($cur_seq, $CpG-1, 2);
		if ($res =~ m/^[^-][^-]$/){
			$numCpG++;
			$mutTpG++ if ($res eq "TG");
			$mutCpA++ if ($res eq "CA");
		}
	}
	foreach my $G(@$G_ref){
		unless($G-1 ~~ @$CpG_ref or $G-1 ~~ @$TpG_ref){
			my $res = substr($cur_seq, $G-1, 1);
			$numG++ if $res =~ m/^[^-]$/;
			$mutG++ if $res eq "A"; 
		}
	}
	foreach my $C(@$C_ref){
		unless($C ~~ @$CpG_ref or $C ~~ @$CpA_ref){
			my $res = substr($cur_seq, $C-1, 1);
			$numC++ if $res =~ m/^[^-]$/;
			$mutC++ if $res eq "T";
		}
	}
	return($numC, $mutC, $numG, $mutG, $numCpG, $mutTpG, $mutCpA);
}

#the backup is not used for calculation now...
sub CpG_rate_backup{ #CpG mutation/all CpG site, CpG to TpA count as 2 mutation
	my $anc_obj = shift;
	my $cur_obj = shift; #both $anc_obj and $cur_obj are aligned Bio::LocatableSeq object
	my ($site, $CpGcriteria, $nonCpGcriteria) = @_;
	$site ||= "all";
	$CpGcriteria ||= "stringent";
	$nonCpGcriteria ||= "all";
	my $pairwise = Bio::SimpleAlign->new(-seqs => [$anc_obj, $cur_obj],
						);
	my $pairwise_aln = $pairwise->remove_gaps("-",1); #flag 1 means only remove all gap columns.
	$anc_obj = $pairwise_aln->get_seq_by_pos(1);
	$cur_obj = $pairwise_aln->get_seq_by_pos(2);
	my $anc_seq = $anc_obj->seq();
	my $cur_seq = $cur_obj->seq();
	my $CpG_ref = find($anc_seq, "CG"); #an array ref of collumn position where "CG" is found.
	my $numCpG = 0;
	my $numTpG = 0;
	my $numCpA = 0;
	my $numTpA = 0;
	my $numTpY = 0;
	my $numRpA = 0;
	foreach my $CpG (@$CpG_ref){
		my $res = $cur_obj->subseq($CpG, $CpG+1);
		$numCpG++ if $res =~ m/[^-][^-]/; #if there is no gap in current seq positions.
		$numTpG++ if $res eq "TG";
		$numCpA++ if $res eq "CA";
		$numTpA++ if $res eq "TA";
		$numTpY++ if $res =~ "T[TC]";
		$numRpA++ if $res =~ "[AG]A";
	}
	my $CGmut; #mutated CpG site.
	if ($site eq "C"){
		$CGmut = ($CpGcriteria eq "stringent")?$numTpG:($numTpG+$numTpA+$numTpY);
	}
	elsif ($site eq "G"){
		$CGmut = ($CpGcriteria eq "stringent")?$numCpA:($numCpA+$numTpA+$numRpA);
	}
	else{
		$CGmut = $numTpG + $numCpA if $CpGcriteria eq "stringent";
		$CGmut = $numTpG + $numCpA + $numTpA + $numTpY + $numRpA if $CpGcriteria eq "one";
		$CGmut = $numTpG + $numCpA + $numTpA*2 + $numTpY + $numRpA if $CpGcriteria eq "two";
	}
	my $lenCpG = ($site eq "C" or $site eq "G")?$numCpG:2*$numCpG; #total length of all CpG sites in ancestor seq.
	my $total_len = $pairwise_aln->length();
	my @anc_char = split //, $anc_seq;
	my @cur_char = split //, $cur_seq;
	my $numC = 0; #number of C in ancestor
	my $mutationC = 0; #number of C mutated in current node
	my $transitionC = 0; #number of C transited to T in current node
	my $numG = 0; #number of G in ancestor
	my $mutationG = 0; #number of G mutated in current node
	my $transitionG = 0; #number of G transited to A in current node
	my $numW = 0; #number of A/T in ancestor
	my $mutationW = 0; #number of A/T mutated in current node
	my $transitionW = 0; #number of A/T transited to G/C in current node
	for (my $i = 0;$i<$total_len;$i++){
		unless($anc_char[$i] eq "-" or $cur_char[$i] eq "-" or $i ~~ @$CpG_ref or $i+1 ~~ @$CpG_ref){
			if ($anc_char[$i] eq "C"){
				$numC++;
				my $mutation_type = which_mutation($anc_char[$i],$cur_char[$i]);
				if($mutation_type){
					$mutationC++;
					$transitionC++ if $mutation_type == 1;
				}
			}
			elsif($anc_char[$i] eq "G"){
				$numG++;
				my $mutation_type = which_mutation($anc_char[$i],$cur_char[$i]);
				if($mutation_type){
					$mutationG++;
					$transitionG++ if $mutation_type == 1;
				}
			}
			else{
				$numW++;
				my $mutation_type = which_mutation($anc_char[$i],$cur_char[$i]);
				if($mutation_type){
					$mutationW++;
					$transitionW++ if $mutation_type == 1;
				}
			}
		}
	}
	my $length; #total length of nonCpG site under current setting;
	my $mut; #total mutation in nonCpG site under current setting;
	if ($site eq "C"){
		$length = $numC;
		$mut = $mutationC if $nonCpGcriteria eq "all";
		$mut = $transitionC if $nonCpGcriteria eq "transition";
	}
	elsif ($site eq "G"){
		$length = $numG;
		$mut = $mutationG if $nonCpGcriteria eq "all";
		$mut = $transitionG if $nonCpGcriteria eq "transition";
	}
	elsif ($site eq "CG"){
		$length = $numC + $numG;
		$mut = $mutationC + $mutationG if $nonCpGcriteria eq "all";
		$mut = $transitionC + $transitionG if $nonCpGcriteria eq "transition";
	}
	else{
		$length = $numC + $numG + $numW;
		$mut = $mutationC + $mutationG + $mutationW if $nonCpGcriteria eq "all";
		$mut = $transitionC + $transitionG + $transitionW if $nonCpGcriteria eq "transition";
	}
	return ($CGmut, $lenCpG, $mut, $length);
}
1;
