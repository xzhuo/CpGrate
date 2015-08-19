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

sub CpG_with_prank{
	use IO::String;
	use Bio::TreeIO;
	#my $prankpath =  "/Users/xiaoyu/bioinfo/prank/bin/";
	#my $prankpath = "/Users/xiaoyu/bioinfo/prank-msa/src/";
	#my $prankpath = "/data2/xiaoyu/test_perlCpG/prank/bin";
	my $prankpath = "/home/xiaoyu/prank/bin/";
	local $ENV{PATH} = "$ENV{PATH}:$prankpath";
	
	my $prank_hash_ref = shift or die $!; #prank parameters
	my $keepfile = shift;
#	my $sub_root_id = shift;
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
#		open $out, ">$domain".".CpG.out" or die "$!";
	}
	else{
		$fasta = $domain.".best.anc.fas";
		$events = $domain.".best.events";
		$tree_file = $domain.".best.anc.dnd";
#		open $out, ">$domain".".best.CpG.out" or die "$!";
	}
	my $str = Bio::AlignIO->new(-file => $fasta,
					-format => 'fasta',
				) or die "\nproblem running prank with $fasta\n";;
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
#	my $sub_root;
#	my @sub_tree_nodes;
#	my @sub_ids;
#	if ($sub_root_id){    #if sub_root_id is provided, do the following things:
#		$sub_root = $tree->find_node($sub_root_id);
#		@sub_tree_nodes = $sub_root->get_all_Descendents;
#		foreach (@sub_tree_nodes) {
#			push @sub_ids, $_->id;
#		}
#	}

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
#			print "$ancestor_id\t$current_node_id\n";
			my $anc_obj = $aln->get_seq_by_id($ancestor_id);
			my $cur_obj = $aln->get_seq_by_id($current_node_id);
#			print "$anc_obj\t$cur_obj\n";
			my ($numC, $mutC, $numG, $mutG, $numCpG, $mutTpG, $mutCpA) = Func::CpG_rate($anc_obj,$cur_obj);
			$all_numC += $numC;
			$all_mutC += $mutC;
			$all_numG += $numG;
			$all_mutG += $mutG;
			$all_numCpG += $numCpG;
			$all_mutTpG += $mutTpG;
			$all_mutCpA += $mutCpA;
#			if($sub_root_id && $current_node_id ~~ @sub_ids){ #if sub_root_id is provided, do the following things:
#				$sub_numC += $numC;
#				$sub_mutC += $mutC;
#				$sub_numG += $numG;
#				$sub_mutG += $mutG;
#				$sub_numCpG += $numCpG;
#				$sub_mutTpG += $mutTpG;
#				$sub_mutCpA += $mutCpA;
#			}

			my $branchratio;
			eval { $branchratio = ($mutTpG/$numCpG)/($mutC/$numC);};
			$branchratio = $@?"NA":sprintf("%.2f",$branchratio); #leave only 2 digits.
			$tree_string =~ s/$current_node_id/$current_node_id $branchratio/;
			#$node->description($branchratio);
#			print $out "$current_node_id\t$numC\t$mutC\t$numG\t$mutG\t$numCpG\t$mutTpG\t$mutCpA\n";
		}
	}
	#my $treeout = Bio::TreeIO->new(-file => ">$domain.CpG.new",
	#				-format => 'newick',
	#			);
	#$treeout->write_tree($tree);
	#print Dumper($tree);
#	print $out "$tree_string\n";
#	print $out "$all_numC\t$all_mutC\t$all_numG\t$all_mutG\t$all_numCpG\t$all_mutTpG\t$all_mutCpA\n";
#	my $ratio;
#	eval { $ratio = ($all_mutTpG/$all_numCpG)/($all_mutC/$all_numC);};
#	$ratio = "NA" if $@;
#	printf $out "%.2f", $ratio;
#	print $out "\n";

#	if($sub_root_id){  #if sub_root_id is provided, do the following things:
#		my $subratio;
#		eval {$subratio = ($sub_mutTpG/$sub_numCpG)/($sub_mutC/$sub_numC);};
#		$subratio = $@?"NA":sprintf("%.2f",$subratio);
#		print $out "$sub_root_id\n$sub_numC\t$sub_mutC\t$sub_numG\t$sub_mutG\t$sub_numCpG\t$sub_mutTpG\t$sub_mutCpA\n";
#		print $out "$sub_root_id subratio is $subratio\n";
#	}
	close ($Fh);
#	close ($out);
	unlink($fasta, $events, $tree_file, "$domain.fas") if $keepfile == 0;
	return ($all_numC,$all_mutC,$all_numG,$all_mutG,$all_numCpG,$all_mutTpG,$all_mutCpA);
}


#how to build the prank_hash: "d" => "some.fas" for -d=some.fas, "F" =>0 for -F in command line.
sub run_prank{ # argument is the prank argument hash ref.
	my $prank = "prank";
	my %inhash = %{shift()};
	my @params = ($prank);
	my $fasta = $inhash{d};
	#rename sequence nome in fasta file to bypass a prank bug.
#	renamefasta($fasta);
#	$inhash{d} = "$fasta.new.fa";
	
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
#	unlink "$fasta.new.fa";
}

#rename fasta seq name to work with phylip alignment file restriction.
sub renamefasta{
	my $infile = shift;
	my %namehash;
	open my $in, "<$infile" or die "$!";
	open my $out, ">$infile.new.fa" or die "$!";
	my $newname = "aaaaaaaaaa";
	while (<$in>) {
		if (/^>/){
			my ($name) = ($_ =~ /^>(.+)\n$/);
			die "too many sequences for the renamefasta sub!!!" if $newname eq "aaaaaaaaaaa";
			$namehash{$newname} = $name;
			print $out ">$newname\n";
			$newname++;
		}
		else{
			print $out "$_";
		}
	}
	return(\%namehash);
}
1;
