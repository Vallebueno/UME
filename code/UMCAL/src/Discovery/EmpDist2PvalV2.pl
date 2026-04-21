$file="/groups/swarts/lab/MAVE/phred/FIN_SUBSET_1.ALL.Merge.db.phredquant";
#$file="/groups/swarts/lab/MAVE/phred/Emp_Dist_sum.txt";

open(LST,"$file") or die;
$out="$file.pval2f";
open(OUT,">$out");
#open(OUTA,">$out.test");
print OUT "#Caller\tPhred\tpval\n";
print "Reading distributions....\n";
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    if($lcol[0] eq ""){next;}

$vv=sprintf("%.2f", $lcol[1]);
    $hd{$lcol[0]}{$vv}=$hd{$lcol[0]}{$vv}+1;
    $ht{$lcol[0]}=$ht{$lcol[0]}+1;
    #push(@{$stat{$lcol[0]}}, $lcol[1]);
}
close(LST);
print "Done!\n";
print "Creating empirical pval table....\n";
foreach my $caller (keys %hd){
    print "$caller\n";
   # @sort= sort { $a <=> $b } @{$stat{$caller}};

    $cont=0;
    foreach my $val (sort {$b <=> $a} keys %{$hd{$caller}}){
        $cont=$cont+$hd{$caller}{$val};
        $nval=$cont/$ht{$caller};
        #print OUTA "$caller\t$val\t<$nval>\t<$cont> <<$hd{$caller}{$val}>>\n";
        print OUT "$caller\t$val\t$nval\n";
    }
}
print "Done!\n";
close(OUT);

