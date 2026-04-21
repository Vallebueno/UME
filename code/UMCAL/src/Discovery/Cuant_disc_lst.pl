$lst="/scratch-cbe/users/miguel.vallebueno/MERGE/SPLIT/disc_ll.lst";


open(LST,"$lst") or die;
$out="$lst.cuant2";
open(OUT,">$out") or die;

print OUT "chr&cof\ttotal\tbialelic\ttrialelic\ttetralelic\n";

while($file=<LST>) {

    @xcol = split(/\./, $file);



    open(FLE,"$file") or die;

    $c=0;
    undef(%hc);

    #print "START <$xcol[4]><$xcol[12]> <$hc{2}>\n";

    while($line=<FLE>) {

        chomp($line);
        $c++;
        @lcol = split(/\t/, $line);
        @col = split(/\s/, $lcol[2]);
        $len=scalar(@col);
        #print "<$len> <$lcol[2]> @lcol"; <STDIN>;
        $hc{$len}=$hc{$len}+1;

    }

    #print "END print t<$c> chr<$xcol[4]> cof<$xcol[12]>  W0<$hc{0}> W1<$hc{1}> W2<$hc{2}> W3<$hc{3}> W4<$hc{4}> W5<$hc{5}> W6<$hc{6}>\n";

    $hfc{$xcol[12]}=$hfc{$xcol[12]}+$c;
    $hfb{$xcol[12]}=$hfb{$xcol[12]}+$hc{2};
    $hft{$xcol[12]}=$hft{$xcol[12]}+$hc{3};
    $hftt{$xcol[12]}=$hftt{$xcol[12]}+$hc{4};

    $hh{$xcol[12]}=1;
    close(FLE);
}


foreach my $aa (sort keys %hh){

    print OUT "$aa\t$hfc{$aa}\t$hfb{$aa}\t$hft{$aa}\t$hftt{$aa}\n";

}




close(OUT);