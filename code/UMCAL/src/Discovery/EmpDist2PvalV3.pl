#$list="/groups/swarts/lab/Resources/Maize/WGS/gVCFS/ALL/files.mergein.lst.Disc.in";
#$file="/scratch-cbe/users/miguel.vallebueno/MERGE/All.tlone.Merge.db.10.poli.db.gz";
$list=$ARGV[0];
$file=$ARGV[1];

open(LST,"$list") or die;

$c=3;

$cfl=10000000;
print "Reading list....\n";
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $name=$lcol[0];
    $caller=$lcol[1];
    $key="$name&$caller";
    $head[$c]=$key;
    $c++;
}
$tot=$c;

close($LST);

open(FIL,"gunzip -c $file |") or die;
$out="$file.pval2f";
open(OUT,">$out");
print OUT "#Caller\tPhred\tpval\n";
print "Reading distributions....\n";

$cfc=0;
while($line=<FIL>) {
    chomp($line);
    @lcol = split(/\t/, $line);

    $c=0;
    foreach my $aa (@lcol){
       if($c<3){$c++;next;}
        if($aa eq ":::"){$c++;next;}
        if($aa =~ ":.:"){$c++;next;}
        if($aa eq ""){$c++;next;}

        @xx= split(/:/, $aa);
        $vv=sprintf("%.1f", $xx[1]);

        $hr=$head[$c];
        @hhrr= split(/&/, $hr);
        $name=$hhrr[0];
        $caller=$hhrr[1];
   #     print "<$c> <$name> <$caller> <$aa> <$vv>\n";<STDIN>;

        $hd{$caller}{$vv}=$hd{$caller}{$vv}+1;
        $ht{$caller}=$ht{$caller}+1;
        $cfc++;
     $c++;
    }
   # print "<$cfc>   <$cfl>\n";

    if($cfc>$cfl){last;}
}
close(FIL);
print "Done!\n";


print "Creating empirical pval table....\n";
foreach my $caller (keys %hd){
    print "$caller\n";
   # @sort= sort { $a <=> $b } @{$stat{$caller}};

    $cont=0;
    foreach my $val (sort {$b <=> $a} keys %{$hd{$caller}}){
        $cont=$cont+$hd{$caller}{$val};
        $nval=$cont/$ht{$caller};

     #   if($val < 20){print  "$caller\t$val\t<$nval>\t<$cont> <<$hd{$caller}{$val}>>\n"; <STDIN>;}

        print OUT "$caller\t$val\t$nval\n";
    }
}
print "Done!\n";
close(OUT);

