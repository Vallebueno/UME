$list=$ARGV[0];  #k
$ref=$ARGV[1];   #r

$out=$list;
if($out =~ /\.lst2$/){$out =~ s/\.lst2$/.ll/;}
elsif($out =~ /\.lst$/){$out =~ s/\.lst$/.ll/;}
else{$out =~ s/lst5/ll/g;}

open(OUT,">$out") or die;

print "FIle:<$list>\n";
print "Reference:<$ref>\n";
open(REF,"$ref") or die;

print "loading ref...\n";
chomp(my @re = <REF>); #magic
$f=0;
foreach my $aa(@re){
    if($f==1){$href{$name}=$aa;$f=0;next;}
    if($f==0){$name=$aa; $name=~ s/>//; print "$name\n";$f=1;$ctrl{$name}=1;next;}
}
close(REF);
undef(@re);
print "Done!...\n";


open(LST,"$list") or die;
print "Sorting alts to VCF genotype...\n";
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $chr=$lcol[0];
    #$chr=10;
    $pos=$lcol[1];
    $chr=~ s/ //g;
    $pos=~ s/ //g;
    $lcol[2]=~ s/^ //;
    if($lcol[2] =~ /N/){next}
    if($lcol[2] =~ /\./){next}
       $lcol[2] =~ s/,/ /;
    @vars=split(/ /,$lcol[2]);
    $gref=substr $href{$chr},$pos-1, 1;

    undef(@alts);
    undef(@altx);
    push(@alts,$gref);
    foreach my $ee (@vars){
        if($ee eq $gref){next;}
        push(@alts,$ee);
        push(@altx,$ee);
    }
    $asys=join(",", @altx);
 #    print "<$chr> <$pos> vars <@vars> ref<$gref>  alts<@alts>  altx<$asys>"; <STDIN>;
    print OUT "$chr\t$pos\t.\t$gref\t$asys\t.\tPASS\t.\tGT\n";

}
close(LST);
close(OUT);
