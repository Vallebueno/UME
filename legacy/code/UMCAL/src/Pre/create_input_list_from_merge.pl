#$alg=$ARGV[0];
#$headf="$alg.header";
$headf=$ARGV[0];
$alg=$headf;
open(OUT,">$alg.keyx1") or die;
open(OUTA,">$alg.keyx2") or die;
print "Extracting headers from input list...\n";
open(HED,"$headf") or die;
chomp(my @LLL = <HED>);
@lcol = split(/\|/, $LLL[0]);
$head=$lcol[0];
$head =~ s/\<\(zcat//g;
$head =~ s/\)//g;
$head =~ s/paste //g;
@heads = split(/\s+/, $head);
foreach my $aa (@heads){  #crate list of unique
    if($aa eq ""){next;}
    $aa =~ s/_mappedAGPv5.sorted.MarkedDup.bam./&/g;
    $aa =~ s/\.vcf.gz.tmpt.gz//g;
    #print "<$aa>\n";<STDIN>;
    @xx = split(/&/, $aa);
    print OUT "$xx[0]\t$xx[1]\n";
    if(!exists $allhd{$xx[0]}){push(@uheads,$xx[0]);}
    $allhd{$xx[0]}=0;
}
foreach my $aa (@uheads){
    print OUTA "$aa\n";
}
print "Done...\n";
