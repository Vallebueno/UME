$file="/groups/swarts/lab/MAVE/MAIZEDB/ID.lst";

print "Loading list...\n";
open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    chop($line);
@col=split(/\t/, $line);


    if($col[1] ne "I1"){next;}

    $col[0] =~ s/SAMN\d+_//;
    $col[0] =~ s/SAMEA\d+_//;
    $col[0] =~ s/NAM_DeNovo_//;
    $col[0] =~ s/HapMap2_//;
    $col[0] =~ s/HapMap3_//;
    $col[0] =~ s/Maize_panel_//;
    $col[0] =~ s/RI_DeNovo_//;
    $col[0] =~ s/282/282set/;
    #print "$col[0] $col[1]";<STDIN>;

    $hl{$col[0]}=$col[1];

}

close(FIL);
print "Done!!!\n";
print "Reading file...\n";
$file="/groups/swarts/lab/MAVE/Bioinfo/Diversity_panels/hmp321_chr1.crossmap.AGPV5.vcf.gz";

$out="$file.HM.IHET.cuantex";


open(FIL,"gunzip -c $file |") or die;
#open(FIL,"$file") or die;
open(OUT,">$out") or die;

$ccc=0;
while($line=<FIL>) {
chomp($line);
    if($line =~ /^#CHROM/){
        @col=split(/\t/, $line);
        $tam=scalar(@col);
        for($i=0;$i<$tam;$i++){ if (exists $hl{$col[$i]}){$hn{$i}=$col[$i]}}
    }
    if($line =~ /^#/){next;}
$ccc++;
    @col=split(/\t/, $line);
    for($i=9;$i<$tam;$i++){
        if(!exists $hn{$i}){next;}
        if( $col[$i] =~ /\./  or $col[$i] =~ /0\/0/ or $col[$i] =~ /1\/1/ or $col[$i] =~ /2\/2/ or $col[$i] =~ /3\/3/ or $col[$i] =~ /4\/4/ ){next;}
            $fhi{$hn{$i}}=$fhi{$hn{$i}}+1;
        }
}

foreach my $aa (keys %fhi){
    $het=$fhi{$aa}/$ccc;
    print OUT "$aa\t$het\n";
}
close(OUT);
close(FIL);
print "Done!!!\n";