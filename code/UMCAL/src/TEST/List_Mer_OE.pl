#$file="1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp.2.IHET.cuantey"
$file=$ARGV[0];
open(OUT, ">$file.eo");
open(EXP, "/groups/swarts/lab/MAVE/Bioinfo/Diversity_panels/hmp321_chr1.crossmap.AGPV5.vcf.gz.HM.IHET.cuantex");
while($line = <EXP>){
    chomp($line);
    #chop($line);
    @lcol=split(/\t/,$line);
    $he{$lcol[0]}=$lcol[1];
}
close (FIL);
open(FIL, "$file");
while($line = <FIL>){
    chomp($line);
    #chop($line);
    @lcol=split(/\t/,$line);
    $n=$lcol[0];
    $n=~ s/SAMN\d+_RI_DeNovo_//;
    $n=~ s/SAMN\d+_NAM_DeNovo_//;
    $n=~ s/SAMEA\d+_NAM_DeNovo_//;
    $n=~ s/SAMN\d+_HapMap2_//;
    $n=~ s/SAMN\d+_HapMap3_//;
    $n=~ s/SAMN\d+_282_Maize_panel_//;
    if(exists $he{$n}){print OUT "$n\t$lcol[0]\t$lcol[1]\t$he{$n}\n";}
}
close (FIL);
close(OUT);