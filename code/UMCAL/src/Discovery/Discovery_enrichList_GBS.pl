# UME makes the union of discovery DB and GBS DB (site summary file)
# Author: Miguel Vallebueno
# Date: 2022-07-29

#$filen=$ARGV[0];
#$chr=$ARGV[1];
#$cof=$ARGV[2];
#$alg=$ARGV[3];
#$qualf=$ARGV[4];
#$file="$filen.$chr.poli.qual$qualf.Union.V8.$alg.dbx.cof$cof.lst";
#$file =~ s/.db.gz//g;
#$file = $filen;
$file=$ARGV[0];
$chr=$ARGV[1];
$GBSf=$ENV{'UME_GBS_SUMMARY'} || "/groups/swarts/lab/MAVE/MAIZEDB/GBS_Kriztinan_filt_summary.txt";

$out=$file;
$out =~ s/.lst/.lst2/g;

open(OUT,">$out") or die;

open(GBS,"$GBSf") or die;
while($line=<GBS>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    if($lcol[2] ne $chr){next;}
    $key="$lcol[2]&$lcol[3]";
    #print "<$lcol[2]> <$lcol[3]>  <$lcol[5]> <$lcol[7]>  <$lcol[11]> <$lcol[15]>  <$lcol[19]> "  ; <STDIN>;
    push(@fin, $lcol[5]);
    $h{$lcol[5]}=0;
    if($lcol[6] ne "N" and $lcol[6] ne "N" and $lcol[6] ne "NA" and !exists $h{$lcol[6]}){push(@fin, $lcol[6]);$h{$lcol[6]}=0;}
    if($lcol[7] ne "N" and $lcol[7] ne "NA" and !exists $h{$lcol[7]}){push(@fin, $lcol[7]);$h{$lcol[7]}=0;}
    if($lcol[11] ne "N" and $lcol[11] ne "NA" and !exists $h{$lcol[11]}){push(@fin, $lcol[11]);$h{$lcol[11]}=0;}
    if($lcol[15] ne "N" and $lcol[15] ne "NA" and !exists $h{$lcol[15]}){push(@fin, $lcol[15]);$h{$lcol[15]}=0;}
    if($lcol[19] ne "N" and $lcol[19] ne "NA" and !exists $h{$lcol[19]}){push(@fin, $lcol[19]);$h{$lcol[19]}=0;}
    $sys=join(" ", @fin);
    $hgbs{$key}=$sys;
    $all{$key}=0;
    undef(%h);
    undef(@fin);
}
close(GBS);
open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $key="$lcol[0]&$lcol[1]";
    $hdb{$key}=$lcol[2];
    $all{$key}=0;
}
close (FIL);

foreach my $aa (sort keys %all){
    @lcol = split(/&/, $aa);
    @adb = split(/ /, $hdb{$aa});
    @agbs = split(/ /, $hgbs{$aa});
    undef(%uniques);
    @uniques{@adb} = @adb x (1);
    @uniques{@agbs} = @agbs x (1);
    @adb = sort keys %uniques;
    print OUT "$lcol[0]\t$lcol[1]\t@adb\n";
}
close(OUT);
