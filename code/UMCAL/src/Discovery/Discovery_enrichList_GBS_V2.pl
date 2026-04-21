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

$GBSf="/groups/swarts/lab/MAVE/Bioinfo/Diversity_panels/GBS/Maize_GBS_Kriztian_build.lst";

$out=$file;
$out =~ s/.lst2/.lst3/g;
$chr=10;
open(OUT,">$out") or die;
print "Loading GBS sites...\n";
open(GBS,"$GBSf") or die;
while($line=<GBS>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    if($lcol[0] ne $chr){next;}
    $key="$lcol[0]&$lcol[1]";
    $sys=$lcol[2];
    $hgbs{$key}=$sys;
    $all{$key}=0;
}
close(GBS);

print "Loading Discovery List sites...\n";
open(FIL,"$file") or die;

while($line=<FIL>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~s/ //g;
    $key="$lcol[0]&$lcol[1]";
    $hdb{$key}=$lcol[2];
    $all{$key}=0;
}
close (FIL);
print "Unifiying sites...\n";
foreach my $aa (sort { $a <=> $b } keys %all){
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