#### scripot to set to missing data positions that do not pass certain Alelic depth AD treshold
$file=$ARGV[0];
print "opening <$file>...\n";
open(FIL,"gunzip -c $file |") or die;
print "creating output <$file.s2m.gz>...\n";
open (my $gout_fh, "| /bin/gzip -c > $file.s2m.gz") or die;
$cof=8; # AD threshold
print "Reading tmpt...\n";
while($slin=<FIL>) {
    $sum=0;
    chomp($slin);
    if($slin eq "./."){print $gout_fh "./.\n";next}
    @lcol = split(/:/, $slin);
    @AD = split(/,/, $lcol[1]);
    foreach (@AD) { $sum += $_ }
    if($sum<$cof){print $gout_fh "./.\n";next}
    print $gout_fh "$lcol[0]\n";
}
close(FIL);
close($gout_fh);
