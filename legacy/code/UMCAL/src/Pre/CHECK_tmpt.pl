$PTH="/path/to/workdir/MERGE";
$key="files.mergein_Discovery_2023_05.lst";

$ct=0;
$cc=0;
open(LST,"$PTH/$key") or die;
open(OUT,">$PTH/$key.rv2") or die;
while($line=<LST>) {
    chomp($line);

    $ct++;
    $filename="$PTH/$line.tmpt.gz";

    if (-e $filename) {next}
    print "$PTH/$line.tmpt.gz"
}

print "<$ct>  <$cc>";
