$PTH="/path/to/gvcf_inputs";
$key="files.mergein_Discovery_2023_05.lst";

$ct=0;
$cc=0;
open(LST,"$PTH/$key") or die;
open(OUT,">$PTH/$key.rv") or die;
while($line=<LST>) {
    chomp($line);
    print "<$line.ll>\n";
    $ct++;
    open(FIL,"$PTH/$line.ll");
    while($lin=<FIL>) {
        chomp($lin);

        @lcol = split(/\t/, $lin);

        if($lcol[0] ne $lcol[1]){next}
        if($lcol[0] eq ""){next}

        print "<$lcol[0]>  <$lcol[1]>\n";

        print OUT "$line\n";
$cc++;


    }


}

print "<$ct>  <$cc>";
