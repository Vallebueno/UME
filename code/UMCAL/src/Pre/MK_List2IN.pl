$lst=$ARGV[0];

open(OUT,">$lst.Disc.lst") or die;
open(LST,"$lst") or die;
while($line=<LST>) {
    @lcol = split(/\./, $line);
    $name=$lcol[0];
    $caller=$lcol[4];
    $name=~ s/_mappedAGPv5//;
    print OUT "$name\t$caller\n";
}
