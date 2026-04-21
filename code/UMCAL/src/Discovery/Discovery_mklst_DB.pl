# UME Make list of known positions based on the assigned cutoff
# Author: Miguel Vallebueno
# Date: 2022-07-29

$filen=$ARGV[0];
$chr=$ARGV[1];
$cof=$ARGV[2];
$alg=$ARGV[3];
$qualf=$ARGV[4];








$file="$filen.$chr.poli.qual$qualf.Union.$alg.dbx";
$file =~ s/.db.gz//g;
#print "<$file>\n";
$OUT="$file.cof$cof.lst";

open(OUT,">$OUT") or die;
#open(DB,"gunzip -c $file |") or die;
open(DB,"$file") or die;
while($line=<DB>) {
    chomp($line);
    if($line =~ /^#chr/){
        @heads = split(/\t/, $line);
    }
    if($line =~ /^#/){next;}
    @lcol = split(/\t/, $line);
    $len=scalar(@lcol);
    $chr=$lcol[0];
    $pos=$lcol[1];
    $A=0;
    $C=0;
    $T=0;
    $G=0;
    undef(@alt);
    for ($i=2;$i<$len;$i++) {
        $samp = $heads[$i];
        $val = $lcol[$i];
        if($i==2){
                if($val=~/A/){$A++}
                if($val=~/C/){$C++}
                if($val=~/T/){$T++}
                if($val=~/G/){$G++}
        }
        @info = split(/:/, $val);
        if($info[1] < $cof){next;}
        if($val=~/A/){$A++}
        if($val=~/C/){$C++}
        if($val=~/T/){$T++}
        if($val=~/G/){$G++}
    }
    if($A>0){push(@alt,"A")}
    if($C>0){push(@alt,"C")}
    if($T>0){push(@alt,"T")}
    if($G>0){push(@alt,"G")}
    if(@alt<2){next}
    print OUT "$chr\t$pos\t@alt\n";
}
close(OUT);