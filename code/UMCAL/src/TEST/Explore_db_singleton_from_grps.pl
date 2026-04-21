$file="/groups/swarts/lab/MAVE/MAIZEDB/ID.lst";

$grp=$ARGV[0];
$prp=10;

print "Loading list...\n";
open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    chop($line);
    @col=split(/\t/, $line);
    $hl{$col[0]}=$col[1];
    $thl{$col[1]}=$thl{$col[1]}+1;
}
close(FIL);
print "Done!!!\n";
print "Reading file...\n";
$file="/scratch-cbe/users/miguel.vallebueno/DB/1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp";


$outa="$file.$grp.sgrp.cuantex";

#open(FIL,"gunzip -c $file |") or die;
open(FIL,"$file") or die;

open(OUTA,">$outa") or die;
$ccc=0;
$sing=0;
while($line=<FIL>) {

    chomp($line);
    if($line =~ /^#C/){
        @col=split(/\t/, $line);
        $tam=scalar(@col);
        for($i=9;$i<$tam;$i++){$hn{$i}=$col[$i];}
    }
    if($line =~ /^#/){next;}

    $ccc++;
    @col=split(/\t/, $line);
    $ala=0;
    $alb=0;
    $alc=0;
    $ald=0;
    $fx=0;

    undef(%als);

    for($i=9;$i<$tam;$i++){
        if($hl{$hn{$i}} ne $grp){next;}
        if ($col[$i] =~ m/.*1.*:.*/) {$als{1}=1;$fx=1}
        if ($col[$i] =~ m/.*2.*:.*/) {$als{2}=1;$fx=1}
        if ($col[$i] =~ m/.*3.*:.*/) {$als{3}=1;$fx=1}
        if ($col[$i] =~ m/.*4.*:.*/) {$als{4}=1;$fx=1}
      }

    if($fx==0){next;}

    for($i=9;$i<$tam;$i++){
        if($hl{$hn{$i}} eq $grp){next;}

        if ($col[$i] =~ m/1\/.*:.*/ and exists $als{1}) {$ala++}
        if ($col[$i] =~ m/.*\/1:.*/ and exists $als{1}) {$ala++}
        if ($col[$i] =~ m/2\/.*:.*/ and exists $als{2}) {$alb++}
        if ($col[$i] =~ m/.*\/2:.*/ and exists $als{2}) {$alb++}
        if ($col[$i] =~ m/3\/.*:.*/ and exists $als{3}) {$alc++}
        if ($col[$i] =~ m/.*\/3:.*/ and exists $als{3}) {$alc++}
        if ($col[$i] =~ m/4\/.*:.*/ and exists $als{4}) {$ald++}
        if ($col[$i] =~ m/.*\/4:.*/ and exists $als{4}) {$ald++}
       }
    if($ala > 0 or $alb > 0 or $alc > 0 or $ald > 0){next;}
    $sing++;
}

$por=$sing/$ccc;

print OUTA "$sing\t$ccc\t$por\n";

close(FIL);
close(OUTA);
print "Done!!!\n";