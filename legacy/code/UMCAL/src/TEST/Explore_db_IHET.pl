$file="/path/to/maizedb/ID.lst";

$cofh=$ARGV[0];


print "Loading list...\n";
open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    chop($line);
@col=split(/\t/, $line);

    $hl{$col[0]}=$col[1];
}
close(FIL);
print "Done!!!\n";
print "Reading file...\n";
$file="/path/to/workdir/DB/1.ALL.Merge.ALLC.poli.qual0.Union.V8.max.production.vcf.gz.smp";

$out="$file.$cofh.IHET.cuantey";


#open(FIL,"gunzip -c $file |") or die;
open(FIL,"$file") or die;
open(OUT,">$out") or die;

$ccc=0;
while($line=<FIL>) {

chomp($line);
    print "$ccc\n";
    if($line =~ /^#C/){
        @col=split(/\t/, $line);
        $tam=scalar(@col);
        for($i=9;$i<$tam;$i++){$hn{$i}=$col[$i];}
    }
    if($line =~ /^#/){next;}

    @col=split(/\t/, $line);
    $ala=0;
    $alb=0;
    $alc=0;
    $ald=0;
    $f=0;
undef(%fhi);
    for($i=9;$i<$tam;$i++){
        if ($col[$i] =~ m/.*1.*:.*/) {$ala++}
        if ($col[$i] =~ m/.*2.*:.*/) {$alb++}
        if ($col[$i] =~ m/.*3.*:.*/) {$alc++}
        if ($col[$i] =~ m/.*4.*:.*/) {$ald++}

  #      if ($col[$i] =~ m/1\/.*:.*/) {$ala++}
  #      if ($col[$i] =~ m/.*\/1:.*/) {$ala++}
  #      if ($col[$i] =~ m/2\/.*:.*/) {$alb++}
  #      if ($col[$i] =~ m/.*\/2:.*/) {$alb++}
  #      if ($col[$i] =~ m/3\/.*:.*/) {$alc++}
  #      if ($col[$i] =~ m/.*\/3:.*/) {$alc++}
  #      if ($col[$i] =~ m/4\/.*:.*/) {$ald++}
  #      if ($col[$i] =~ m/.*\/4:.*/) {$ald++}



        if($hl{$hn{$i}} eq "I1" and $col[$i] !~ /\./  and $col[$i] !~ /0\/0/ and $col[$i] !~ /1\/1/ and $col[$i] !~ /2\/2/ and $col[$i] !~ /3\/3/ and $col[$i] !~ /4\/4/ ){
            $fhi{$hn{$i}}=1;
            $f++;
            print "yyy $hn{$i}\n";
        }
    }
    if($f>=$cofh){next;}
    if($ala<1 and $alb<1 and $alc<1 and $ald<1 ){next;}
    $ccc++;

    foreach my $aa (keys %fhi){
        $fhic{$aa}=$fhic{$aa}+1;
    }
}

foreach my $aa (keys %fhic){
    $het=$fhic{$aa}/$ccc;
    print OUT "$aa\t$het\n";
}

close(OUT);
close(FIL);

print "Done!!!\n";
