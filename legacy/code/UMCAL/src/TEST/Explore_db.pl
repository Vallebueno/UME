$file="/path/to/maizedb/ID.lst";

$cofh=$ARGV[0];
$prp=10;
$sf=0;
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

$out="$file.$cofh.cuantex";
$outa="$file.$cofh.mins.cuantex";

#open(FIL,"gunzip -c $file |") or die;
open(FIL,"$file") or die;
open(OUT,">$out") or die;
open(OUTA,">$outa") or die;
$ccc=0;
while($line=<FIL>) {
    $ccc++;
chomp($line);
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
    $fhi=0;
    for($i=9;$i<$tam;$i++){
        if ($col[$i] =~ m/1\/.*:.*/) {$ala++}
        if ($col[$i] =~ m/.*\/1:.*/) {$ala++}
        if ($col[$i] =~ m/2\/.*:.*/) {$alb++}
        if ($col[$i] =~ m/.*\/2:.*/) {$alb++}
        if ($col[$i] =~ m/3\/.*:.*/) {$alc++}
        if ($col[$i] =~ m/.*\/3:.*/) {$alc++}
        if ($col[$i] =~ m/4\/.*:.*/) {$ald++}
        if ($col[$i] =~ m/.*\/4:.*/) {$ald++}

        #if ($col[$i] =~ m/.*1.*:.*/) {$ala++}
        #if ($col[$i] =~ m/.*2.*:.*/) {$alb++}
        #if ($col[$i] =~ m/.*3.*:.*/) {$alc++}
        #if ($col[$i] =~ m/.*4.*:.*/) {$ald++}
        if($hl{$hn{$i}} eq "I1" and $col[$i] !~ /\./  and $col[$i] !~ /0\/0/ and $col[$i] !~ /1\/1/ and $col[$i] !~ /2\/2/ and $col[$i] !~ /3\/3/ and $col[$i] !~ /4\/4/ ){
            $fhi++;
        }
    }
    if($fhi>=$cofh){$sfh++;next;}
    if($ala>1 or $alb>1 or $alc>1 or $ald>1 ){$sf++;}
    my $rs = int(rand(100));  ### random sampler
    if ($rs > $prp) {next;}
undef(@x);
    if($ala ne "0"){ push(@x, $ala)}
    if($alb ne "0"){ push(@x, $alb)}
    if($alc ne "0"){ push(@x, $alc)}
    if($ald ne "0"){ push(@x, $ald)}
        my ($min,$max) = (sort {$a <=> $b} @x)[0,-1];
    if($min==""){next}
    $minf=$min/(($tam-9)*2);
    print OUTA "$minf\n";
    #print  "<$min>\n";
}
$asfh=$sfh/$ccc;
$asf=$sf/$ccc;
print OUT "$asfh\t$asf\n";
close(OUT);
close(FIL);
close(OUTA);
print "Done!!!\n";
