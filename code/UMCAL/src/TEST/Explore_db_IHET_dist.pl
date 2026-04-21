$file="/groups/swarts/lab/MAVE/MAIZEDB/ID2.lst";

$grp="I1";

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
$file="/scratch-cbe/users/miguel.vallebueno/DB/1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf";


$out="$file.$cofh.DHET.cuantex";


#open(FIL,"gunzip -c $file |") or die;
open(FIL,"$file") or die;
open(OUT,">$out") or die;

$ccc=0;
while($line=<FIL>) {

chomp($line);
  #  print "$ccc\n";
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
        if($hl{$hn{$i}} ne $grp){next;}

        if ($col[$i] =~ m/.*1.*:.*/) {$ala++}
        if ($col[$i] =~ m/.*2.*:.*/) {$alb++}
        if ($col[$i] =~ m/.*3.*:.*/) {$alc++}
        if ($col[$i] =~ m/.*4.*:.*/) {$ald++}

        if($hl{$hn{$i}} eq "I1" and $col[$i] !~ /\./  and $col[$i] !~ /0\/0/ and $col[$i] !~ /1\/1/ and $col[$i] !~ /2\/2/ and $col[$i] !~ /3\/3/ and $col[$i] !~ /4\/4/ ){
            $fhi{$hn{$i}}=$col[$i];
            $f++;
            #print "yyy $hn{$i}  $col[$i]  <<$f>>\n"; <STDIN>;
        }
    }

    undef(%ho);
    undef(@cun);

    $ccx=0;

    foreach my $aa (keys %fhi){
        foreach my $bb (keys %fhi){
            $sumx=0;
            $key="$aa&$bb";
            if(exists $ho{$key}){next;}
            if($aa eq $bb){next;}
            @acol=split(/:/, $fhi{$aa});
            @bcol=split(/:/, $fhi{$bb});
            @acl=split(/\//, $acol[0]);
            @bcl=split(/\//, $bcol[0]);
            if($acl[0] ne $bcl[0] and $acl[0] ne $bcl[1]){$sumx=$sumx+0.5;}
            if($acl[1] ne $bcl[0] and $acl[1] ne $bcl[1]){$sumx=$sumx+0.5;}
      #      if($acl[0] eq $bcl[0] or $acl[0] eq $bcl[1] and $acl[1] eq $bcl[0] or $acl[1] eq $bcl[1]){$ccx++;}
            push(@cun, $sumx);
          # print "xxx $aa $bb $acol[0] $bcol[0]   <@cun> <<$f>>\n";<STDIN>;
            $key="$bb&$aa";
            $ho{$key}=1;
        }
    }

    $mean=mean(@cun);
   # print "XXXX<@cun>   <<$mean>>\n";<STDIN>;
    print OUT "$f\t$mean\n";
}

close(OUT);
close(FIL);

print "Done!!!\n";
sub mean {
    my @vals = @_;
    my $result;
    if (scalar @vals == 0){return 0}
    foreach (@vals) { $result += $_ }
    return $result / @vals;
}