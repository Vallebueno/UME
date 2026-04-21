# UME Union Discovery
# Author: Miguel Vallebueno
# Date: 2022-07-29
$filen=$ARGV[0];
$chr=$ARGV[1];
$file="$filen";

$alg=$ARGV[2];
$lexp="$ARGV[3].pval2f";
$list="$ARGV[3].Disc.lst";
print "<$ARGV[0]><$ARGV[1]><$ARGV[2]><$ARGV[3].pval2f><$ARGV[4]>\n";


#$filen=$ARGV[0];
#$chr=$ARGV[1];
#$filen=~s/.gz//;
#$file="$filen.$chr.poli.db.gz";
#$alg=$ARGV[2];
#$lexp="$filen.pval2f";
#$list="$ARGV[3].Disc.lst";



$qualf=$ARGV[4]; ### quality filter in Phred score  WARNING::: should remain in 0 under current branch code
print "<$ARGV[0]><$ARGV[1]><$ARGV[2]><$ARGV[3].pval2f><$ARGV[4]>\n";

open(TST,"$lexp") or die;
close(TST);
open(TST,"$list") or die;
close(TST);
open(TST,"$file") or die;
close(TST);

#### posible to change but not recomended just for testing porpuses
$qualf=0; ### quality filter in Phred score  WARNING::: should remain in 0 under current branch code


$out="$file.qual$qualf.Union.$alg.dbx";
$out =~ s/.db.gz//g;

open (OUT, "> $out") or die;
open (OUU, "> $out.q.lst") or die;

#open (my $gout_fh, "| /bin/gzip -c > $out.gz") or die;
#open (my $gouu_fh, "| /bin/gzip -c > $out.q.lst.gz") or die;

print OUT "#0=missing info\n";
print OUT "#1=reference allele\n";
print OUT "#Mapped against B73-NAM AGPV5\n";

############################## start

#$hh{$pos}="$geno:$qual:$hformat{DP}:$hformat{AD}:$hformat{PL}";
$ccvc=0;
$ccvt=0;


print "Extracting headers from input list...\n";
open(LST,"$list") or die;
while($line=<LST>) {
    chomp($line);
    @llo= split(/\t/, $line);
    push(@heads,$llo[0]);
    push(@hcal,$llo[1]);
}
$lenh=scalar(@heads);
$lenh=$lenh+3;
foreach my $aa (@heads){  #crate list of unique
    if(!exists $allhd{$aa}){push(@finh,$aa);}
    $allhd{$aa}=0;
}
close(LST);
print "Done!\n";


###norm lookup table
print "Loading norm lookup table...\n";
open(EXP,"$lexp") or die;
while($line=<EXP>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $hpv{$lcol[0]}{$lcol[1]}=$lcol[2];
    push(@{$HAVA{$lcol[0]}}, $lcol[1]);
    push(@{$HAVB{$lcol[0]}}, $lcol[2]);
}
close (EXP);
print "Done!\n";

if ($alg eq "sum") {
    print "Loading norm SUM lookup table...\n";
    $fsum = $ENV{'UME_SUM_PVAL_TABLE'} || "";
    if ($fsum eq "") {
        die "UME_SUM_PVAL_TABLE must be set when using the sum aggregation mode\n";
    }
    open(EXP, "$fsum") or die;
    while ($line = <EXP>) {
        chomp($line);
        @lcol = split(/\t/, $line);
        $spv{$lcol[0]}{$lcol[1]} = $lcol[2];
        push(@{$SAVA{$lcol[0]}}, $lcol[1]);
        push(@{$SAVB{$lcol[0]}}, $lcol[2]);
    }
    close(EXP);
    print "Done!\n";
}

$sys=join("\t",@finh);

print OUT "#chr\tpos\tref\t$sys\n";

print "Opening file <$file>...\n";
#open(DB,"$file") or die;
open(DB, "gunzip -c $file |") or die "gunzip $file: $!";
print "Union...\n";
$aas=0;

$cw=0;
while($line=<DB>) {
    $aas++;
    chomp($line);
    if($line eq ""){next;}
    @lcol = split(/\t/, $line);
    $chr=$lcol[0];
    $pos=$lcol[1];
    $ref=$lcol[2];
    $len=scalar(@lcol);
    if($cw > 100){print "LETAL WARN: maximum number of warnings reached.... there is something wrong with the INPUT files." and die}
    if($len != $lenh){print "WARNING!!:38945939  length of header <$lenh> do not match the line length <$len>\n";$cw++;next}
    undef(%HoA);
    undef(%HoB);
    for ($i=3;$i<$len;$i++){
        $h=$i-3;
        $samp=$heads[$h];
        $call=$hcal[$h];
        $val=$lcol[$i];
        push(@{$HoA{$samp}}, $val);
        push(@{$HoB{$samp}}, $call);
    }
    $cct=0;
    foreach my $aa (@finh){
        @check=@{$HoA{$aa}};
        @acall=@{$HoB{$aa}};
        $sys=join("",@check);
        undef (%h);
        %h=();
        %hq=();
        %hd=();
        %hi=();
        %posq=();
        $ctr=0;
        foreach my $va (@check){
            $caller=$acall[$ctr];
            @ll = split(/:/, $va);
            $qual=sprintf("%.1f", $ll[1]);
            if($ll[2]>$hd{$ll[0]}){$hd{$ll[0]}=$ll[2];$hi{$ll[0]}=$ll[3];}
            if($qual eq "." or $qual == 0){if($ll[0] eq "$ref$ref"){$h{$ll[0]}=$h{$ll[0]}+1} ;$ctr++;  next;}
            if(exists $hpv{$caller}{$qual}){$norm=$hpv{$caller}{$qual};
            }
            else{
                $bin = bubble_search(\@{$HAVA{$caller}},$qual);
                $norm=@{$HAVB{$caller}}[$bin];
                $cht=@{$HAVA{$caller}}[$bin];
                $hpv{$caller}{$qual}=$norm;####add value to the pnorm so does not look again for val
            }
            $ctr++;
            $h{$ll[0]}=$h{$ll[0]}+1;
            push(@{$posq{$ll[0]}}, $norm);
        }
        $funfg=0;
        foreach my $qx (keys %h){
            @sarray = sort { $a <=> $b } @{$posq{$qx}};
            if($alg eq "max"){$max=$sarray[-1];$hq{$qx}=$max;}
            if($alg eq "min"){$min=$sarray[0];$hq{$qx}=$min;}
            if($alg eq "mean"){$mean=mean(@sarray);$hq{$qx}=$mean;}
            if($alg eq "median"){$median=median(@sarray);$hq{$qx}=$median;}
            if($alg eq "sum"){
                $sumx=sum(@sarray);
                $sum=sprintf("%.1f", $sumx);
                if(exists $spv{"SUM"}{$sum}){$snorm=$spv{"SUM"}{$sum};
                }
                else{
                    $bin = bubble_search(\@{$SAVA{"SUM"}},$sum);
                    $snorm=@{$SAVB{"SUM"}}[$bin];
                    $cht=@{$SAVA{"SUM"}}[$bin];
                    $spv{"SUM"}{$sum}=$snorm;####add value to the pnorm so does not look again for val
                }
                $hq{$qx}=$snorm;
            }
            $a=$hq{$qx};
            $hq{$qx}=n2phred($a);
        }
        $f=0;
        $fn=0;
        $fx=0;
        $chs=0;
        $gan="";
        $rfin="";
        $win=0;
        $ccvt++;
  $ffof=0;
        $hwq=0;
        foreach my $bb (keys %h){
            if($bb eq "$ref$ref"){$f=1;$rfin="$bb:$hq{$bb}:$hd{$bb}:$hi{$bb}";next;}
            if($bb eq ""){$fn=1;next}
            if($h{$bb}>$chs ){$fx=1;$chs=$h{$bb};$gan="$bb:$hq{$bb}:$hd{$bb}:$hi{$bb}";$hwq=$hq{$bb};$hqx=$hq{$bb}}
            if( $h{$bb} == $chs and $hq{$bb} >= $hwq){$fx=1;$chs=$h{$bb};$gan="$bb:$hq{$bb}:$hd{$bb}:$hi{$bb}";$hwq=$hq{$bb};$hqx=$hq{$bb}}
        }
       if($fx==1){$win=$gan; print OUU "$hqx\n";}
        elsif($f==1){$win=$rfin;}
        elsif($fn==1){$win=0;}
        push(@winers, $win);
    }
    $sys=join("\t",@winers);
    print OUT "$chr\t$pos\t$ref\t$sys\n";
    undef (@winers);
}
close(OUT);
close(OUU);
close(DB);

print "Done!\n";
#close($goou_fh);
#https://perlmaven.com/binary-search-in-perl-array
#my $res = binary_search($name, \@planets);
#https://www.perlmonks.org/?node_id=260548
sub bubble_search {
    my @array = @{$_[0]};
    my $key=$_[1];
    my $low = 0;
    $high = scalar @array -1;
    while($low < $high){
        $mid = int(( $low + $high ) / 2);
        if($low==$high or $low == $high-1){$index=$low;last}
        if($key==$array[$mid]){$index=$mid; last;}
        elsif($key>$array[$mid]){$low=$mid;}
        elsif($key<$array[$mid]){$high=$mid;}
    }
    return $index;
}
sub median {
    my @vals = sort {$a <=> $b} @_;
    my $len = @vals;
    if($len%2) #odd?
    {
        return $vals[int($len/2)];
    }
    else #even
    {
        return ($vals[int($len/2)-1] + $vals[int($len/2)])/2;
    }
}
sub mean {
    my @vals = @_;
    my $result;
    if (scalar @vals == 0){return 0}
    foreach (@vals) { $result += $_ }
    return $result / @vals;
}
sub sum {
   my @vals = @_;
    my $result;
    foreach my $i (@vals){
        $result = $result + (1-$i);
    }
return $result;
}
sub n2phred{
    my $n = $_[0];
    $res=-10*log10($n);
    return $res;
}
sub log10 {
    my $n = shift;
    if($n == 0){return "."}
    return log($n)/log(10);
}
