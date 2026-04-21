$file="/groups/swarts/lab/MAVE/MAIZEDB/ID2.lst";

$grp="I1";
$fxx=0;
print "Loading list...\n";
open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    chop($line);
    @col=split(/\t/, $line);
    $hl{$col[0]}=$col[1];

    if ($col[1] eq $grp){$fxx++;}
   # print "$col[0]  $hl{$col[0]} ";<STDIN>;
}
close(FIL);
print "Done!!!\n";

print "$grp = $fxx\n";


print "Reading file...\n";
$file="/scratch-cbe/users/miguel.vallebueno/DB/1.ALL.Merge.SAMP5.poli.qual0.Union.V8.max.production.vcf_IBS.dist.txt";

$out="$file.DHET.cuantex";


open(FIL,"$file") or die;
open(OUT,">$out") or die;

$c=1;
while($line=<FIL>) {
    chomp($line);
    if ($line =~ /^#/) {next;}
    if ($line =~ /^485/) {next;}
    @col=split(/\t/, $line);
    $hnames[$c]=$col[0];
    $c++;
}
close (FIL);


open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    if ($line =~ /^#/) {next;}
    if ($line =~ /^485/) {next;}
    @col=split(/\t/, $line);
    #print "$col[0]\n"; <STDIN>;
    $c=0;
      foreach my $aa (@col){
#print "<$col[0]>  <$hl{$col[0]}>   <$aa>   <$c>  <$hnames[$c]>  <$hl{$hnames[$c]}> ";<STDIN>;

          $key="$hnames[$c]&$col[0]";

     if($hl{$col[0]} eq "I1" and $hl{$hnames[$c]} eq "I1" and $col[0] ne $hnames[$c] and !exists $hobs{$key}){$hf{$key}=$aa;}#print OUT "$aa\n";}

          $key="$col[0]&$hnames[$c]";
          $hobs{$key}=1;
          $key="$hnames[$c]&$col[0]";
          $hobs{$key}=1;

          $c++;
      }
}
close (FIL);


##### Random sampler dist

$sam[0]=2;
$sam[1]=4;
$sam[2]=8;
$sam[3]=16;
$sam[4]=32;
$sam[5]=64;
$sam[6]=96;
$sam[7]=112;

$reps=100;

for ($i=1;$i<=$reps;$i++){

foreach my $aa (@sam){

    undef(@cun);
    for ($j=1;$j<=$aa;$j++){

        my $random_value = $hf{(keys %hf)[rand keys %hf]};

        push(@cun, $random_value);
        #print "<$i> <$j> <$aa> <$random_value>  <@cun>\n"; <STDIN>;

    }

    $mean=mean(@cun);

print OUT "$aa\t$mean\n";

    }
}


#https://stackoverflow.com/questions/8547642/selecting-a-random-key-from-a-hash

print "Done!!!\n";
sub mean {
    my @vals = @_;
    my $result;
    if (scalar @vals == 0){return 0}
    foreach (@vals) { $result += $_ }
    return $result / @vals;
}

close(OUT);