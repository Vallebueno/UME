$file="/path/to/diversity_panels/hmp321_chr1.crossmap.AGPV5.vcf.gz";
$ref="/path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL";

print "Reference:<$ref>\n";
open(REF,"$ref") or die;

print "loading ref...\n";
chomp(my @re = <REF>); #magic
$f=0;
foreach my $aa(@re){
    if($f==1){$href{$name}=$aa;$f=0;next;}
    if($f==0){$name=$aa; $name=~ s/>//; print "$name\n";$f=1;$ctrl{$name}=1;next;}
}
close(REF);
undef(@re);
print "Done!...\n";



$GBSf="/path/to/maizedb/GBS_Kriztinan_filt_summary.txt";
$chr=1;
open(GBS,"$GBSf") or die;
while($line=<GBS>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    if($lcol[2] ne $chr){next;}
    $key="$lcol[3]";
    # print "<$lcol[2]> <$lcol[3]>  <$lcol[5]> <$lcol[7]>  <$lcol[11]> <$lcol[15]>  <$lcol[19]> "  ; <STDIN>;
    push(@fin, $lcol[5]);
    $h{$lcol[5]}=0;
    if($lcol[6] ne "N" and $lcol[6] ne "N" and $lcol[6] ne "NA" and !exists $h{$lcol[6]}){push(@fin, $lcol[6]);$h{$lcol[6]}=0;}
    if($lcol[7] ne "N" and $lcol[7] ne "NA" and !exists $h{$lcol[7]}){push(@fin, $lcol[7]);$h{$lcol[7]}=0;}
    if($lcol[11] ne "N" and $lcol[11] ne "NA" and !exists $h{$lcol[11]}){push(@fin, $lcol[11]);$h{$lcol[11]}=0;}
    if($lcol[15] ne "N" and $lcol[15] ne "NA" and !exists $h{$lcol[15]}){push(@fin, $lcol[15]);$h{$lcol[15]}=0;}
    if($lcol[19] ne "N" and $lcol[19] ne "NA" and !exists $h{$lcol[19]}){push(@fin, $lcol[19]);$h{$lcol[19]}=0;}
    $sys=join("", @fin);
    $sys =~ s/\-//;
    $sys =~ s/\+//;
    $hgbs{$key}=$sys;
    undef(%h);
    undef(@fin);
}
close(GBS);

print "Done!...\n";




$ccg=0;
$cch=0;
$cctg=0;
$ccth=0;
$cct=0;
$cr=0;


#$gref=substr $href{$chr},$i-1, 1;

print "Loading HM list...\n";
open(FIL,"gunzip -c $file |") or die;
while($line=<FIL>) {
    if($line =~ /^#/){next;}
    @col = split(/\t/, $line);
    $col[4]=~ s/<DEL>//;
    $col[4]=~ s/<INS>//;
   # print " <$col[0]> <$col[1]>   <$col[3]> <$col[4]>";<STDIN>;
    $hm{$col[1]}="$col[3]$col[4]";
    $hm{$col[1]}=~ s/,//g;
}
close(FIL);

print "Done!!!\n";

print "Loading Disc list...\n";

$file="1.ALL.Merge.1.poli.qual0.Union.V8.max.dbx.cof5.lst";
open(FIL,"$file") or die;
open(OUT,">$file.ccnLHG") or die;
while($line=<FIL>) {
    if ($line =~ /^#/) {next;}
    chomp($line);
    @col = split(/\t/, $line);

$cct++;
    if(exists $hgbs{$col[1]}){
        $cctg++;
        @cyy = split(//, $hgbs{$col[1]});
        $gref=substr $href{1},$col[1]-1, 1;
        #if($col[2] =~ /$cyy[0]/ and $col[2] =~ /$cyy[1]/){$ccg++;}
        if($col[2] !~ /$cyy[0]/ or $col[2] !~ /$cyy[1]/){$ccg++;}


    }### end GBS


    if(!exists $hm{$col[1]}){$nhm++;next;}
    $ccth++;
    Liubliana, Eslovenia
    @cxx = split(//, $hm{$col[1]});

    $gref=substr $href{1},$col[1]-1, 1;

    $col[2]=~ s/ //g;

    if(exists $hm{$col[1]}){
     #   if($col[2] =~ /$cxx[0]/ and $col[2] =~ /$cxx[1]/){$cch++;}
        if($col[2] !~ /$cxx[0]/ or $col[2] !~ /$cxx[1]/){$cch++;}
        if($col[2] =~ /$gref/) {$cr++;}
    }
}
print "Done!!!\n";
close(FIL);

$ph=$cch/$ccth;
$pg=$ccg/$cctg;
$pr=$cr/$ccth;

print OUT "$cct\t$cch\t$ccth\t$ph\t$ccg\t$cctg\t$pg\t$cr\t$pr\n";
close(OUT);
