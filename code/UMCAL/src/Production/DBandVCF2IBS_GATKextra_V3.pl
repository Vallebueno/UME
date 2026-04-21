$list="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone.Merge.db.7.poli.qual0.Union.V7-1.max.db.lst";



$ccox=20;
open(LST,"$list") or die;
print "Loading sites...\n";
while($line=<LST>) {
    my $color = int(rand(100));  ### random sampler
    if ($color > $ccox) {next;}
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]";
    @lx = split(/ /, $lcol[2]);
    $sys=join("/", @lx);
    $hlst{$key}=$sys;
}
close(LST);
print "DONE!!!!!\n";

print "extracting snps from GATK control...\n";

$filex="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone.Merge.db.7.poli.db.gz";
$list="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone";

print "Extracting headers from input list...\n";
open(LST,"$list") or die;
chomp(my @LLL = <LST>);
@lcol = split(/\|/, $LLL[0]);
$head=$lcol[0];

#print "<$head>";<STDIN>;
$head =~ s/\.vcf.tmpt.gz//g;
$head =~ s/\.7.gvcf.tmpt.gz//g;
$head =~ s/\.bam.vcf.tmpt.gz//g;
$head =~ s/\.vcf.gz.tmpt.gz//g;
$head =~ s/\.bam//g;
$head =~ s/\.7//g;
$head =~ s/\.samtools//g;
$head =~ s/\.deepvar//g;

$head =~ s/freebayes//g;
$head =~ s/\<\(zcat//g;
$head =~ s/\)//g;
$head =~ s/\.vcf.gz.tmpt.gz//g;
$head =~ s/\.sorted.MarkedDup.bam//g;
$head =~ s/paste  //g;
$head =~ s/mappedAGPv5//g;
$head =~ s/_.I//g;
$head =~ s/x./x/g;
$head =~ s/> files.mergein3.lst.tlone.Merge.db//g;

@heads = split(/\s+/, $head);

$c=0;
foreach my $aa(@heads){  #crate list of unique
    if($aa=~ /30xgatk/){$hgatk{$c}=$aa; push(@allnames,$aa);}
    if($aa=~ /0.5xgatk/){$hgatk{$c}=$aa;push(@allnames,$aa);}
    $c++;
}

close(LST);
print "Done!\n";

print "opening gatks\n";
open(DB,"gunzip -c $filex |") or die;
while($line=<DB>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $chr = $lcol[0];
    $chr =~ s/chr//;
    $chr =~ s/ //;
    $pos = $lcol[1];

    $key="$chr&$pos";

    if(!exists $hlst{$key}){next;}

    $ref = $lcol[2];
    $len=scalar(@lcol);
    undef(%HoA);
    undef(%HoB);
    for ($i = 3; $i < $len; $i++) {
        $h = $i - 3;
        $samp = $hgatk{$h};
        $val = $lcol[$i];


        undef(@ran);
        @XA = split(/:/, $lcol[$i]);
        @GT = split(//, $XA[0]);
        push(@ran, $GT[0]);
        push(@ran, $GT[1]);
        $gen=$ran[rand @ran];
        if($gen eq ""){next;}

        if(exists $hgatk{$h}){$H{$key}{$samp}=$gen; } #print "<$h> <$key> <$samp> <$val><$gen>\n";<STDIN>

    }
}
close (DB);

print "Done!\n";

print "extracting snps from VCF production calling...\n";


$file="/scratch-cbe/users/miguel.vallebueno/Downsampling2/mpileup.lst.production.db.vcf";
$out="$file.hmp";


open(VCF, "$file") or die;
$r=0;
while($line=<VCF>) {
    if ($line =~ /^##/) {next;}
    chomp($line);
    @lcol = split(/\t/, $line);

    if ($line =~ /#chr/ | $line =~ /#CHR/) {
        $f = 0;
        if ($r == 1) {next;}
        $r = 1;
        $ccfs = 0;
        foreach my $aa (@lcol) {

                $ht{$aa} = $f;
                if($aa eq "#CHROM"){$f++;next;}
                if($aa eq "POS"){$f++;next;}
            if($aa eq "ID"){$f++;next;}
            if($aa eq "REF"){$f++;next;}
            if($aa eq "ALT"){$f++;next;}
            if($aa eq "QUAL"){$f++;next;}
            if($aa eq "FILTER"){$f++;next;}
            if($aa eq "INFO"){$f++;next;}
            if($aa eq "FORMAT"){$f++;next;}
            #print "<$aa>  <$f>\n"; <STDIN>;
                push(@at, $f);
                $ccfs++;
                $hDB{$f} = $aa;
                $inx{$aa} = 1;
              push(@allnames,$aa);
            $f++;
        }
        print "inds found in file that match list <$ccfs>\n";
        next;
    } #end header

    $chr=$lcol[0];
    $chr=~ s/chr//;
    $chr=~ s/ //;
    $pos=$lcol[1];
    $key="$chr&$pos";
    $ref=$lcol[3];
    if(!exists $hlst{$key}){next;}
    #### alts array
    undef(@alts);
    @vars=split(/,/,$lcol[4]);
    push(@alts,$ref);
    foreach my $ee (@vars){
        push(@alts,$ee);
    }
#print "<@alts>  $line\n"; <STDIN>;
    foreach my $aa (@at) {
        $namex=$hDB{$aa};
        undef(@ran);
        @info = split(/:/, $lcol[$aa]);
        if ($info[0] eq "./."){ next;}
        @GT = split(/\//, $info[0]);
        $geno="$alts[$GT[0]]$alts[$GT[1]]";
        push(@ran, $alts[$GT[0]]);
        push(@ran, $alts[$GT[1]]);
        $gen=$ran[rand @ran];
      #print "XXXXXXXXX pos<$key>  pind<$aa> <$namex> ran<@ran>  ind<$namx>  infoin<$lcol[$aa]>  GT1<$GT[0]> GT2<$GT[1]> geno<$geno> <$gen>   <$algx><$cutoff>\n"; <STDIN>;
        $H{$key}{$namex}=$gen;
    }
} #end VCF

open(OUT, ">$out") or die;
print "Done!!!\n";

$gg=join "\t", @allnames;

print OUT "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t$gg\n";


$cc=1;
foreach my $key (sort keys %hlst){
    undef(@vv);
    foreach my $name (@allnames){
        $var="$H{$key}{$name}$H{$key}{$name}";
        if($var eq "" ){$var="NN";}
        push(@vv, $var);
    }
    @ii = split(/&/, $key);
   $gg=join("\t",@vv);
    print OUT "$cc\t$sys\t$ii[0]\t$ii[1]\t+\tAGPv5\tGMI\t\t\tTEST_Panel\t\t$gg\n";
    $cc++;
}







exit;
print "IBS...\n";
foreach my $pos (sort { $H{$a} <=> $H{$b} } keys %H) {
    #print "<$pos>\n";
    foreach my $indA (keys %inx){
        foreach my $indB (keys %inx){
            if(!exists $H{$pos}{$indA} or !exists $H{$pos}{$indB} ){next;}
            if($H{$pos}{$indA} ne $H{$pos}{$indB}){
                $HMD{$indA}{$indB}=$HMD{$indA}{$indB}+1;
            }
            $HMT{$indA}{$indB}=$HMT{$indA}{$indB}+1;
        }#end foreacha matrx2
    }#end foreach matrx1
}#end pos
#print "<$pos> <$ind><$H{$pos}{$ind}>  <$ind2><$H{$pos}{$ind2}> <$HMD{$ind}{$ind2}>";
print "Done!!!\n";
print "Printing matrix...\n";
undef(@nax);

push(@nax, "IND");
foreach my $aa (sort keys %inx){
    push(@nax, $aa);
}
$sys=join("\t",@nax);

print OUT "$sys\n";

foreach my $ind (sort keys %inx){
    undef(@nax);
    push(@nax, $ind);

    foreach my $ind2 (sort keys %inx) {
        if($HMD{$ind}{$ind2} == 0){ push(@nax, 0);next;}
        $dist=$HMD{$ind}{$ind2}/$HMT{$ind}{$ind2};
        push(@nax, $dist);
    }
    $sys=join("\t",@nax);
    print OUT "$sys\n";
}