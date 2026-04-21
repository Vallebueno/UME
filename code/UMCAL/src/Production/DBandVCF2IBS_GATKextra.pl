#for i in {max,min,median,mean,sum}; do echo $i; time perl $caronte/Modules/PopGene/DBandVCF2IBS.pl $i; done;



#read DB with raw snps to extract GATK sites
$ccox=20;

$filex="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone.Merge.db.7.poli.db.gz";
$list="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone";
my @colors = (1,2,3,4,5,6,7,8,9,10);




$file = "/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone.Merge.db.7.poli.qual0.Union.V7-1.max.db";

open(VCF, "$file") or die;


print "Sampling sites...\n";
$c=0;
while($line=<VCF>) {
    if ($line =~ /^#/) {next;}
    chomp($line);
    @lcol = split(/\t/, $line);
    $chr = $lcol[0];
    $chr =~ s/chr//;
    $chr =~ s/ //;
    $pos = $lcol[1];
    $key = "$chr&$pos";

    my $color = int(rand(100));
    if ($color > $ccox) {next;}
$c++;
    $poss{$key}=1;
}

print "sampled <$c>sites\n";
close(VCF);
print "Done\n";
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
    if($aa=~ /30xgatk/){$inx{$aa}=1;$hgatk{$c}=$aa;}
    if($aa=~ /0.5xgatk/){$inx{$aa}=1;$hgatk{$c}=$aa;}
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
 #   print "<$chr><$pos> <$line>\n";<STDIN>;
    $ref = $lcol[2];
    $len=scalar(@lcol);
    undef(%HoA);
    undef(%HoB);
    for ($i = 3; $i < $len; $i++) {
        $h = $i - 3;
        $samp = $hgatk{$h};
        $val = $lcol[$i];
        $key="$chr&$pos";
        if(!exists $poss{$key}){next;}
undef(@ran);
        @XA = split(/:/, $lcol[$i]);
        @GT = split(//, $XA[0]);
        push(@ran, $GT[0]);
        push(@ran, $GT[1]);
        $gen=$ran[rand @ran];
        if($gen eq ""){next;}

        if(exists $hgatk{$h}){$H{$key}{$samp}=$gen;}

    }
}
close (DB);

print "Done!\n";




#### read Union


print "start!!!\n";
$algx="max";


$out="MATRIX.$algx.testsubplusGATK.V7.IBS.txt4";
open(OUT,">$out") or die;



print "Done!!!\n";

undef(%ht);
undef(@at);


#############################################################     DB
$ilst="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/inds_2plot.lst";

#$ilst="/groups/swarts/lab/MAVE/DB/Union/inds_2plot.lst";

print "loading Individuals list...\n";
open(IL,"$ilst") or die;
$cct=0;
$cin=0;
while($line=<IL>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    if ($line =~ /##/) {next;}

    $ht{$lcol[0]}=1;
    $cin++;
}
print "inds to be considered <$cin>\n";
close (IL);
print "Done!!!\n";

#$algs[0]="sum";
#$algs[1]="min";
#$algs[2]="max";
#$algs[3]="mean";
#$algs[4]="median";

#$cof[0]=0;
#$cof[1]=1;
#$cof[2]=2;
#$cof[3]=3;
#$cof[4]=4;
#$cof[5]=5;
#$cof[6]=6;
#$cof[7]=7;
#$cof[8]=8;
#$cof[9]=9;
#$cof[10]=10;

$cof[0]=0;
$cof[1]=3;
$cof[2]=5;
$cof[3]=10;

#foreach my $algx (@algs) {


foreach my $cutoff (@cof){

    print "Reading VCF cutoff $cutoff...\n";

     $file = "/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone.Merge.db.7.poli.qual0.Union.V7-1.max.db";

    #$file = "/groups/swarts/lab/MAVE/DB/Union/1.ALL.Merge.2.poli.qual0.Union.V7-1.max.db.gz";
    #open(VCF, "gunzip -c $file |") or die;
    open(VCF, "$file") or die;
    $r=0;
    while($line=<VCF>) {
        if ($line =~ /^##/) {next;}
        chomp($line);
        @lcol = split(/\t/, $line);

        if ($line =~ /#chr/) {
            $f = 0;
            if ($r == 1) {next;}
            $r = 1;
            $ccfs = 0;
            foreach my $aa (@lcol) {
                if (exists $ht{$aa}) {
                    #print "<$aa>  <$ht{$aa}> <$f>\n";<STDIN>;
                    $ht{$aa}=$f;
                    push(@at,$f);
                    $ccfs++;
                    $hDB{$f}=$aa;
                    $namx="$hDB{$f}&$cutoff";
                    $inx{$namx}=1;
                } #save pos of array for neded inds
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
        $ref=$lcol[2];

        if(!exists $poss{$key}){next;}



        $cflt = 0;
        foreach my $aa (@at) {
            undef(@ran);
            $namx="$hDB{$aa}&$cutoff";
            @ii = split(/_/, $namx);
            @rin = split(/&/, $ii[4]);

            if ($lcol[$aa] eq "0"){ next;}
            # if ($lcol[$aa] eq "0" and exists $H{$pos}{}) {$lcol[$aa] = "NN"}
            #if($lcol[$aa] =~ "$ref$ref"){next;}

            @info = split(/:/, $lcol[$aa]);
            $qualx = $info[1];

            if ($qualx <= $cutoff) { next;}
            @GT = split(//, $info[0]);
            $geno="$GT[0]$GT[1]";
            push(@ran, $GT[0]);
            push(@ran, $GT[1]);
            $gen=$ran[rand @ran];
            #print "XXXXXXXXX pos<$key> ran<@ran>  pind<$aa>  ind<$namx>  infoin<$lcol[$aa]>  GT1<$GT[0]> GT2<$GT[1]> geno<$geno> <$gen>   <$algx><$cutoff>\n";
            $H{$key}{$namx}=$gen;
        }
    } #END VCF
    close(VCF);
}#END COFS
#} #END ALGS

print "Done!!!\n";

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