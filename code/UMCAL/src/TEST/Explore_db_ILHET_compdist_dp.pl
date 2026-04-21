$file="/groups/swarts/lab/MAVE/MAIZEDB/ID_forHet_error.txt";

print "Loading list...\n";
open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    chop($line);
    @col=split(/\t/, $line);

    if($col[1] ne "I1"){next;}
    $hl{$col[0]}=$col[1];


    # print "$col[0] $col[1]";<STDIN>;
}
close(FIL);
print "Done!!!\n";
print "Reading file...\n";
$file="/scratch-cbe/users/miguel.vallebueno/Production_ALL_Data_2023_10/MaizeDB_UME_production_20231010.vcf.gz.smp.vcf";

$out="$file.IHET.cuantex";

$outa="$file.IHET.HSFS.cuantex";
$outb="$file.ISFS.cuantex";

$filt=0;

#open(FIL,"gunzip -c $file |") or die;
open(FIL,"$file") or die;
open(OUT,">$out") or die;
open(OUTA,">$outa") or die;
open(OUTB,">$outb") or die;
print OUT "sample\tshet\tshom\tahet\tahom\tpshet\tpshom\n";
print OUT "$aa\t$phet\t$phom\t$p2het\t$p2hom\t$ppshet\t$ppshom\n";
$hff{1}=0;
$hff{2}=0;

$chetx=0;
$chomx=0;
$cct=0;

$gg=0;
$cbss=0;
$cthg=0;
$cuhg=0;

$sigm1=0.6826;
$sigm2=1.49;
while($line=<FIL>) {

    chomp($line);

    if($line =~ /^#C/){
        @col=split(/\t/, $line);
        $tam=scalar(@col);
        $ct=0;
        for($i=9;$i<$tam;$i++){
            if(!exists $hl{$col[$i]}){next}
            #print "ee $i $col[$i] $hl{$col[$i]}\n";<STDIN>;
            $hall{$i}=0;
            $hn{$i}=$col[$i];
            $ct++;
        }
    }
    if($line =~ /^#/){next;}

    @col=split(/\t/, $line);


    ##################################  unfiltered
    undef(%hre);
    $cmis=0;
    undef(%hf);
    foreach my $aa (sort keys %hall){
        @colx=split(/:/, $col[$aa]);
        if($colx[0] eq "./."){$cmis++;next}
        $hre{$colx[0]}++;
        $hf{$colx[0]}=1;
    } # end foreach %hall


    ################################## filtered

    undef(%hreff);
    $cmisff0;
    undef(%hfff);
    foreach my $aa (sort keys %hall){
        $gg++;
        @colx=split(/:/, $col[$aa]);
        @dp = split(/,/, $colx[1]);
        $sumx = sum(@dp);



        if($colx[0] ne "0/0" and $colx[0] ne "1/1" and $colx[0] ne "2/2" and $colx[0]ne "3/3"){
            $cthg++;
        ### balance btween het read counts
            @hgen=split(/\//, $colx[0]);

           $dp1=$dp[$hgen[0]];
            $dp2=$dp[$hgen[1]];


            if($dp1 == 0 or $dp2 == 0){next;}
            $rel=$dp1/$dp2;


            #print "$col[$aa]   $dp1 $dp2  @hgen   <$rel>";<STDIN>;
            if($rel < $sigm1 or $rel >$sigm2){$cuhg++;}

            # $sigm1=0.6826;
            # $sigm2=1.49;
    }




        if ($sumx < 9) {$colx[0] = "./."}
        if($colx[0] eq "./."){$cmisff++;next}#print "set to miss <$col[0]><$col[1]>   <$aa>  <$col[$aa]>  <$colx[0]>\n";
        $hreff{$colx[0]}++;
        $hfff{$colx[0]}=1;
    } # end foreach %hall

    #print "comparing dists...\n";
    $allg="X";
    undef(@pbadg);
    $fbad=0;


    foreach my $aa (keys %hf){


        $ff=0;
        if($aa eq "0/0" or $aa eq "1/1" or $aa eq "2/2" or $aa eq "3/3"){$ff=1;$allg=$allg.$aa;} ### can improve this now is only homos maybe other hets can be taken? dificult to code
        elsif($aa ne "0/0" and $aa ne "1/1" and $aa ne "2/2" and $aa ne "3/3"){$ff=2;}
        else{print "ERROR 390847r <<$aa>>\n";<STDIN>}

        $hfx{$hf{$aa}}{$ff}=$hfx{$hf{$aa}}{$ff}+1;
        $hs{$hf{$aa}}=0;



      #  print "zsx <$allg> <$col[0]><$col[1]>   <$aa> exists in hid<$hfff{$aa}> rawval in hid<$hreff{$aa}> rawval<$hre{$aa}> mis<$cmis> norm<<$hf{$aa}>>  cuant<$hfx{$hf{$aa}}{$ff}> flag<$ff>\n"; <STDIN>;

        if($ff == 2 and $hre{$aa} eq 1){
            #   print "ttt  <$col[0]><$col[1]>   <$aa> exists in hid<$hfff{$aa}> rawval in hid<$hreff{$aa}>  rawval<$hre{$aa}> mis<$cmis> norm<<$hf{$aa}>>  cuant<$hfx{$hf{$aa}}{$ff}> flag<$ff>\n";
            push(@pbadg,$aa);
            $fbad=1;
        }
    }

    if($fbad==1){
        $good=0;
        foreach my $cc (@pbadg){
            #  print "233453454 <$cc>\n";
            @colg=split(/\//, $cc);

            # print "aa  <$col[0]><$col[1]>   <$aa> exists in hid<$hfff{$aa}> rawval in hid<$hreff{$aa}>   rawval<$hre{$aa}> mis<$cmis> norm<<$hf{$aa}>>  cuant<$hfx{$hf{$aa}}{$ff}> flag<$ff>\n";
            # print "<@colg>   <$allg>\n";
            if($allg =~ /$colg[1]/ and $allg =~ /$colg[0]/){ $good=1; }
            if($good==1){next;}
            $cbss++;
            $fc=$cbss/$gg;

            print "bad <$cc>  rawval in hid<$hreff{$cca}>  rawval<$hre{$cc}>           <$cbss> total<$gg>  fc<$fc>\n";

        }

    }
}

$pert=$cuhg/$cthg;
print "proportion disbalanced Hets    <$cuhg> <$cthg>    <$pert>\n";

print "";
close(OUT);
close(OUTA);
close(FIL);

print "Done!!!\n";

sub sum {
    my @vals = @_;
    my $result;
    foreach my $i (@vals){
        $result = $result + ($i);
    }
    return $result;
}