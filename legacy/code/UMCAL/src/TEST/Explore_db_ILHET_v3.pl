$file="/path/to/maizedb/ID_forHet_error.txt";

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
$file="/path/to/workdir/Production_ALL_Data_2023_10/MaizeDB_UME_production_20231010.vcf.gz.smp.vcf";

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


undef(%hre);
$cmis=0;
foreach my $aa (sort keys %hall){
    @colx=split(/:/, $col[$aa]);


    #################Alele freq filter

if($filt==1) {
    @dp = split(/,/, $colx[1]);
    $sumx = sum(@dp);
    if ($sumx < 9) {$colx[0] = "./."}
}
    ##############

#   print "ii <$col[0]><$col[1]>  <$aa> <$colx[0]> <$col[$aa]> <$hn{$aa}>   <$hl{$hn{$aa}}>  <@dp> <$sumx>\n";<STDIN>;

    if($colx[0] eq "./."){$cmis++;next}
    $hre{$colx[0]}++;
   #print "ii <$col[0]><$col[1]>  <$aa> <$colx[0]> <$ff>  <<$hre{$colx[0]}>>  <$col[$aa]> <$hn{$aa}>   <$hl{$hn{$aa}}>\n";<STDIN>;
}
    ###norm
    undef(%hf);
    foreach my $aa (keys %hre){
        $tt=$ct-$cmis;
        $norm=$hre{$aa}/$tt;
        $hf{$aa}=$norm;
    }

##freq

    $fuf=0;
    $fux=0;
    $nohet=0;

    foreach my $aa (keys %hf){

        $ff=0;
        if($aa eq "0/0" or $aa eq "1/1" or $aa eq "2/2" or $aa eq "3/3"){$ff=1;}
        elsif($aa ne "0/0" and $aa ne "1/1" and $aa ne "2/2" and $aa ne "3/3"){$ff=2;}
        else{print "ERROR 390847r <<$aa>>\n";<STDIN>}

        $hfx{$hf{$aa}}{$ff}=$hfx{$hf{$aa}}{$ff}+1;
        $hs{$hf{$aa}}=0;

        if($ff == 2 and $hre{$aa} eq 1){ $fuf=1;}
        if($ff == 1 and $hre{$aa} eq 1){ $fux=1;}

    # print "zx  <$col[0]><$col[1]>   <$aa>  rawval<$hre{$aa}> mis<$cmis> norm<<$hf{$aa}>>  cuant<$hfx{$hf{$aa}}{$ff}> flag<$ff>\n";<STDIN>;
if($ff==2){print OUTA "$hf{$aa}\n"; print OUTB "$hf{$aa}\n";$nohet=1;$gg++}

    }

    if($nohet == 0){print OUTB "0\n";$gg++}
$cct++;
    foreach my $aa (sort keys %hall){
        @colx=split(/:/, $col[$aa]);


        #################Alele freq filter
        if($filt==1) {
            @dp = split(/,/, $colx[1]);
            $sumx = sum(@dp);
            if ($sumx < 9) {$colx[0] = "./."}
        }
        ##############

        # print "ii <$col[0]><$col[1]>  <$aa> <$colx[0]> <$col[$aa]> <$hn{$aa}>   <$hl{$hn{$aa}}>\n";<STDIN>;
        if($colx[0] eq "./."){$sfmis{$hn{$aa}}++;next}

        if($colx[0] eq "0/0" or $colx[0] eq "1/1" or $colx[0] eq "2/2" or $colx[0] eq "3/3"){$ff=1;}
        elsif($colx[0] ne "0/0" and $colx[0] ne "1/1" and $colx[0] ne "2/2" and $colx[0] ne "3/3"){$ff=2;}
        else{print "ERROR 390847r <<$aa>>\n";<STDIN>}

      #  print "iix <$col[0]><$col[1]>  <$aa> <$colx[0]> <$ff>  <<$hre{$colx[0]}>>  <$col[$aa]> <$hn{$aa}>   <$hl{$hn{$aa}}>\n";#<STDIN>;

        if($ff==2 and  $hre{$colx[0]} ==1){$sfhet{$hn{$aa}}++;}#print "het bad $ff   $colx[0]  $hn{$aa}  <$sfhet{$hn{$aa}}>\n"

        if($ff==1 and  $hre{$colx[0]} ==1){$sfhom{$hn{$aa}}++; }#print "hom bad $ff   $colx[0]  $hn{$aa}  <$sfhom{$hn{$aa}}>\n"


        if($ff==2 ){$sahet{$hn{$aa}}++;}
        if($ff==1 ){$sahom{$hn{$aa}}++;}

    }

} # end file


print "<<<$gg>>>\n";

foreach my $aa (sort keys %sfhet){

$phet=$sfhet{$aa}/($cct-$sfmis{$aa});
    $phom=$sfhom{$aa}/($cct-$sfmis{$aa});

    $p2het=$sahet{$aa}/($cct-$sfmis{$aa});
    $p2hom=$sahom{$aa}/($cct-$sfmis{$aa});

    #print "fin <$aa> <$sfhet{$aa}> <$sfhom{$aa}>  <$sfmis{$aa}>  <$cct>\n";<STDIN>;

   $ppshet=$sfhet{$aa}/$sahet{$aa};
    $ppshom=$sfhom{$aa}/$sahom{$aa};

print OUT "$aa\t$phet\t$phom\t$p2het\t$p2hom\t$ppshet\t$ppshom\n";

}


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
