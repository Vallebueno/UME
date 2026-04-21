$file="/path/to/maizedb/ID.lst";



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


#open(FIL,"gunzip -c $file |") or die;
open(FIL,"$file") or die;
open(OUT,">$out") or die;


$hff{1}=0;
$hff{2}=0;

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
   # print "ii <$col[0]><$col[1]>  <$aa> <$colx[0]> <$col[$aa]> <$hn{$aa}>   <$hl{$hn{$aa}}>\n";<STDIN>;
    if($colx[0] eq "./."){$cmis++;next}
    $hre{$colx[0]}++;
  #  print "ii <$col[0]><$col[1]>  <$aa> <$colx[0]> <$ff>  <<$hre{$colx[0]}>>  <$col[$aa]> <$hn{$aa}>   <$hl{$hn{$aa}}>\n";<STDIN>;
}
    ###norm
    undef(%hf);
    foreach my $aa (keys %hre){
        $tt=$ct-$cmis;
        $norm=$hre{$aa}/$tt;
        $hf{$aa}=$norm;
    }

##freq
    foreach my $aa (keys %hf){

        $ff=0;
        if($aa eq "0/0" or $aa eq "1/1" or $aa eq "2/2" or $aa eq "3/3"){$ff=1;}
        elsif($aa ne "0/0" and $aa ne "1/1" and $aa ne "2/2" and $aa ne "3/3"){$ff=2;}
        else{print "ERROR 390847r <<$aa>>\n";<STDIN>}

        $hfx{$hf{$aa}}{$ff}=$hfx{$hf{$aa}}{$ff}+1;
        $hs{$hf{$aa}}=0;

       print "zx  <$col[0]><$col[1]>   <$aa>  rawval<$hre{$aa}> mis<$cmis> norm<<$hf{$aa}>>  cuant<$hfx{$ff}{$hf{$aa}}> flag<$ff>\n";<STDIN>;

    }



}

foreach my $aa (sort keys %hs){
    foreach my $ff (sort keys %hff){

        #if(!exists $hfx{$aa}){}

      print "<$ff>  <$aa>  <$hfx{$hf{$aa}}{$ff}>\n";<STDIN>;


    }

}




close(OUT);
close(FIL);

print "Done!!!\n";
