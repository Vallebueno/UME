# UME makes the production call generating final VCF file
# Author: Miguel Vallebueno
# Date: 2022-07-29

$filen=$ARGV[0];
$filen =~ s/.db.gz//g;
$chr=$ARGV[1];
$cof=$ARGV[2];
$alg=$ARGV[3];
$qualf=$ARGV[4];
$file="$filen.$chr.poli.qual$qualf.Union.V8.$alg.dbx";
$list="$filen.$chr.poli.qual$qualf.Union.V8.$alg.dbx.cof$cof.lst2";
$out="$filen.$chr.poli.cof$cof.Union.V8.$alg.production.vcf";

open(OUT,">$out") or die;
print "Loading DB of known variants...\n";
open(LST,"$list") or die;
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]";
    $h{$key}=$lcol[2];
    #print "<$key> <$h{$key}>\n"; <STDIN>;
}
close(LST);
print "Done!!!\n";
print "Production calling...\n";

open(FIL,"$file") or die;
$cc=0;
while($line=<FIL>) {
    chomp($line);
    @lcol = split(/\t/, $line);


    if($line =~ /^#chr/){
        $tam=scalar(@lcol);
        for($i=3;$i<$tam+1;$i++){push(@heads,$lcol[$i])}
        $sys=join("\t", @heads);
        print OUT "##fileformat=VCFv4.2\n";
        print OUT '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n';
        print OUT '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n';
        print OUT '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n';
        print OUT "##reference=$ref\n";
        print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sys\n";
    }


    #if($cc==0){$tam=scalar(@lcol);$cc++;if($tam != $lenh){print "ERROR:HEADER. mismatch between header <$lenh> and number of columns in DB file <$tam>";die}}


    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]";
    if(!exists $h{$key}){next;}
    $ref=$lcol[2];
    undef(@alts);
    undef(@altx);
    push(@alts,$ref);
    @vars = split(/ /, $h{$key});
    foreach my $ee (@vars){  #extract vector of alternative aleles
        if($ee eq $ref){next;}
        push(@alts,$ee);
        push(@altx,$ee);
    }
    $cvc=0;
    undef(%halt);
    foreach my $ee (@alts){ #transform vector of alleles into ordered nums
        $halt{$ee}=$cvc;
        $cvc++;
    }
    undef(@array);
    for($i=3;$i<$tam;$i++){
        #print "<$i>  <$tam>  <$key> <$lcol[$i]> <$h{$key}>";<STDIN>;
        if($lcol[$i] eq ":::"){push(@array, ".");next;}
        if($lcol[$i] eq ""){push(@array, ".");next;}
        if($lcol[$i] eq "0"){push(@array, ".");next;}
        @col = split(/:/, $lcol[$i]);
        @xcol = split(//, $col[0]);
        @ycol = split(/,/, $col[3]);
        $xcc=0;
        foreach my $xcx (@ycol){
            if($xcx ne $xcx+0){$ycol[$xcc]=1;$xcc++;}
        }

        $fa=0;$fb=0;
        if($lcol[$i] =~/\+/){push(@array, ".");next;}
        if($h{$key} =~ /$xcol[0]/){$fa=1}
        if($h{$key} =~ /$xcol[1]/){$fb=1}
         $gt=".";
           if($fa == 1 and $fb == 1){$gt="$halt{$xcol[0]}/$halt{$xcol[1]}:$ycol[0],$ycol[1]"}
        elsif($fa == 1 and $fb == 0){$gt="$halt{$xcol[0]}/$halt{$xcol[0]}:$ycol[0],$ycol[0]";}
        elsif($fa == 0 and $fb == 1){$gt="$halt{$xcol[1]}/$halt{$xcol[1]}}:$ycol[1],$ycol[1]";}
        elsif($fa == 0 and $fb == 0){push(@array, ".");next;}
        else{print "ERROR 49o857\n";print "wfrvv <$h{$key}> <$fa> <$fb>  <$lcol[$i]>  <$gt>";<STDIN>;die}
        push(@array, $gt);
    }
    $sys=join("\t", @array);
    $asys=join(",", @altx);
    print OUT "$lcol[0]\t$lcol[1]\t.\t$ref\t$asys\t.\tPASS\t.\tGT:AD\t$sys\n";
}
print "Done!!!\n";
close (DB);
close (OUT);
print "END\n";

