###########################Transpose VCF into a single line
#
###los getops estaran en el main.sh

$ref=$ARGV[0];
$file=$ARGV[1];
$bed=$ARGV[2];
$odir=$ARGV[3];
print "File:<$file>\n";
if ($file =~ /vcf$/){ open (FIL, "$file") or die ;}
elsif ($file =~ /vcf.gz$/){open (FIL, "gunzip -c $file |") or die;}
else{print "Input must be either .vcf or .vcf.gz\n";}
print "Checking files exist...\n";
print "BED file:<$bed>\n";
open(BED,"$bed") or die;
print "Reference:<$ref>\n";
open(REF,"$ref") or die;
$out="$odir/$file.tmpt";
#$ffile="$odir/$file.tmpt.gz";
print "OUT file:<$out>\n";
#if(-e $ffile){exit}
open(OUT,">$out") or die;
print "Done!...\n";
$outa="$file.ll";
print "OUT file:<$outa>\n";
open(LL,">$outa") or die;
print "Done!...\n";
#print "<$bed><$ref><$file><$out>\n";
print "loading bed...\n";
while($slin=<BED>) {
    chomp($slin);
    @lcol = split(/\t/, $slin);
    $hbed{$lcol[0]}="$lcol[2]";
    $ccll=$ccll+$lcol[2];
}
close(BED);
print "Done!...\n";
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
#######CONTROL for INFO to get use a different array pos each time you add one
$i2g[0]="DP=";
$h2g[0]="DP";
$i2g[1]="MLEAF=";
$h2g[1]="MLEAF";
#########################
print "Reading VCF...\n";
$f=0;
$achr=1;
$ccpp=0;
while($slin=<FIL>) {
    if ($slin =~ /^\n/) {next;}
    chomp($slin);
    if ($slin =~ /^\#/) {
        if ($f==0 and $slin ne "##fileformat=VCFv4.2"){print "ERROR format: looking for fileformat=VCFv4.2 file have $slin\n" and die;}
        $f=1;
        if ($slin =~ /^\#CHROM/) {
            @lcol = split(/\t/, $slin);
            if($lcol[0] ne "#CHROM"){print "ERROR format:Header wrong format<#CHROM><$slin>\n" and die;}
            if($lcol[1] ne "POS"){print "ERROR format:Header wrong format<POS><$slin>\n" and die;}
            if($lcol[2] ne "ID"){print "ERROR format:Header wrong format<ID><$slin>\n" and die;}
            if($lcol[3] ne "REF"){print "ERROR format:Header wrong format<REF><$slin>\n" and die;}
            if($lcol[4] ne "ALT"){print "ERROR format:Header wrong format<ALT><$slin>\n" and die;}
            if($lcol[5] ne "QUAL"){print "ERROR format:Header wrong format<QUAL><$slin>\n" and die;}
            if($lcol[6] ne "FILTER"){print "ERROR format:Header wrong format<FILTER><$slin>\n" and die;}
            if($lcol[7] ne "INFO"){print "ERROR format:Header wrong format<INFO><$slin>\n" and die;}
            if($lcol[8] ne "FORMAT"){print "ERROR format:Header wrong format <FORMAT><$slin>\n" and die;}
        }#end if #chrom
        next;
    } #end of satring with #
    @lcol = split(/\t/, $slin);
    $chr=$lcol[0];
    $chr=~ s/chr//;
    $pos=$lcol[1];
    $ref = $lcol[3];
    $alt = $lcol[4];
    $alt =~ s/\./$ref/;
    $alt =~ s/<NON_REF>//;
    $qual=sprintf("%.4f", $lcol[5]);
    $info = $lcol[7];
    $gref=substr $href{$chr},$pos-1, 1;

    #### Print chromosome to file and empty hash
    if($achr ne $chr) {
        print "$achr...\n";
        #print "$slin...\n";<STDIN>;
        for ($i = 1; $hbed{$achr} >= $i; $i++) {
         if(exists $hh{$i}){print OUT "$hh{$i}\n"; $ccpp++;$chpp{$achr}++;}
        else{print OUT "::::\n"; $ccpp++;$chpp{$achr}++;}
        }
        undef(%hh);
        $achr=$chr;
    }

    #INDELS  ### this needs to be recoded
    if(length($ref) > 1 ){next}

    #check ref
    if($ref eq "."){$ref=$gref}
    if($info !~ /END=/ and $ref !~ /^$gref/){print "ERROR reference: expecting <$gref>  you have <$ref> line: <$slin>" and die; } ### this is done due deepvar ref deletions

    ### deal with different subformats in FORMAT column
    @format=split(/:/,$lcol[8]);
    @samp=split(/:/,$lcol[9]);
    if(@format ne @aformat){
        $cc=0;
        #print "1 @format   @aformat\n";<STDIN>;
        foreach my $aa (@format){
            $cformat{$aa}=$cc;
            #  print "2) $aa $cformat{$aa}\n";<STDIN>;
            $cc++;
        }
    }
    @aformat=@format;
    ###### Extract info from format
    foreach my $aa (keys %cformat){
        $hformat{$aa}=$samp[$cformat{$aa}]; # print "$aa $hformat{$aa}"; <STDIN>;
    }
    ###########
    #### deal with different subformats in INFO column
    #https://stackoverflow.com/questions/6237968/how-to-extract-the-text-between-two-patterns-using-regex-perl
    $cc=0;
    foreach my $vv (@i2g){
        $ss=$vv;
        $ee=";";
        $info =~ /$ss(.*?)$ee/;
        $hi{$h2g[$cc]}=$1;
        $cc++;
    }
    #### get alts
    undef(@alts);
    $alts[0]=$ref;
    @posb=split(/,/,$alt);
    push(@alts, @posb);

    #####FILTERS
    # if ($qual eq "."){next;}
    # if ($qual < 6){next;}
    if ($hi{DP} < 1 and  $hformat{DP} < 1 ){next;}

    #################RECODE INDELS
    if(length($alts[0])>1){$alts[0]="-"}
    $cc=0;
    foreach my $xx (@alts){
        if($cc==0){$cc++;next;}
        #print "111 $cc $xx  $alts[$cc]\n";
        if($xx =~ /<NON_REF>/){$cc++; next;}
        if($xx =~ /INS/){$alts[$cc]="+"}
        if($xx =~ /DEL/){$alts[$cc]="-"}
        if(length($xx)>1){$alts[$cc]="+"}
        #print "222 $cc $xx $alts[$cc]\n";
        $cc++;
    }
    $all=join ",", @alts;
    #### recode genotype
    chomp($hformat{GT});
    $hformat{GT} =~ s/\./$ref/g;
    @GT=split(/\//,$hformat{GT});
    $est = $hformat{GT};
    $geno="$alts[$GT[0]]$alts[$GT[1]]";

   # print "aaaax est<$est>  gt<@GT> geno<$geno>  ref<$ref> alt<$alt> gref<$gref>  alts<@alts>  line<$slin>\n";<STDIN>;

    #checks genotype
    if($geno =~ /<NON_REF>/){print "Warning <NON_REF> allele in $slin\nTEST<$lcol[3]><@posb><@alts><$geno><$GT[0]><$GT[1]>\n";}
    if(length($geno) > 2 ){ print "ERROR format: geno <$geno> more than 2 aleles. this asumes diploid  <$slin>" and die;}
    if($est =~ /\.\/\./) {next;print "ERROR 235: geno <$geno> should not be ./. <$slin>" and die;}


    if($info !~ /END=/){
        $hh{$pos}="$geno:$qual:$hformat{DP}:$hformat{AD}:$hformat{PL}";
      #  print "XXX <$pos> <$hh{$pos}> $slin...\n";<STDIN>;
    } #end if SNV not END


    if($info =~ /END=/){
        @aaa=split(/END=/,$info);
        @bbb=split(/;/,$aaa[1]);
        $END=$bbb[0];
        $hfor="$hformat{DP},$hformat{DP}";
        for($i=$pos;$i<=$END;$i++){
            $gref=substr $href{$chr},$i-1, 1;
            $geno=$gref.$gref;
            $hh{$i}="$geno:$qual:$hformat{DP}:$hfor:$hformat{PL}";
            #print "yyy <$i> <$hh{$i}> $slin...\n";<STDIN>;
        }
    } #end if SNR END


} #end While FIL  vcf

$chr="LAST_CHR";

if($achr ne $chr) {
    print "$achr...\n";

    for ($i = 1; $hbed{$achr} >= $i; $i++) { ### this is what makes print all the chromosomes sinche the VCF is sorted the last line has the last chromosome
        # $gref = substr $href{$achr}, $i - 1, 1;  #print " <$achr> <$hbed{$achr}> <$ctrl{$achr}> <$i> <$gref>\n";<STDIN>;
        if(exists $hh{$i}){print OUT "$hh{$i}\n"; $ccpp++;$chpp{$achr}++;}
        else{print OUT ":::\n"; $ccpp++;$chpp{$achr}++; }
        $ctrl{$achr}=$i;
    }

    undef(%hh);
    $achr=$chr;
}

####################### CHECK metrics of alignment #########################

if($ccll ne $ccpp){print "CRITICAL ERROR: Aligment mismatch!!!!!  there is something in the VCF format not concidered in the script, therfore it creates diferent positions than needed in the genome aligment... FIX it  before proceeding to downstream analysis";  }

foreach my $zz (keys %chpp){
    if($chpp{$zz} ne $hbed{$zz}){ print " CRITICAL ERROR:  chr <$zz>  you have: <$chpp{$zz}>  should be :<$hbed{$zz}> Aligment mismatch!!!!!  there is something in the VCF format not concidered in the script, therfore it creates diferent positions than needed in the genome aligment... FIX it  before proceeding to downstream analysis";  }
}

print LL "$ccll\t$ccpp\n";

print "Done!...\n";

close(OUT);
close(LL);

sub log10 {
    my $n = shift;
    return log($n)/log(10);
}