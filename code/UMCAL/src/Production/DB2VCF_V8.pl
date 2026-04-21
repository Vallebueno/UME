$head=$ARGV[0];  #f
$head="$head.tlone";

$list=$ARGV[1];  #k
$ref=$ARGV[2];   #r


$file="All.tlone.Merge.db.gz";
#$file="TEST.gz";

#$head="/scratch-cbe/users/miguel.vallebueno/mpilp2db/mpileups.lst.tlone";


print "Extracting headers from input list...\n";
open(HED,"$head") or die;
chomp(my @LLL = <HED>);
@lcol = split(/\|/, $LLL[0]);
$head=$lcol[0];
$head =~ s/paste //g;
$head =~ s/ \S+\// /;
$head =~ s/ //;
$head =~ s/\<\(zcat//g;
$head =~ s/\)//g;
$head =~ s/.7.mpileup.ol//g;
$head =~ s/.sorted.MarkedDup.bam.ALL.mpileup.gz.tmpt.gz//g;
$head =~ s/.sorted.MarkedDup.bam.mpileup.gz.tmpt.gz//g;

$head =~ s/ \S+\// /g;
$head =~ s/> \S+//g;
@heads = split(/\s+/, $head);

#print "@heads"; <STDIN>;

close (HED);

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


open(LST,"$list") or die;
print "Loading sites...\n";
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]&$lcol[2]";
    push(@hh,$key);
}
close(LST);
print "DONE!!!!!\n";
print "DB2VCF...\n";
open(DB,"gunzip -c $file |") or die;
$out="$file.vcf";
open(OUT,">$out") or die;
$cx=0;
$sys=join("\t", @heads);
print OUT "##fileformat=VCFv4.2\n";
print OUT '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths for the ref and alt alleles in the order listed">\n';
print OUT '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n';
print OUT '##INFO=<ID=END,Number=1,Type=Integer,Description="Stop position of the interval">\n';
print OUT "##reference=$ref\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sys\n";
while($line=<DB>) {
    chomp ($line);
    @lcol = split(/\t/, $line);
    $da=$hh[$cx];
    @data = split(/&/, $da);
    $chr=$data[0];
        $pos=$data[1];
    @vars=split(/ /,$data[2]);
    $gref=substr $href{$chr},$pos-1, 1;
    undef(@alts);
    undef(@altx);
    push(@alts,$gref);
    foreach my $ee (@vars){
        if($ee eq $gref){next;}
        push(@alts,$ee);
        push(@altx,$ee);
    }
    $cvc=0;
    undef(%halt);
    $halt{"N"}=".";
    $halt{" "}=".";
    foreach my $ee (@alts){
        $halt{$ee}=$cvc;
        $cvc++;
    }
    $cx++;
    if($chr ne $achr){print "$chr\n"; $achr=$chr}
    $len=scalar(@lcol);
    undef(@ll);
    for ($i=0;$i<$len;$i++){
        @inf = split(/:/, $lcol[$i]);
        @GT = split(//, $inf[0]);
        $GT[0] =~ s/$GT[0]/$halt{$GT[0]}/;
        $GT[1] =~ s/$GT[1]/$halt{$GT[1]}/;
        if($GT[0] !~ /\d/ and  $GT[1] !~ /\d/){$GT[0]=".";$GT[1]=".";}
        elsif($GT[0] !~ /\d/ ){$GT[0]=$GT[1]}
        elsif($GT[1] !~ /\d/ ){$GT[1]=$GT[0]}
        $ll="$GT[0]/$GT[1]";
        if($GT[0] > $GT[1]){$ll="$GT[1]/$GT[0]"; }
        $geno[$i]="$ll";
        if($geno[$i] eq "/"){$geno[$i]="./."}
    }
    $sys=join("\t", @geno);
    $asys=join(",", @altx);
    print OUT "$chr\t$pos\t.\t$gref\t$asys\t.\tPASS\t.\tGT\t$sys\n";
}
print "END\n";
close (DB);