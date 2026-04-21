$file="/path/to/workdir/ALL/Test1_.ALL.Merge.10.poli.Union.db";

$out="$file.vcf";

open(OUT,">$out") or die;
#open(DB,"gunzip -c $file |") or die;
open(DB,"$file") or die;

$datestring = localtime();

print '##fileformat=VCFv4.2
##fileDate=$datestring
##reference=file:////path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa
##Mapped against B73-NAM AGPV5
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AD,Number=A,Type=Float,Description="Allele Depth">
';

print "Trasforming into VCF...\n";

while($line=<DB>) {
    chomp ($line);

    if($line =~ /^#/){
        if ($line =~ /^#0=missing info/){next}
        if ($line =~ /^#1=reference allele/){next}
        if ($line =~ /^#Mapped against B73-NAM AGPV5/){next}
        @heads = split(/\t/, $line);
        $len=scalar(@heads);
        next;
    }

    undef (@alts);
    push(@alts, $ref);

    if($line=~ /A/ and $ref ne "A"){push(@alts, "A");}
    if($line=~ /C/ and $ref ne "C"){push(@alts, "C");}
    if($line=~ /T/ and $ref ne "T"){push(@alts, "T");}
    if($line=~ /G/ and $ref ne "G" ){push(@alts, "G");}
    if($line=~ /\+/){push(@alts, "X"); $line=~ s/\+/X/g;}
    if($line=~ /^\-/){push(@alts, "Y"); $line=~ s/^\-/Y/g;}
    if($line=~ /,\-/ ){push(@alts, "Y"); $line=~ s/,\-/Y/g;}

    @lcol = split(/\t/, $line);
    $chr=$lcol[0];
    $pos=$lcol[1];
    $ref=$lcol[2];

    if($chr ne $achr){print "$chr\n"; $achr=$chr}

    undef(%halt);

    $c=0;
    foreach my $xx (@alts){
        $halt{$xx}=$c;
        $c++;
    }
    print "<$chr><$pos>\n";

    undef(@ll);
    for ($i=3;$i<$len;$i++){

        if ($lcol[$i] eq "0"){$ll[$i]="./."; next}
        @inf = split(/:/, $lcol[$i]);
        @GT = split(//, $inf[0]);
        $GT[0] =~ s/$GT[0]/$halt{$GT[0]}/;
        $GT[1] =~ s/$GT[1]/$halt{$GT[1]}/;

        if($GT[0] > $GT[1]){$ll[$i]="$GT[1]/$GT[0]";}
        elsif($GT[1] > $GT[0] or $GT[1] == $GT[0]){$ll[$i]="$GT[0]/$GT[1]";}

        print "<$i> orig<$lcol[$i]> <$heads[$i]>   ref<$ref> alts<@alts>    <@inf[0]>  transf<$ll[$i]>\n";

        if ($ll[$i] eq "1/0"){print "ERROR 284u4\n";<STDIN>;}

    }

    $sys=join("", @ll);
    $asys=join(",", @alts);

    print "$chr\t$pos\t.\t$ref\t$asys\t$qual\t.\t";

}

print "END\n";

close (DB);
