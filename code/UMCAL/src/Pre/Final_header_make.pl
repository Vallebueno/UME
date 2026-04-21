$list=$ARGV[0];
$ref=defined $ARGV[1] ? $ARGV[1] : "";
$out="$list.hh";

print "making header...\n";

open(LST,"$list") or die;
while($line=<LST>) {

    chomp($line);
    @lcol = split(/\t/, $line);

    $lcol[0] =~ s/_mappedAGPv5.sorted.MarkedDup.bam.ALL.mpileup.gz//;
    $lcol[0] =~ s/_mappedAGPv5.sorted.MarkedDup.bam.mpileup.gz//;
    $lcol[0] =~ s/.tmpt.gz//;
    push(@ff,$lcol[0]);
}

$sys=join("\t", @ff);

#print "$sys"; <STDIN>;

open(OUT,">$out") or die;

print OUT "##fileformat=VCFv4.2\n";
print OUT '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n';
print OUT "##reference=$ref\n";
print OUT "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t$sys\n";

close(OUT);
close(LST);
