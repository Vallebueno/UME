# UME makes the production call from mpileups generating one column files guided by the Discovery list file. need to be joined afterwise into vcf
# Author: Miguel Vallebueno
# Date: 2023-01-10

#$ref="/path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL";
#$file="/path/to/workdir/mpilp2db/SM8_NoUDG_Vallebueno_SM8_mappedAGPv5.sorted.MarkedDup.bam.ALL.mpileup.gz";
#$list="/path/to/workdir/mpilp2db/1.ALL.Merge.10.poli.qual0.Union.V8.max.dbx.cof5.lst";

$IDIR=$ARGV[0];
$filen=$ARGV[1];
$ref=$ARGV[2];
$list=$ARGV[3];
$Odir=$ARGV[4];

$file="$IDIR/$filen";

open(TST,"$ref") or die;
close(TST);
open(TST,"$list") or die;
close(TST);

$outx="$Odir/$filen.tmpt.gz";

if (-e $outx) {
    print "the file exists\n";
 #   exit;
}

#open(TST,"$file") or die;
#close(TST);

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
    $key="$lcol[0]&$lcol[1]";
    $h{$key}=$lcol[2];
  #  $gref=substr $href{$lcol[0]},$lcol[1]-1, 1;
}
close(LST);
print "Done!!!\n";




print "reading mpileup <$file>...\n";
open(FIL,"gunzip -c $file |") or die;
print "Checkpoint 1 <$file>...\n";

while($line=<FIL>){
    $dp=0;
    chomp($line);
    if($line =~ /^#/){next;}
    @col = split(/\t/, $line);
    $chr=$col[0];
    #if($chr !~ /^10/){next;}
    if($line =~/INDEL/){next;}
    $pos=$col[1];
    $key="$col[0]&$col[1]";
    $info=$col[7];

    if($info =~ /END=/){
        @aaa=split(/END=/,$info);
        @bbb=split(/;/,$aaa[1]);
        $END=$bbb[0];
        @xaa=split(/MinDP=/,$info);
        @xbb=split(/;/,$xaa[1]);
        $mdp=$xbb[0];
        for($i=$pos;$i<=$END;$i++){
            $key="$chr&$i";
            if(exists $h{$key}){
                $gref=substr $href{$chr},$i-1, 1;
                if($h{$key} !~ /$gref/){next;}
                $geno=$gref.$gref;
                $hall{$key}="$geno";
            }
        }
        next;
    } #end if SNR END
    if(exists $h{$key}) {
        ####weird exceptions without AD
        if ($info !~ /AD/ and $col[4] eq "<*>") {
            if (exists $h{$key}) {
                $geno = $col[3].$col[3];
                $hall{$key}="$geno";
            }
            next;
        }
        if ($info =~ /AD=0/ and $col[4] eq "<*>") {
            $geno = $col[3].$col[3];
            $hall{$key}="$geno";
            next;
        }
        if($info !~ /AD/){
            $geno = $col[3].$col[3];
            $hall{$key}="$geno";
            next;
        }
        #end exists weird exceptions
        ###############polimorphic
        undef(@alts);
        push(@alts,$col[3]);
        @al=split(/,/,$col[4]);
        push(@alts,@al);
        undef(%hformat);
        @format=split(/:/,$col[8]);
        @samp=split(/:/,$col[9]);
        if(@format ne @aformat){
            $cc=0;
            foreach my $aa (@format){
                $cformat{$aa}=$cc;
                $cc++;
            }
        }
        @aformat=@format;
        ###### Extract info from format
        foreach my $aa (keys %cformat){
            $hformat{$aa}=$samp[$cformat{$aa}]; # print "$aa $hformat{$aa}"; <STDIN>;
        }
        $size=scalar(@alts);
        @dps=split(/,/,$hformat{AD});
        undef(@fin);
        undef(@find);

        for($i=0;$i<$size;$i++){
            if($dps[$i]==0){next;}
            if($h{$key} !~ /$alts[$i]/){next;}
            push(@fin,$alts[$i]);
            push(@find,$dps[$i]);
        }
        $ff=0;
        if(scalar(@fin)==1){$ff=1;$geno="$fin[0]$fin[0]"; }
        if(scalar(@fin)==2){$ff=1;$geno="$fin[0]$fin[1]"; }
        if(scalar(@fin)>2){
            $ff=1;
            $ct1=0;
            $max=0;
            $smax=0;
            $min=1000000000000;
            foreach my $cg (@find){
                if($cg>$max){$max=$cg;$mmx=$ct1;}
                if($cg<$min){$min=$cg;$mmy=$ct1;}
                $ct1++;
            }
            $ct1=0;
            foreach my $cg (@find){
                if($cg==$max and $mmx !~ $ct1){$smax=$cg;$mms=$ct1;}
                if($cg>$smax and $cg>=$min and $cg!=$max ){$smax=$cg;$mms=$ct1;}
                $ct1++;
            }
            $geno="$fin[$mmx]$fin[$mms]";
        }
        if($geno !~ /A|T|G|C|N/){print "WARNING sds <> <$geno> <@fin> <@find>\n"; $hall{$key}="NN";next;}
        if($ff==1){$hall{$key}="$geno";}
        #####END polimorphic
    } ###end IF exists
} ### end FIL
close(FIL);
open(LST,"$list") or die;


print "writing sites...\n";
$out="$Odir/$filen.tmpt";
#open (OUT, ">$out");
open (my $gout_fh, "| /bin/gzip -c > $out.gz") or die;


while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]";
    if(exists $hall{$key}){
        print $gout_fh "$hall{$key}\n";
    }
    if(!exists $hall{$key}){
        print $gout_fh "NN\n";
    }
}
close(LST);
close($gout_fh);
print "Final Done!!!\n";
