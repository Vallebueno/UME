$file=$ARGV[0];
$name=$ARGV[1];
$lt=$ARGV[2];
$ref=$ARGV[3];



$file="/groups/swarts/lab/MAVE/Experiments/Hydra_coverage/OL/files.mergein3.lst.tlone.Merge.db.7.poli.qual0.Union.V7-1.max.db.lst";
$ref="/groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL";

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

open(LST,"$file") or die;
print "Loading sites...\n";
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
$key="$lcol[0]&$lcol[1]";

    #  print "<$lcol[0]><$lcol[1]><$lcol[2]>";<STDIN>;
    $h{$key}=$lcol[2];
}
close(LST);
print "DONE!!!!!\n";

$file="/scratch-cbe/users/miguel.vallebueno/Downsampling/mpileup.lst";

open(LST,"$file") or die;

while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    push(@indsx,$lcol[1]);
    print "loading <$lcol[1]> $lcol[0]\n";
    open(FIL,"$lcol[0]") or die;
    while($line=<FIL>){
        chomp($line);
        if($line =~ /^#/){next;}
       @col = split(/\t/, $line);
        $chr=$col[0];
        if($chr !~ /^7/){next;}
        $pos=$col[1];
        $key="$col[0]&$col[1]";
        $info=$col[7];
#print "key <$key> <$info>";<STDIN>;
        if($info =~ /END=/){
            @aaa=split(/END=/,$info);
            @bbb=split(/;/,$aaa[1]);
            $END=$bbb[0];
            for($i=$pos;$i<=$END;$i++){
                $key="$chr&$i";
                if(exists $h{$key}){
                    $gref=substr $href{$chr},$i-1, 1;
                    if($h{$key} !~ /$gref/){next;}
                    $geno=$gref.$gref;
                    $hall{$key}{$lcol[1]}=$geno;
                }
            }
            next;
        } #end if SNR END
        if($info !~ /AD/ and $col[4] eq "<*>"){
            if(exists $h{$key}){
                $geno=$col[3].$col[3];
                $hall{$key}{$lcol[1]}=$geno;
            }
            next;
        }
        if($info !~ /AD/){print "WARNING: 345367 no AD <$line>\n";next; }
        undef(@alts);
        if(exists $h{$key}){
            push(@alts,$col[3]);
            #print "exists <$h{$key}>\n";<STDIN>;
            @al=split(/,/,$col[4]);
            push(@alts,@al);
          #  print "<@alts><@al>";<STDIN>;
            undef(%hformat);
            @format=split(/:/,$col[8]);
            @samp=split(/:/,$col[9]);
            if(@format ne @aformat){
                $cc=0;
               # print "1) @format   @aformat\n";<STDIN>;
                foreach my $aa (@format){
                    $cformat{$aa}=$cc;
                   #   print "2) $aa $cformat{$aa}\n";<STDIN>;
                    $cc++;
                }
            }
            @aformat=@format;
            ###### Extract info from format
            foreach my $aa (keys %cformat){
                $hformat{$aa}=$samp[$cformat{$aa}]; # print "$aa $hformat{$aa}"; <STDIN>;
            }
             #   print "yyy <$key> <$hformat{AD}> <$line>\n";<STDIN>;
            $size=scalar(@alts);
            @dps=split(/,/,$hformat{AD});
            undef(@fin);
            for($i=0;$i<$size;$i++){
                    #print "xxxx   <$key> <@alts><$hformat{AD}>   <$i> <$alts[$i]> <$dps[$i]>\n";<STDIN>;
                    if($dps[$i]==0){next;}
                    if($h{$key} !~ /$alts[$i]/){next;}
                    push(@fin,$alts[$i]);
                    #print "zzz   <$h{$key}><$alts[$i]>   <@fin>\n";<STDIN>;
                }
            my $color = $fin[ rand @fin ];
            $geno=$color.$color;
            #print "<$key> <$h{$key}> <@alts>  <$hformat{AD}>  <@fin>  chs<$geno>\n";<STDIN>;
            #if($size == 1){$geno=$fin[0].$fin[0]}
            $hall{$key}{$lcol[1]}=$geno;
           # $hall{$key}{$lcol[1]}=;
        } # end if exists
    }
    close(FIL);
}
print "printing hmp file....\n";
open(OUT,">$file.hmp") or die;


$gg=join "\t", @indsx;
print OUT "rs#\talleles\tchrom\tpos\tstrand\tassembly#\tcenter\tprotLSID\tassayLSID\tpanelLSID\tQCcode\t$gg\n";
$cc=0;
foreach my $aa (sort keys %hall){
    $cc++;
    @ii=split(/&/,$aa);
    undef(@alts);
    undef(@geno);
    $fa=0;
    $fc=0;
    $ft=0;
    $fg=0;
    $fI=0;
    foreach my $bb (@indsx){

        $gen=$hall{$aa}{$bb};
        if(!exists $hall{$aa}{$bb}){$gen="NN"};
        push(@geno,$gen);

        if($gen =~ /A/ and $fa==0){$fa=1; push(@alts, "A")}
        if($gen =~ /T/ and $fc==0){$fc=1; push(@alts, "T")}
        if($gen =~ /C/ and $ft==0){$ft=1; push(@alts, "C")}
        if($gen =~ /G/ and $fg==0){$fg=1; push(@alts, "G")}
    }
    $ll=scalar(@alts);
    if($ll != 2){next}

    $gg=join "\t", @geno;
    $sys=join "/", @alts;
    print OUT "$cc\t$sys\t$ii[0]\t$ii[1]\t+\tAGPv5\tGMI\t\t\tTEST_Panel\t\t$gg\n";

}

print "Final Done!!!\n";