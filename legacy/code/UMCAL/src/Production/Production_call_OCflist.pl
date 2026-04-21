$ref="/path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL";
$file="/path/to/workdir/Downsampling2/SAMN06685743_RI_DeNovo_PHG35_mappedAGPv5.I_DSN0.5x.bam_Samtools/SAMN06685743_RI_DeNovo_PHG35_mappedAGPv5.I_DSN0.5x.bam.7.mpileup";
$list="/path/to/project_inputs/files.mergein3.lst.tlone.Merge.db.7.poli.qual0.Union.V7-1.max.db.lst";

#$file=$ARGV[0];
#$ref=$ARGV[1];
#$list=$ARGV[2];

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
   # print "<$key>\n"; <STDIN>;
}
close(LST);


print "reading mpileup <$file>...\n";
    open(FIL,"$file") or die;

$out="$file.ol";
open (OUT, ">$out");
    while($line=<FIL>){
        $dp=0;
        chomp($line);
        if($line =~ /^#/){next;}
       @col = split(/\t/, $line);
        $chr=$col[0];
        if($chr !~ /^7/){next;}
        if($line =~/INDEL/){next;}
        $pos=$col[1];
        $key="$col[0]&$col[1]";
        $info=$col[7];
#print "key <$key> <$info> \n";<STDIN>;
        if($info =~ /END=/){
            @aaa=split(/END=/,$info);
            @bbb=split(/;/,$aaa[1]);
            $END=$bbb[0];

            @xaa=split(/MinDP=/,$info);
            @xbb=split(/;/,$xaa[1]);
            $mdp=$xbb[0];

            for($i=$pos;$i<=$END;$i++){
                $geno="XX";
                $key="$chr&$i";

                if(exists $h{$key}){
                    #print "exists <$key> <$h{$key}> <>\n"; <STDIN>;

                    $gref=substr $href{$chr},$i-1, 1;
                    if($h{$key} !~ /$gref/){next;}
                    $geno=$gref.$gref;
                    $hall{$key}="$geno:$mdp,0";
                    #print "$hall{$key}\n";<STDIN>;
                    #if($hall{$key} !~ /A|T|G|C|N/){print "WARNING:1111  <$line>\n";<STDIN>;
                }
            }
            next;
        } #end if SNR END


        ####weird exceptions
        if(exists $h{$key}) {

            if ($info !~ /AD/ and $col[4] eq "<*>") {
                if (exists $h{$key}) {
                    $geno = $col[3].$col[3];
                    $hall{$key}="$geno:1,0";
                    #if ($hall{$key} !~ /A|T|G|C|N/) {print "WARNING:2222  <$line>\n";<STDIN>;}
                }
                next;
            }

            if ($info =~ /AD=0/ and $col[4] eq "<*>") {


                $geno = $col[3].$col[3];
                $ $hall{$key}="$geno:1,0";

                next;
                #if ($hall{$key} !~ /A|T|G|C|N/) {print "WARNING:555  <$line>\n";<STDIN>;}

            }
            if($info !~ /AD/){
                $geno = $col[3].$col[3];
                $hall{$key}="$geno:1,0";
                next;
                #if ($hall{$key} !~ /A|T|G|C|N/) {print "WARNING: 345367 no AD <$line>\n";<STDIN>;}
            }

        } #end exists weird exceptions

        #######################

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
            undef(@find);
            for($i=0;$i<$size;$i++){
                    #print "xxxx   <$key> <@alts><$hformat{AD}>   <$i> <$alts[$i]> <$dps[$i]>\n";<STDIN>;
                    if($dps[$i]==0){next;}
                if($dps[$i]>$dp){$dp=$dps[$i]}
                    if($h{$key} !~ /$alts[$i]/){next;}
                    push(@fin,$alts[$i]);
                    push(@find,$dps[$i]);
                   # print "zzz   <$h{$key}><$alts[$i]>   <@fin> <@find>\n";<STDIN>;
                }
            if(scalar(@fin)==1){$geno="$fin[0]$fin[0]:$find[0],0";}
            if(scalar(@fin)==2){$geno="$fin[0]$fin[1]:$find[0],$find[1]";}
            if(scalar(@fin)==3){
                if($find[0]>$find[2] and $find[1]>$find[2]){$geno="$fin[0]$fin[1]:$find[0],$find[1]";}
                if($find[0]>$find[1] and $find[2]>$find[1]){$geno="$fin[0]$fin[2]:$find[0],$find[2]";}
                if($find[1]>$find[0] and $find[1]>$find[0]){$geno="$fin[1]$fin[2]:$find[1],$find[2]";}
            }
            if(scalar(@fin)==4){$geno="$fin[0]$fin[1]:$find[0],$find[1]";print "Warning 42o4854\n"}
            if($geno !~ /A|T|G|C|N/){print "WARNING sds <$geno> <@fin>";<STDIN>; $hall{$key}="NN:0,0";next;}
            $hall{$key}="$geno";

        } # end if exists
    }


open(LST,"$list") or die;
print "writing sites...\n";
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]";

    if(exists $hall{$key}){
        #if($hall{$key} !~ /A|T|G|C|N/){print "WARNING:444  <$key><$hall{$key}>\n";<STDIN>;}
           print OUT "$hall{$key}\n";
    }

    if(!exists $hall{$key}){
        print OUT "NN:0,0\n";
    }

}
close(LST);


print "Final Done!!!\n";
