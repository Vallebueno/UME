# UME makes the production call from mpileups generating one column files guided by the Discovery list file. need to be joined afterwise into vcf
# Author: Miguel Vallebueno
# Date: 2022-11-08

#$ref="/groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL";
#$file="/scratch-cbe/users/miguel.vallebueno/mpilp2db/SM8_NoUDG_Vallebueno_SM8_mappedAGPv5.sorted.MarkedDup.bam.ALL.mpileup.gz";
#$list="/scratch-cbe/users/miguel.vallebueno/mpilp2db/1.ALL.Merge.10.poli.qual0.Union.V8.max.dbx.cof5.lst";


$file=$ARGV[0];
$ref=$ARGV[1];
$list=$ARGV[2];
$Odir=$ARGV[3];

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

    $gref=substr $href{$lcol[0]},$lcol[1]-1, 1;
   #  print "<$lcol[0]>  <$lcol[1]>    <$key><$h{$key}>   <$gref>\n"; <STDIN>;
}
close(LST);
print "Done!!!\n";


print "reading mpileup <$file>...\n";

open(FIL,"gunzip -c $file |") or die;

$out="$file.ol";

open (OUT, ">$out");
while($line=<FIL>){
    $dp=0;
    chomp($line);
    if($line =~ /^#/){next;}
    @col = split(/\t/, $line);
    $chr=$col[0];
    if($chr !~ /^10/){next;}
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
            $key="$chr&$i";
            if(exists $h{$key}){
                #print "exists <$key> <$h{$key}> <>\n"; <STDIN>;
                $gref=substr $href{$chr},$i-1, 1;
                if($h{$key} !~ /$gref/){next;}
                $geno=$gref.$gref;
                $hall{$key}="$geno:$mdp,0";
                #    print "res ref <$hall{$key}>\n";<STDIN>;
                #if($hall{$key} !~ /A|T|G|C|N/){print "WARNING:1111  <$line>\n";<STDIN>;
            }
        }
        next;
    } #end if SNR END



    if(exists $h{$key}) {
        # print "key <$key> <$info> <$h{$key}>  <$line>\n";<STDIN>;

        ####weird exceptions without AD
        if ($info !~ /AD/ and $col[4] eq "<*>") {
            if (exists $h{$key}) {
                $geno = $col[3].$col[3];
                $hall{$key}="$geno:1,0";
                #   print "AA1<$hall{$key}>\n";<STDIN>;
            }
            next;
        }

        if ($info =~ /AD=0/ and $col[4] eq "<*>") {
            $geno = $col[3].$col[3];
            $hall{$key}="$geno:1,0";
            #  print "AA2<$hall{$key}>\n";<STDIN>;
            next;
        }

        if($info !~ /AD/){
            $geno = $col[3].$col[3];
            $hall{$key}="$geno:1,0";
            # print "AA3<$hall{$key}>\n";<STDIN>;
            next;
        }
        #end exists weird exceptions


        ###############polimorphic
        undef(@alts);

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
        #   print "yyy <$key> <$hformat{AD}> <$line>\n";<STDIN>;
        $size=scalar(@alts);
        @dps=split(/,/,$hformat{AD});
        undef(@fin);
        undef(@find);

        for($i=0;$i<$size;$i++){
            #print "xxxx   <$key> <@alts><$hformat{AD}>   <$i> <$alts[$i]> <$dps[$i]>\n";<STDIN>;
            if($dps[$i]==0){next;}
            #if($dps[$i]>$dp){$dp=$dps[$i]}
            if($h{$key} !~ /$alts[$i]/){next;}
            push(@fin,$alts[$i]);
            push(@find,$dps[$i]);
            # print "zzz   <$h{$key}><$alts[$i]>   <@fin> <@find>\n";<STDIN>;
        }

        $ff=0;
        if(scalar(@fin)==1){$ff=1;$geno="$fin[0]$fin[0]:$find[0],0"; }
        if(scalar(@fin)==2){$ff=1;$geno="$fin[0]$fin[1]:$find[0],$find[1]"; }
        if(scalar(@fin)>2){
            $ff=1;
            $ct1=0;
            $max=0;
            $smax=0;
            $min=100000000000;
            foreach my $cg (@find){
                if($cg>$max){$max=$cg;$mmx=$ct1;}
                if($cg<$min){$min=$cg;$mmy=$ct1;}
                #   print "d1 <$ct1>  <$cg>   max<$max> mmx<$mmx>  min<$min> mmy<$mmy>"; <STDIN>;
                $ct1++;
            }

            $ct1=0;

            foreach my $cg (@find){
                if($cg==$max and $mmx !~ $ct1){$smax=$cg;$mms=$ct1;}
                if($cg>$smax and $cg>=$min and $cg!=$max ){$smax=$cg;$mms=$ct1;}
                #  print "d2 <$ct1>  <$cg>   max<$max> mmx<$mmx>  min<$min> mmy<$mmy>    smax<$smax>  mms<$mms>"; <STDIN>;
                $ct1++;
            }
            $geno="$fin[$mmx]$fin[$mms]:$find[$mmx],$find[$mms]";
            #print "Warning 42o4854  <$h{$key}> <$geno> <@fin> <@find>   max<$max>  smax<$smax>  <$geno>\n";<STDIN>;
        }
        if($geno !~ /A|T|G|C|N/){print "WARNING sds <> <$geno> <@fin> <@find>";<STDIN>; $hall{$key}="NN:0,0";next;}

        if($ff==1){$hall{$key}="$geno";}

        #####END polimorphic

    } ###end IF exists



} ### end FIL
close(FIL);
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