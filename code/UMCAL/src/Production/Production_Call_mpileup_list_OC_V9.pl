# UME makes the production call from mpileups generating one column files
# guided by the discovery .ll file. The resulting per-sample columns are merged
# later into the final multi-sample production output.
# Author: Miguel Vallebueno
# Date: 2023-01-10

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
# exit;
}


open(LST,"$list") or die;
print "Loading sites...\n";
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]";
    @col = split(/,/, $lcol[4]);
    undef(@ff);
    push(@ff,$lcol[3]);
    push(@ff,@col);
    $sys=join(" ", @ff);
    $h{$key}=$sys;
}
close(LST);
print "Done!!!\n";


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

print "reading mpileup <$file>...\n";
open(FIL,"gunzip -c $file |") or die;
print "Checkpoint 1 <$file>...\n";
#if($chr !~ /^10/){next;}
while($line=<FIL>){
    $dp=0;
    chomp($line);
    if($line =~ /^#/){next;}
    @col = split(/\t/, $line);
    $chr=$col[0];
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
                $hall{$key}="0/0:$mdp,0";
            }
        }
        next;
    } #end if SNR END
    if(exists $h{$key}) {
        #weird exceptions without AD
        if ($info =~ /DP=0/ ) {next;}
        elsif ($info !~ /AD/ and $col[4] eq "<*>") {
            @xaa=split(/DP=/,$info);
            @xbb=split(/;/,$xaa[1]);
            $mdp=$xbb[0];
            $hall{$key}="0/0:$mdp,0";
            next;
        }
        elsif ($info =~ /AD=0;/ and $col[4] eq "<*>") {
            @xaa=split(/DP=/,$info);
            @xbb=split(/;/,$xaa[1]);
            $mdp=$xbb[0];
            $hall{$key}="0/0:$mdp,0";
            next;
        }
        elsif ($info =~ /AD=0,0;/ and $col[4] eq "<*>") {
            @xaa=split(/DP=/,$info);
            @xbb=split(/;/,$xaa[1]);
            $mdp=$xbb[0];
            $hall{$key}="0/0:$mdp,0";
            next;
        }
        #polimorphic
        @altx = split(/ /, $h{$key});
        $grc=0;
        foreach my $aa (@altx){
            $genor{$aa}=$grc;
            $grc++;
        }
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
        # Extract info from format
        foreach my $aa (keys %cformat){
            $hformat{$aa}=$samp[$cformat{$aa}];
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
        if(scalar(@fin)==1){
            $val=$genor{$fin[0]};
            if($val == 0){$hall{$key}="0/0:$find[0],0,0,0,0,0";next;}
            if($val == 1){$hall{$key}="1/1:0,$find[0],0,0,0,0";next;}
            if($val == 2){$hall{$key}="2/2:0,0,$find[0],0,0,0";next;}
            if($val == 3){$hall{$key}="3/3:0,0,0,$find[0],0,0";next;}
            print "WARNING 2o3u4o23u4\n";
        }
        elsif(scalar(@fin)==2){
            $val1=$genor{$fin[0]};
            $val2=$genor{$fin[1]};
            @ADAS=(0,0,0,0,0,0);
            $ADAS[$val1]=$find[0];
            $ADAS[$val2]=$find[1];
            $sys=join(",", @ADAS);
            if($val1<$val2){$hall{$key}="$val1/$val2:$sys";next;}
            $hall{$key}="$val2/$val1:$sys";
            next;
        }
        elsif(scalar(@fin)==0){next;}
        $ddf=scalar(@fin);
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
        $val1=$genor{$fin[$mmx]};
        $val2=$genor{$fin[$mms]};
        @ADAS=(0,0,0,0,0,0);
        $ADAS[$val1]=$find[$mmx];
        $ADAS[$val2]=$find[$mms];
        $sys=join(",", @ADAS);
        if($val1<$val2){$hall{$key}="$val1/$val2:$sys";next;}
        $hall{$key}="$val2/$val1:$sys";
        #####END polimorphic
    } ###end IF exists
} ### end FIL
close(FIL);
open(LST,"$list") or die;
print "writing sites...\n";
$out="$Odir/$filen.tmpt";
#open (OUT, ">$out");
open (my $gout_fh, "| gzip -c > $out.gz") or die;
while($line=<LST>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0]=~ s/ //;
    $key="$lcol[0]&$lcol[1]";
    if(!exists $hall{$key}){print $gout_fh "./.\n";next;}
    $geno=$hall{$key};
    print $gout_fh "$geno\n";
}
close(LST);
close($gout_fh);
print "Final Done!!!\n";
#print "$line\n";
#print "zzzz lst<@altx>   det<@fin> AD<@find>   tran<$val>  value<$fin[0]>\n";
#print "despues key<$key>  indgeno<$geno>  alts<@alts>\n";<STDIN>;

#print "antes key<$key>  indgeno<$geno>  alts<@alts>\n";<STDIN>;
