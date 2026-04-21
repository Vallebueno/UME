$ref=$ARGV[0];
$file=$ARGV[1];
$bed=$ARGV[2];
$list=$ARGV[3];

$MaO=5;

#####################################################################################
print "Extracting headers from input list...\n";
open(LST,"$list") or die;
chomp(my @LLL = <LST>);
@lcol = split(/\|/, $LLL[0]);
$head=$lcol[0];
$head =~ s/\<\(zcat//g;
$head =~ s/\)//g;
$head =~ s/.vcf.gz.tmpt.gz//g;
$head =~ s/tmpt.gz//g;
$head =~ s/.sorted.MarkedDup.bam//g;
$head =~ s/paste  //g;
@heads = split(/\s+/, $head);
close(LST);

open(REF, $ref);
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

open(DB,"pigz -cdk $file |") or die;

#open (OUT, ">:gzip $file.Filt.gz" ) or die;
open (my $gzip_fh, "| /bin/gzip -c > $file.filt5.gz") or die;

$cc=0;
while($line=<DB>) {
$cc++;
    @lcol = split(/\t/, $line);
    $lcol[0] =~ s/ //;
    $chr=$lcol[0];
    $pos=$lcol[1];

    $ref=substr $href{$chr},$pos-1, 1;

  #  print "<$ref> <$lcol[0]> <$lcol[1]> <$lcol[2]> \n"; <STDIN>;



    undef %ex;
    undef (@alts);
    # push(@alts, $ref);
    undef %c;
    %c = ();

    if($line=~ /A/ ){$tst="A"; $c{$tst} = () = $line =~ /$test/g; }
    if($line=~ /C/ ){$tst="C"; $c{$tst} = () = $line =~ /$test/g;}
    if($line=~ /T/){$tst="T"; $c{$tst} = () = $line =~ /$test/g;}
    if($line=~ /G/ ){$tst="G"; $c{$tst} = () = $line =~ /$test/g; }
    if($line=~ /\+/ ){$tst="+"; $c{$tst} = () = $line =~ /$test/g; }
    if($line=~ /^\-/ or /,\-/){$tst="-"; $c{$tst} = () = $line =~ /$test/g;}

  #  print "count A<$c{A}> C<$c{C}> T<$c{T}> G<$c{G}>  ref <$ref>"; <STDIN>;

    if($line=~ /A/ and $c{"A"} > $MaO){push(@alts, "A");}
    if($line=~ /C/ and $c{"C"} > $MaO){push(@alts, "C");}
    if($line=~ /T/ and $c{"T"} > $MaO){push(@alts, "T");}
    if($line=~ /G/ and $c{"G"} > $MaO){push(@alts, "G");}
    if($line=~ /\+/ and $c{"+"} > $MaO){push(@alts, "+");}
    if($line=~ /^\-/ and $c{"-"} > $MaO){push(@alts, "-");}
    if($line=~ /,\-/ and $c{"-"} > $MaO){push(@alts, "-");}

    $length = @alts;


    if($length == 0){next;$NN++;}
    if($length < 2 ){next;} #####main filter
   # print "PASS $chr\t$pos\t$ref\t$line";<STDIN>; $cp++;

    $cp++;
    $pp=$cp/$cc;


#if($pos > 1000000){print "$chr\t$pos\t$ref\t<$pp>\n"; last}
    print $gzip_fh "$chr\t$pos\t$ref\t$line\n";

}

print "Done!...\n";
close($gzip_fh);




