$file=$ARGV[0];
$list=$ARGV[1];


$bed="/path/to/reference/Zm-B73-REFERENCE-NAM-5.0.bed";
$ref="/path/to/reference/Zm-B73-REFERENCE-NAM-5.0.fa.OL";



$MaO=3;
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

#print "<$head>";<STDIN>;


@heads = split(/\s+/, $head);

#foreach my $aa(@heads){
#print "<$aa>";<STDIN>;
#}


#open(DB,"pigz -cdk $file |") or die;
open(DB,"gunzip -c $file |") or die;
print "Done!...\n";

print "loading bed...\n";
open(BED,"$bed") or die;
while($slin=<BED>) {
    chomp($slin);
    @lcol = split(/\t/, $slin);
    $hbed{$lcol[0]}=$lcol[2];
	print "$lcol[0] $lcol[2]\n";
    $ccll=$ccll+$lcol[2];
}
close(BED);
print "Done!...\n";

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

############################################################

print "Trasforming into mono & poli DB...\n";
$chr=1;
$pos=0;
$cp=0;
$achr=1;
$mono=0;
$NN=0;
### to start w chromosome 1
	$out="$file.$chr.poli.db";
		#$outm="$file.$chr.mono.db";
		open(OUT,">$out") or die;
		#open(OUTM,">$outm") or die;

while($line=<DB>) {
    $pos++; 
    if($pos>$hbed{$chr}){$chr++;$pos=1;}
	if($chr ne $achr){
        print "$achr...\n";
		close(OUT);
		#close(OUTM);
		$out="$file.$achr.poli.db";
		#$outm="$file.$achr.mono.db";
		
		#system ("sbatch /path/to/legacy/Caronte/code/Modules/RAND/GZIPER.sh $out");
		#system ("sbatch /path/to/legacy/Caronte/code/Modules/RAND/GZIPER.sh $outm");
		
		$out="$file.$chr.poli.db";
		#$outm="$file.$chr.mono.db";
		open(OUT,">$out") or die;
		#open(OUTM,">$outm") or die;
		print "$achr\n"; $achr=$chr
		}
    $ref=substr $href{$chr},$pos-1, 1;
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
	#print "<$c{A}><$c{C}> <$c{T}><$c{G}>  ref <$ref>"; <STDIN>;
 
    if($line=~ /A/ and $c{"A"} > $MaO){push(@alts, "A");}
    if($line=~ /C/ and $c{"C"} > $MaO){push(@alts, "C");}
    if($line=~ /T/ and $c{"T"} > $MaO){push(@alts, "T");}
    if($line=~ /G/ and $c{"G"} > $MaO){push(@alts, "G");}
    if($line=~ /\+/ and $c{"+"} > $MaO){push(@alts, "+");}
    if($line=~ /^\-/ and $c{"-"} > $MaO){push(@alts, "-");}
    if($line=~ /,\-/ and $c{"-"} > $MaO){push(@alts, "-");}
	
    $length = @alts;

    #print "aaa <$chr> <$achr> <$pos>  <$ref> $length <$line>"; <STDIN>;

#tell apart all N from ALL else
#if($length == 1 ){print OUTM "$chr\t$pos\t$ref\t$line";$mono++;}
if($length == 0){next;$NN++;}	
  if($length < 2 ){next;} #####main filter
  print OUT " $chr\t$pos\t$ref\t$line"; $cp++;
  $cc{$length}=$cc{$length}+1;
# print " bbb $pos PASS filter <$c{A}><$c{C}> <$c{T}><$c{G}> ref <$ref> $length <$cc{$length}>"; <STDIN>;
} #end while DB

### to end w last chromosome

close(OUT);
#close(OUTM);


		
print "total polimorphic sites  <$pos>   pass filter $MaO <$cp>\n";
foreach my $aa (sort keys %cc){
#print "<$aa> <$cc{$aa}>\n";<STDIN>;
}
print "monomorphic count <$mono>\n";
print "NNN count <$NN>\n";
print "END\n";
close (DB);

