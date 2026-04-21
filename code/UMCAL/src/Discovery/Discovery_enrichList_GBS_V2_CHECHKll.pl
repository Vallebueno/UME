# UME makes the union of discovery DB and GBS DB (site summary file)
# Author: Miguel Vallebueno
# Date: 2022-07-29

#$filen=$ARGV[0];
#$chr=$ARGV[1];
#$cof=$ARGV[2];
#$alg=$ARGV[3];
#$qualf=$ARGV[4];
#$file="$filen.$chr.poli.qual$qualf.Union.V8.$alg.dbx.cof$cof.lst";
#$file =~ s/.db.gz//g;
#$file = $filen;
use List::Util qw(uniq);
$file="/scratch-cbe/users/miguel.vallebueno/Production_ALL_Data_2023_10/Discovery_2023_08_cof10_ALLC_V2.ll";

$GBSf="/groups/swarts/lab/MAVE/Bioinfo/Diversity_panels/GBS/Maize_GBS_Kriztian_build.lst";

$out="/scratch-cbe/users/miguel.vallebueno/Production_ALL_Data_2023_10/Discovery_2023_08_cof10_ALLC_V3.ll";


$ref="/groups/swarts/lab/Resources/ReferenceGenomes/Zm-B73-REFERENCE-NAM-5.0.fa.OL";
print "loading ref...\n";
open(REF,"$ref") or die;
chomp(my @re = <REF>); #magic
$f=0;
foreach my $aa(@re){
    if($f==1){$href{$name}=$aa;$f=0;next;}
    if($f==0){$name=$aa; $name=~ s/>//; print "$name\n";$f=1;$ctrl{$name}=1;next;}
}
close(REF);
undef(@re);
print "Done!...\n";


open(OUT,">$out") or die;
print "Loading GBS sites...\n";
open(GBS,"$GBSf") or die;

while($line=<GBS>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $key="$lcol[0]&$lcol[1]";
    $sys=$lcol[2];
    $hgbs{$key}=$sys;
    if($lcol[0] eq "Chromosome"){next;}
    $all{$lcol[0]}{$lcol[1]}=0;

}
close(GBS);

print "Loading Discovery List sites...\n";
open(FIL,"$file") or die;
while($line=<FIL>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    $lcol[0] =~ s/ //g;
    $key = "$lcol[0]&$lcol[1]";
    $hdb{$key} = $line;
    $all{$lcol[0]}{$lcol[1]}=1;
    if($lcol[4] eq ""){$all{$lcol[0]}{$lcol[1]}=0;}

}
close (FIL);

    print "Unifiying sites...\n";
    foreach my $chr (sort { $a <=> $b } keys %all){
      #  print "<<$chr>>";<STDIN>;
        foreach my $pos (sort { $a <=> $b } keys %{$all{$chr}}) {
            # print "<<$chr>>  <$pos>";<STDIN>;
            $key = "$chr&$pos";
            if (!exists $hgbs{$key} and $all{$chr}{$pos} == 1) {
                print OUT "$hdb{$key}\n";
                next;
            }
            if (!exists $hgbs{$key} and $all{$chr}{$pos} == 0) {next;}
            if (!exists $hdb{$key}) {
                $gref = substr $href{$chr}, $pos - 1, 1;

                @agbs = split(/ /, $hgbs{$key});

                undef(@array);
                foreach my $aa (@agbs) {
                    if ($aa ne $gref) {push @array, $aa}
                }
                $sys = join(",", @array);
                print OUT "$chr\t$pos\t.\t$gref\t$sys\t.\tPASS\t.\tGT\n";
                #print "nex $chr\t$pos\t.\t$gref\t$sys\t.\tPASS\t.\tGT\n";
                next;
            }

            if (exists $hgbs{$key} and exists $hdb{$key}) {
                $gref = substr $href{$chr}, $pos - 1, 1;

                @agbs  = split(/ /, $hgbs{$key});

                undef(@adb);
                @lcol  = split(/\t/, $hdb{$key});
                @xcol  = split(/,/, $lcol[4]);
                push(@adb,$gref);
                push(@adb,@xcol);


                @a=sort(@adb);
                @b=sort(@agbs);


                if (@adb eq @agbs){print OUT "$hdb{$key}\n"; next;}
                $gref = substr $href{$chr}, $pos - 1, 1;

                undef(@sums);
                push (@sums, @adb);
                push (@sums, @agbs);

                @uniques= uniq(@sums);

                #print "tyt <@adb>   <@agbs>   <@uniques>\n"; <STDIN>;

                undef(@array);
                foreach my $aa (@uniques) {
                    if ($aa ne $gref) {push @array, $aa}
                }
                $sys = join(",", @array);
                print OUT "$chr\t$pos\t.\t$gref\t$sys\t.\tPASS\t.\tGT\n";
               # print "tyt <@adb>   <@agbs>   <@uniques>\n"; <STDIN>;
                #print  "$chr\t$pos\t.\t$gref\t$sys\t.\tPASS\t.\tGT\n";<STDIN>;




                #print OUT "$hdb{$key}\n";

                #@adb   = split(/ /, $lcol[4]);
                #@agbs  = split(/ /, $hgbs{$key});

                #print " <$lcol[3]> <@adb>   <@agbs>"; <STDIN>;}

            }
        }}
    close(OUT);



