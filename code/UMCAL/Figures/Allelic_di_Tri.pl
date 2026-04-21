
# Author: Miguel Vallebueno
# Date: 2022-07-29


open(OUT,">cuant_di_tri.txt") or die;

for ($i=1;$i<11;$i++){


    $file="1.ALL.Merge.$i.poli.qual0.Union.V8.max.dbx.cof5.lst2";



open(FIL,"$file") or die;

while($line=<FIL>) {
    chomp($line);
    @lcol = split(/\t/, $line);
    @col = split(/\s/, $lcol[2]);
   $tam=scalar(@col);
    if($tam==1){$h{$i}{1}++;}
    if($tam==2){$h{$i}{2}++;}
    if($tam==3){$h{$i}{3}++;}
    if($tam==4){$h{$i}{4}++;}
    if($tam==5){$h{$i}{5}++;}
    if($tam==6){$h{$i}{6}++;}
    if($tam>6){$h{$i}{7}++;}

}
close (FIL);

foreach my $aa (sort keys %h){

    foreach my $bb (sort keys %{ $h{$aa} } ){

        print OUT "$aa\t$bb\t$h{$aa}{$bb}\n";

    } #END AA
    }   #END BB
} #end for
close(OUT);