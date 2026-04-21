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

#$chr=$ARGV[0];


#$lfile="/scratch-cbe/users/miguel.vallebueno/MERGE/files_$chr.lst";
$lfile=$ARGV[0];
$cf=$ARGV[1];

$out=$lfile;
$out =~ s/.lst/.$cf.lst2/g;


open(OUT,">$out") or die;

open(LST,"$lfile") or die;

chomp(my @re = <LST>);

foreach my $aa (@re){

    $aa =~ s/cof5/cof$cf/g;

    print "Loading Discovery List sites <$aa>...\n";
    open(FIL,"$aa") or die;
$ffx=0;
    while($line=<FIL>) {
        chomp($line);
        @lcol = split(/\t/, $line);
        $chr=$lcol[0];
        if($ffx==0){$achr=$chr;$ffx=1}
        if($chr ne $achr){print"FATAL:43fj"; die;}
        @co = split(/ /, $lcol[2]);
        foreach my $aa (@co){
            $h{$lcol[1]}{$aa}=1;
        }
    }
    close (FIL);
}

foreach my $aa (sort { $a <=> $b } keys %h){
    undef(@pa);
    foreach my $bb (keys  %{$h{$aa}}){
        push(@pa,$bb);
    }

    print OUT "$chr\t$aa\t@pa\n";
}

close(OUT);