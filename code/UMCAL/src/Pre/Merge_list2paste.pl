
$KEY=$ARGV[0];

open(KK,"$KEY") or die;
chomp(my @lk = <KK>);
$out="$KEY.tlone";

open(OUT,">$out") or die;
$step=3;
$size=3;

$esize=$size+$step;
$cc=1;

for ($i=$step; $i <= $esize; $i=$i+$step){
        $ff=$i-$step+1;

        print "<$cc>\n";

    undef(@spl);

        for ($j=$ff; $j<=$i; $j++) {
            if ($j == $size+1){last;}
            $ff=$lk[$j-1];
            $ffil="<(zcat $ff.tmpt.gz)";

            push(@spl, $ffil);
        }

    $all=join " ", @spl;

    print OUT "paste $all > $cc.ALL.Merge.db\n";$cc++;


#system ("sbatch $caronte/Modules/VCF/Merge_columns2DB.sh $all");

}
