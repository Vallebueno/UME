file=MasterLandraceTeoInbredGBS_collapseDist0.02.vcf
awk '!/^#/' $file |cut -f 1,2,4,5 | awk '{ nbed=$1"\t"$2"\t.\t"$3"\t"$4"\t.\tPASS\t.\tGT"; print nbed }' | sed 's/,N//g' | sed 's/N//g'> $file.pos.ll

cat $file.pos | sed 's/,N//g' | sed 's/N//g'> $file.pos.ll


1      17703032



