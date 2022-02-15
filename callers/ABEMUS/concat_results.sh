

for i in /data/dilutions_chr22/*sorted.bam ;

do echo $i ;
mkdir /data/abemus_outdir/abemus_outdir_chr22/$(basename $i .bam)

for j in {2..3} ;
	do ls /data/abemus_outdir/abemus_outdir_chr22_*/Results/$(basename $i .bam)/pmtab_F${j}_*.tsv
awk '
    FNR==1 && NR!=1 { while (/^<header>/) getline; }
    1 {print}
' /data/abemus_outdir/abemus_outdir_chr22_*/Results/$(basename $i .bam)/pmtab_F${j}_*.tsv  >/data/abemus_outdir/abemus_outdir_chr22/$(basename $i .bam)/pmtab_F${j}_$(basename $i .bam).tsv
echo /data/abemus_outdir/abemus_outdir_chr22/$(basename $i .bam)/pmtab_F${j}_$(basename $i .bam).tsv

done
done
