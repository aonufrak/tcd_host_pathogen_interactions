spack load py-htseq@0.11.2%gcc@8.4.1
for bam in ../analyses/02_star/geosmithia_morbida_mapping/*.out.bam
do
name=`basename $bam|sed s'/.Aligned.sortedByCoord.out.bam//g'`
echo Counting reads $name
htseq-count\
 --format=bam\
 --order=pos\
 --stranded=no\
 --type=gene\
 --idattr=ID\
 $bam\
 ../raw_data/GCF_012550715.1_ASM1255071v1_genomic.gff\
 >../analyses/03_htseq_count/g.morbida_counts/${name}.htseq.counts.txt
done
