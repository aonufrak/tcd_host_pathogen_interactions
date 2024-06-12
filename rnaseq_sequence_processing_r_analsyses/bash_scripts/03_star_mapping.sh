for forward in ../analyses/01_fastp/*r1.fastq
do
name=`basename $forward|sed 's/r1.fastq//'`
reverse=`basename $forward|sed 's/r1.fastq/r2.fastq/'`
echo Mapping files $forward and $reverse
/pickett_centaur/software/STAR/STAR\
 --runThreadN 10\
 --genomeDir /pickett_centaur/project/tcd_control_rnaseq/analyses/02_star/genome_directory_gmorbida\
 --readFilesIn $forward /pickett_centaur/project/tcd_control_rnaseq/analyses/01_fastp/$reverse\
 --outSAMtype BAM SortedByCoordinate\
 --outSAMunmapped Within\
 --outFileNamePrefix ../analyses/02_star/geosmithia_morbida_mapping/${name}geosmithia.
done
