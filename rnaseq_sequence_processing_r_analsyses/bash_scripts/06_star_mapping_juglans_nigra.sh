for forward in /pickett_centaur/project/tcd_control_rnaseq/analyses/02_star/gm_unmapped_fastq/*1.fastq
do
name=`basename $forward|sed 's/.gm.unmapped.1.fastq//'`
reverse=`basename $forward|sed 's/.1.fastq/.2.fastq/'`
echo Mapping files $forward and $reverse
/pickett_centaur/software/STAR/STAR\
 --runThreadN 10\
 --genomeDir /pickett_centaur/project/tcd_control_rnaseq/analyses/02_star/genome_directory_juglans_nigra\
 --readFilesIn $forward /pickett_centaur/project/tcd_control_rnaseq/analyses/02_star/gm_unmapped_fastq/$reverse\
 --outSAMtype BAM SortedByCoordinate\
 --outSAMunmapped Within\
 --outFileNamePrefix ../analyses/02_star/juglans_nigra_mapping/${name}.jnigra.
done
