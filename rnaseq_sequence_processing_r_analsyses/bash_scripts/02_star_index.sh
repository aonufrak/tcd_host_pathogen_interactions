/pickett_centaur/software/STAR/STAR \
--genomeSAindexNbases 11 \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /pickett_centaur/project/tcd_control_rnaseq/analyses/2_star/genome_directory_gmorbida \
--genomeFastaFiles /pickett_centaur/project/tcd_control_rnaseq/raw_data/GCF_012550715.1_ASM1255071v1_genomic.fna \
--genomeChrBinNbits=10 \
--sjdbGTFfile /pickett_centaur/project/tcd_control_rnaseq/raw_data/GCF_012550715.1_ASM1255071v1_genomic.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 149
