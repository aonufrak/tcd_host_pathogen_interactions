/pickett_centaur/software/STAR/STAR \
--genomeSAindexNbases 11 \
--runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /pickett_centaur/project/tcd_control_rnaseq/analyses/2_star/genome_directory_juglans_nigra \
--genomeFastaFiles /pickett_centaur/project/tcd_control_rnaseq/raw_data/jni1_wg.fa \
--genomeChrBinNbits=10 \
--sjdbGTFfile /pickett_centaur/project/tcd_control_rnaseq/raw_data/jni1_wg.gff \
--sjdbGTFtagExonParentTranscript Parent \
--sjdbOverhang 149

