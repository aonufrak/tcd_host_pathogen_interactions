for forward in ../raw_data/*1.fq.gz
do
name=`basename $forward|sed 's/_1.fq.gz//'`
reverse=`basename $forward|sed 's/_1.fq.gz/_2.fq.gz/'`
echo Trimming files $forward and ../raw_data/$reverse
fastp\
 -i $forward\
 -I ../raw_data/$reverse\
 -o ../analyses/01_fastp/${name}.r1.fastq\
 -O ../analyses/01_fastp/${name}.r2.fastq\
 --adapter_sequence AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT\
 --adapter_sequence_r2 GATCGGAAGAGCACACGTCTGAACTCCAGTCACGGATGACTATCTCGTATGCCGTCTTCTGCTTG\
 -g \
 -x \
 -w 16\
 -j ../analyses/01_fastp/${name}.report
done
