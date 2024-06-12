touch /pickett_centaur/project/tcd_control_rnaseq/analyses/02_star/juglans_nigra_mapping/star_juglans_alignment_parse_summary.txt
echo 'Tree Total_Reads Unique_Mapped Percent_Unique Multiple_Loci Percent_Multi_Loci Many_Loci Percent_Many_Loci Many_Mismatches Percent_Many_Mismatches Short Percent_Short Other_Unmapped Percent_Other_Unmapped' >> /pickett_centaur/project/tcd_control_rnaseq/analyses/2_star/juglans_nigra_mapping/star_juglans_alignment_parse_summary.txt

for file in /pickett_centaur/project/tcd_control_rnaseq/analyses/02_star/juglans_nigra_mapping/*jnigra.Log.final.out
do
base=`basename $file | sed s/.Log.final.out//`
echo $base 
seqstats=`python ./star_alignment_final_log_parser.py $file`
echo $base $seqstats >> /pickett_centaur/project/tcd_control_rnaseq/analyses/02_star/juglans_nigra_mapping/star_juglans_alignment_parse_summary.txt
done
