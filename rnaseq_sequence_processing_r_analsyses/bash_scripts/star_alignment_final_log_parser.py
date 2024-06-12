import re
import sys
import os.path

summary_file=sys.argv[1]
with open(summary_file,'r') as summary:
	with open('star_alignment_parse_summary.txt','a') as parsed_summary:
		for line in summary.readlines():
			if re.search('Number of input reads',line):
				filtered_reads=re.sub('                          Number of input reads \|	','',line)
			if re.search('Uniquely mapped reads number',line):
				unique_map_reads=re.sub('                   Uniquely mapped reads number \|	','',line)
			if re.search('Uniquely mapped reads %',line):
				percent_unique_map_reads=re.sub('                        Uniquely mapped reads % \|	','',line)
			if re.search('Number of reads mapped to multiple loci',line):
				multiple_loci_mapped=re.sub('        Number of reads mapped to multiple loci \|	','',line)
			if re.search('% of reads mapped to multiple loci',line):
				percent_multiple_loci=re.sub('             % of reads mapped to multiple loci \|	','',line)
			if re.search('Number of reads mapped to too many loci',line):
				number_too_many_loci=re.sub('        Number of reads mapped to too many loci \|	','',line)
			if re.search('% of reads mapped to too many loci',line):
				percent_too_many_loci=re.sub('             % of reads mapped to too many loci \|	','',line)
			if re.search('Number of reads unmapped: too many mismatches',line):
				number_mismatches=re.sub('  Number of reads unmapped: too many mismatches \|	','',line)
			if re.search('% of reads unmapped: too many mismatches',line):
				percent_mismatches=re.sub('       % of reads unmapped: too many mismatches \|	','',line)
			if re.search('Number of reads unmapped: too short',line):
				number_short=re.sub('            Number of reads unmapped: too short \|	','',line)
			if re.search('% of reads unmapped: too short',line):
				percent_short=re.sub('                 % of reads unmapped: too short \|	','',line)
			if re.search('Number of reads unmapped: other',line):
				number_other=re.sub('                Number of reads unmapped: other \|	','',line)
			if re.search('% of reads unmapped: other',line):
				percent_other=re.sub('                     % of reads unmapped: other \|	','',line)
				print(filtered_reads.strip() + '\t' + unique_map_reads.strip() + '\t' +  percent_unique_map_reads.strip()+ '\t' + multiple_loci_mapped.strip()+  '\t' + percent_multiple_loci.strip()+ '\t' + number_too_many_loci.strip() + '\t' + percent_too_many_loci.strip() + '\t' + number_mismatches.strip() + '\t' + percent_mismatches.strip() + '\t' + number_short.strip() + '\t' + percent_short.strip()+ '\t' + number_other.strip() + '\t' + percent_other.strip() +  '\n')
