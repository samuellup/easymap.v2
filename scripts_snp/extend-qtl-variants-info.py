import argparse

# Parse command arguments
parser = argparse.ArgumentParser()
parser.add_argument('-cand', action="store", dest='cand_list', required=True)
parser.add_argument('-af-info', action="store", dest='af_info', required=True)
parser.add_argument('-out', action="store", dest='output', required=True)
args = parser.parse_args()

dic_af = dict()
with open(args.af_info, "r") as info: 
	for line in info: 
		sp = line.split("\t")
		va_id = str(sp[1]).lower() + "-" + str(sp[0]).lower()
		dic_af[va_id] = [sp[9], sp[10], sp[11]]

with open(args.output, "w") as out: 
	out.write('@type\tcontig\tposition\tref_base\talt_base\talt_quality\tref_test_count\talt_test_count\talt_allele_freq_test\talt_allele_freq_control\tdAF\tdist_to_qtl_center\thit\tmrna_start\tmrna_end\tstrand\tgene_model\tgene_element\taa_pos\taa_ref\taa_alt\tgene_funct_annot\tf_primer\ttm_f_primer\tr_primer\ttm_r_primer\tupstream_sequence\tdownstream_sequence\n')
	with open(args.cand_list, "r") as candidates: 
		for line in candidates: 
			if not line.startswith("@"):
				sp = line.split("\t")
				if len(str(sp[3])) == 1 and len(str(sp[4])) == 1: mut_type = "snp"
				else: mut_type = "indel"
				pos, cont = str(sp[2]).lower(), str(sp[1]).lower()
				va_id = pos + "-" + cont
				try: 
					af_mut, af_wt, daf = str(dic_af[va_id][0])[:6], str(dic_af[va_id][1])[:6], str(dic_af[va_id][2])[:6].strip()
					out.write(
						mut_type + "\t" + cont + "\t" + pos + "\t" + str(sp[3]) + "\t" + str(sp[4]) + "\t" + str(sp[5]) + "\t" +  str(sp[6]) + "\t" + str(sp[7]) + "\t" + af_mut + "\t" +
						af_wt + "\t" + daf + "\t" + str(sp[9]) + "\t" + str(sp[10]) + "\t" + str(sp[11]) + "\t" + str(sp[12]) + "\t" + str(sp[13]) + "\t" + str(sp[14]) + "\t" + str(sp[15]) + "\t" + str(sp[16]) + "\t" + str(sp[17]) + "\t" +
						str(sp[18]) + "\t" + str(sp[19]) + "\t" + str(sp[20]) + "\t" + str(sp[21]) + "\t" + str(sp[22]) + "\t" + str(sp[23]) + "\t" + str(sp[24]) + "\t" + str(sp[25]).strip("\n") + "\n"
					)
				except: pass


'''
candidate variants input columns: 
@type   contig  position        ref_base        alt_base        quality ref_count       alt_count       alt_allele_freq   dist_to_selected_pos     hit     mrna_start      mrna_end        strand  gene_model      gene_element    aa_pos  aa_ref aa_alt gene_funct_annot f_primer        tm_f_primer     r_primer        tm_r_primer     upstream        downstream 
	0	1		2				3				4				5		6				7				8					9						10		11				12				13		14				15				16		17		18		19				20				21				22				23				24				25			

output columns: 
@type\tcontig\tposition\tref_base\talt_base\talt_quality\tref_test_count\talt_test_count\talt_allele_freq_test\talt_allele_freq_control\tdAF\tdist_to_qtl_center\thit\tmrna_start\tmrna_end\tstrand\tgene_model\tgene_element\taa_pos\taa_ref\taa_alt\tgene_funct_annot\tf_primer\ttm_f_primer\tr_primer\ttm_r_primer\tupstream_sequence\tdownstream_sequence
    0	1		2			3		4			5			6				7				8						9					10		11					12	13			14			15		16			17			18		19		20		21					22		23			24			25			26					27
'''


# Saving this block for testing
'''
try: 
	af_mut, af_wt, daf = str(dic_af[va_id][0])[:6], str(dic_af[va_id][1])[:6], str(dic_af[va_id][2])[:6]
	out.write(
		mut_type + "\t" + cont + "\t" + pos + "\t" + sp[3] + "\t" + sp[4] + "\t" + sp[5] + "\t" +  sp[6] + "\t" + sp[7] + "\t" + af_mut + "\t" +
		af_wt + "\t" + daf + "\t" + sp[9] + "\t" + sp[10] + "\t" + sp[11] + "\t" + sp[12] + "\t" + sp[13] + "\t" + sp[14] + "\t" + sp[15] + "\t" + sp[16] + "\t" + sp[17] + "\t" +
		sp[18] + "\t" + sp[19] + "\t" + sp[20] + "\t" + sp[21] + "\t" + sp[22] + "\t" + sp[23] + "\t" + sp[24] + "\t" + sp[25] + "\n"
	)
except: pass
'''