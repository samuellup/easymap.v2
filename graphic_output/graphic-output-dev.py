# This script calls on functions from draw.py to draw the images required for each annalysis


# python2 graphic-output-dev.py -my_mut qtl -asnp F2_control_comparison.va -gff SL3.gff  -iva varanalyzer_output.txt -rrl 150 -pname project  -bsnp genome.fa  -cr_file candidate_region.txt  -mut_type EMS 
# python2 graphic-output-dev.py -my_mut qtl -asnp F2_control_comparison.va -gff SL3.gff  -iva varanalyzer_output.txt -rrl 150 -pname project  -bsnp genome.fa  -cr_file candidate_region.txt  -mut_type EMS


#from __future__ import division
import argparse
from drawdev import fa_vs_pos, insertions_overview_and_histograms, gene_plot, legend, candidates_zoom, gene_legend, dens_graphs, dens_ovw, qtl_plot, qtl_plot_points
parser = argparse.ArgumentParser()

#INPUT VARIABLES FOR SNP
parser.add_argument('-my_mut', action="store", dest = 'my_mut')			#snp or lin
parser.add_argument('-asnp', action="store", dest = 'input_snp')		
parser.add_argument('-bsnp', action="store", dest = 'input_f_snp')		#Fasta genome input
parser.add_argument('-snp_analysis_type', action="store", dest='my_snp_analysis_type')
parser.add_argument('-cross', action="store", dest='my_cross')
parser.add_argument('-interval_width', action="store", dest='interval_width')

#INPUT VARIABLES FOR LIN
parser.add_argument('-a', action="store", dest = 'input')		
parser.add_argument('-b', action="store", dest = 'input_f')				#Fasta genome input
parser.add_argument('-gff', action="store", dest = 'gff')				#Genome feature file
parser.add_argument('-m', action="store", dest = 'mode', default = 'pe')
parser.add_argument('-ins_pos', action="store", dest = 'ins_pos')

# INPUT VARIABLES FOR DENS
parser.add_argument('-1', action="store", dest = 'F2')		
parser.add_argument('-2', action="store", dest = 'F2_filtered')				
parser.add_argument('-3', action="store", dest = 'EMS')
parser.add_argument('-4', action="store", dest = 'EMS_hz')	
parser.add_argument('-5', action="store", dest = 'HZ')		
parser.add_argument('-mut_type', action="store", dest = 'mut_type', default="EMS")		
parser.add_argument('-cr_file', action="store", dest = 'cr_file')		

#SHARED VARIABLES
parser.add_argument('-iva', action="store", dest = 'input_va')	 		#Output de varanalyzer
parser.add_argument('-rrl', action="store", dest = 'rrl') 				#Regulatory region lenght
parser.add_argument('-pname', action="store", dest='project_name')

args = parser.parse_args()
project = args.project_name


if args.my_mut == 'qtl':
	qtl_plot()
	qtl_plot_points()

	#gene_plot()	
