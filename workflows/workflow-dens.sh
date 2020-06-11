#!/bin/bash
#
# This is the command sent by 'master.sh':
# ./workflow-x.sh $my_log_file $project_name $workflow $data_source $lib_type $ins_seq $read_s $read_f $read_r $gff_file $ann_file
#
# This command is always the same, regardless of the workflow
#
# Fields that can be equal to 'n/p' (= not provided) are the following: $ins_seq, $read_s, $read_f, $read_r, $ann_file
# If $ins_seq = n/p, that is because $workflow = snp. Therefore, no insertion seq is needed and $ins_seq is ignored by hte program.
# If $read_s = n/p, that is because $lib_type = pe, so it is ignored by the program.
# If $read_f and $read_r = n/p, that is because $lib_type = se, so it is ignored by the program.
# If $ann_file = n/p, this is because user did no have it. This program has the deal with this: if data not provided, simply do not
# include gene annotated info to the report.
# 
# my_log_file		>	$1
# project_name		>	$2
# workflow			>	$3
# data_source		>	$4
# lib_type			>	$5
# ins_seq			>	$6
# read_s			>	$7
# read_f			>	$8
# read_r			>	$9
# gff_file			>	${10}
# ann_file			>	${11}
# read_s_par							>	${12}
# read_f_par							>	${13}
# read_r_par							>	${14}
# cross_type							>	${15} 
# is_ref_strain							>	${16} 
# parental_reads_provided				>	${17}
# snp_analysis_type [par/f2wt]			>	${18}
# lib_type_control						>	${19}
# stringency							>	${20}


# Some initial parameters
start_time=`date +%s`	
exit_code=0 				# Set 'exit_code' (flag variable) to 0
my_log_file=$1 				# Set location of log file
export location="$PWD" 			#Save path to hisat2-build and hisat2 in variable BT2

# Create input variables
my_log_file=$1
project_name=$2
my_sample_mode=$5 												#[pe, se], paired/single  
my_control_mode=${19}											#
my_rd=$7											 			#reads (single)
my_rf=$8 														#forward reads
my_rr=$9												 		#reverse reads 			
my_p_rd=${12}											 		#reads (single) control	
my_p_rf=${13} 													#forward reads control	
my_p_rr=${14}											 		#reverse reads control	
my_gs=gnm_ref_merged/genome.fa 									#genome sequence
my_ix=genome_index 							
my_gff=${10}													#Genome feature file
my_ann=${11}													
my_rrl=250 														#Regulatory region length
my_log_file=$1
my_mut=snp  													#my_mut takes the values 'snp' in this workflow and 'lin' in the large insertions workflow, for the execution of the graphic output module
snp_analysis_type=${18}
stringency=${20}
exp_mut_type=${21}

#Set number of maximum CPU for steps compatible with multithreading, default = 1 
threads=2

# Set internal variables according to the SNP validation stringency chosen by the user
if [ $stringency == high_stringency ]; then
	problemSample_bowtie_mp="--mp 6,2"
	problemSample_mpileup_C="-C50"
	problemSample_snpQualityTheshold="120"
else
	problemSample_bowtie_mp="--mp 3,2"
	problemSample_mpileup_C=""
	problemSample_snpQualityTheshold="30"
fi

# Define the folders in the easymap directory 
f0=user_data
f1=$project_name/1_intermediate_files
f2=$project_name/2_logs
f3=$project_name/3_workflow_output

# Write PID to status file
my_status_file=$f2/status
echo 'pid workflow '$$ >> $my_status_file

# Check genome size to set interval_width
{
	interval_width=`python2 $location/scripts_snp/set-interval.py -a $f1/$my_gs`
} || {
	interval_width=4000007
	echo $(date "+%F > %T")': set-interval.py failed.' >> $my_log_file
}
echo $(date "+%F > %T")': set-interval.py finished, interval set at: '$interval_width   >> $my_log_file

# Set variant density mapping window width and step according to genome size
if [ $interval_width -lt 4000009 ]; then
	dens_width=500000
	dens_step=250000
elif [ $interval_width -gt 4000009 ] && [ $interval_width -lt 11000000 ]; then
	dens_width=1000000
	dens_step=1000000
elif [ $interval_width -gt 4000009 ]; then
	dens_width=3000000
	dens_step=3000000
fi

##################################################################################################################################################################################
#																																												 #
#																																												 #
#																	F2 FQ PROCESSING FUNCTION																					 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

#Run hisat2-build on genome sequence 
{
	$location/hisat2/hisat2-build $f1/$my_gs $f1/$my_ix 1> $f2/hisat2-build_std1.txt 2> $f2/hisat2-build_std2.txt

} || {
	echo $(date "+%F > %T")': HISAT2-build on genome sequence returned an error. See log files.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': HISAT2-build finished.' >> $my_log_file

function get_problem_va {  
	if [ $my_sample_mode == se ] 
	then
		#Run hisat2 unpaired to align raw reads to genome 
		{
			$location/hisat2/hisat2 -p $threads -x $f1/$my_ix -U $my_rd -S $f1/alignment1.sam 2> $f2/hisat2_problem-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': HISAT2 returned an error during the aligment of F2 reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': HISAT2 finished the alignment of F2 reads to genome.' >> $my_log_file
	fi

	if [ $my_sample_mode == pe ] 
	then
		#Run hisat2 paired to align raw reads to genome 
		{
			$location/hisat2/hisat2  -p $threads  -x $f1/$my_ix -1 $my_rf -2 $my_rr -S $f1/alignment1.sam 2> $f2/hisat2_problem-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': HISAT2 returned an error during the aligment of F2 reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': HISAT2 finished the alignment of F2 reads to genome.' >> $my_log_file
	fi

	#SAM to BAM
	{
		$location/samtools1/samtools sort  -@ $threads  $f1/alignment1.sam > $f1/alignment1.bam 2> $f2/sam-to-bam_problem-sample_std2.txt
		rm -rf ./user_projects/$project_name/1_intermediate_files/alignment1.sam

	} || {
		echo 'Error transforming SAM to BAM.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': SAM to BAM finished.' >> $my_log_file

	#Variant calling
	{

		$location/samtools1/samtools mpileup  -t DP,ADF,ADR $problemSample_mpileup_C -uf $f1/$my_gs $f1/alignment1.bam 2> $f2/mpileup_problem-sample_std.txt | $location/bcftools-1.3.1/bcftools call -mv -Ov > $f1/raw_variants.vcf 2> $f2/call_problem-sample_std.txt
		# -B: Disables probabilistic realignment for the computation of base alignment quality (BAQ). Applying this argument reduces the number of false negatives during the variant calling
		# -t DP,ADF,ADR: output VCF file contains the specified optional columns: read depth (DP), allelic depths on the forward strand (ADF), allelic depths on the reverse strand (ADR)
		# -uf: uncompressed vcf output / fasta imput genome file
		# -mv: include only polymorphic sites in output
		# -Ov: uncompressed VCF output file 
		# -C50: reduce the effect of reads with excessive mismatches. This aims to fix overestimated mapping quality

	} || {
		echo $(date "+%F > %T")': Error during variant-calling of F2 data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': F2 data variant calling finished.' >> $my_log_file

	#Groom vcf
	{
		python2 $location/scripts_snp/vcf-groomer.py -a $f1/raw_variants.vcf -b $f1/F2_raw.va  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of vcf-groomer.py with F2 data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF grooming of F2 data finished.' >> $my_log_file

	#RD graphics
	depth_alignment $f1/alignment1.bam $f3/frequence_depth_alignment_distribution_sample.png

	#Run vcf filter
	if [ $my_cross == bc ]; then mut_type=all ; fi		#SDL V2 Update, prev: then mut_type=EMS
	if [ $my_cross == oc ]; then mut_type=all ; fi

	if [ $av_rd -gt 25 ]; then dp_min=15 ; fi
	if [ $av_rd -le 25 ]; then dp_min=10 ; fi

	dp_max=$(($av_rd * 3))
	if [ $dp_max -le 40 ]; then dp_max=100 ; fi


	{
		python2 $location/scripts_snp/variants-filter.py -a $f1/F2_raw.va -b $f1/F2_filtered.va -step 3 -fasta $f1/$my_gs -dp_min $dp_min -dp_max $dp_max -qual_min $problemSample_snpQualityTheshold -mut_type all  2>> $my_log_file

	} || {
		echo 'Error during execution of variants-filter.py with F2 data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': First VCF filtering step of F2 data finished.' >> $my_log_file

	#Intermediate files cleanup
	rm -f $f1/*.sam
	rm -f $f1/*.vcf

}


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																	control FQ PROCESSING FUNCTION																				 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

function get_control_va { 
	if [ $my_control_mode == se ] 
	then
		#Run hisat2 unpaired to align raw reads to genome 
		{
			$location/hisat2/hisat2  -p $threads  -x $f1/$my_ix -U $my_p_rd -S $f1/alignment1P.sam 2> $f2/hisat2_control-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': HISAT2 returned an error during the aligment of control reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': HISAT2 finished the alignment of control reads to genome.' >> $my_log_file
	fi

	if [ $my_control_mode == pe ] 
	then
		#Run hisat2 paired to align raw reads to genome 
		{
			$location/hisat2/hisat2  -p $threads -x $f1/$my_ix -1 $my_p_rf -2 $my_p_rr -S $f1/alignment1P.sam 2> $f2/hisat2_control-sample_std2.txt

		} || {
			echo $(date "+%F > %T")': HISAT2 returned an error during the aligment of control reads. See log files.' >> $my_log_file
			exit_code=1
			echo $exit_code
			exit
		}
		echo $(date "+%F > %T")': HISAT2 finished the alignment of control reads to genome.' >> $my_log_file
	fi

	#SAM to BAM
	{
		$location/samtools1/samtools sort  -@ $threads $f1/alignment1P.sam > $f1/alignment1P.bam 2> $f2/sam-to-bam_control-sample_std2.txt

		rm -rf ./user_projects/$project_name/1_intermediate_files/alignment1P.sam

	} || {
		echo $(date "+%F > %T")': Error transforming SAM to BAM' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': SAM to BAM finished' >> $my_log_file

	#Variant calling
	{

		$location/samtools1/samtools mpileup  -t DP,ADF,ADR  -uf $f1/$my_gs $f1/alignment1P.bam 2> $f2/mpileup_control-sample_std.txt | $location/bcftools-1.3.1/bcftools call -mv -Ov > $f1/raw_p_variants.vcf 2> $f2/call_control-sample_std.txt

	} || {
		echo $(date "+%F > %T")': Error during variant-calling of control data' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': Control data variant calling finished' >> $my_log_file

	#Groom vcf
	{
		python2 $location/scripts_snp/vcf-groomer.py -a $f1/raw_p_variants.vcf -b $f1/control_raw.va  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of vcf-groomer.py with control data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': VCF grooming of control data finished.' >> $my_log_file

	#RD graphics
	depth_alignment $f1/alignment1P.bam $f3/frequence_depth_alignment_distribution_control.png

	#Run vcf filter

	if [ $av_rd -gt 25 ]; then dp_min=15 ; fi
	if [ $av_rd -le 25 ]; then dp_min=10 ; fi

	dp_max=$(($av_rd * 3))
	if [ $dp_max -le 40 ]; then dp_max=100 ; fi

	{
		python2 $location/scripts_snp/variants-filter.py -a $f1/control_raw.va -b $f1/control_filtered.va -step 3 -fasta $f1/$my_gs -dp_min 10 -dp_max $dp_max -qual_min 20  2>> $my_log_file

	} || {
		echo $(date "+%F > %T")': Error during execution of variants-filter.py with control data.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
	echo $(date "+%F > %T")': First VCF filtering step of control data finished.' >> $my_log_file

	#Intermediate files cleanup
	rm -f $f1/*.sam
	rm -f $f1/*.vcf
}


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																			DEPTH ALIGNMENT ANALYSIS FUNCTION																	 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

function depth_alignment {
	{
		python2 $location/scripts_snp/depth_measures_generation.py -genome $f1/$my_gs -bam $1 -out $f1/coverage_alignment1.txt  2>> $my_log_file
		rm -rf ./user_projects/$project_name/1_intermediate_files/alignment1.bam

	} || {
		echo $(date "+%F > %T")': Error during obtaining of alignment depth .' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}

	{
		av_rd=`python2 $location/graphic_output/graphic-alignment.py -coverages $f1/coverage_alignment1.txt   -out $2  2>> $my_log_file `

	} || {
		echo $(date "+%F > %T")': Error during Graphic_alignment execution in sample alignment.' >> $my_log_file
		exit_code=1
		echo $exit_code
		exit
	}
}


##################################################################################################################################################################################
#																																												 #
#																																												 #
#																			DATA ANALYSIS 																						 #
#																																												 #
#																																												 #
##################################################################################################################################################################################

#_________________________________________________________________________________________________________________________________________________________________________________

#_________________________________________________________________________________________________________________________________________________________________________________


# (1) Get problem and control VA files
get_problem_va 
get_control_va


# (2) Run VA operations: Remove control SNPs from problem file
my_operation_mode=A
{
	python2 $location/scripts_snp/variants-operations.py -a $f1/F2_filtered.va -b $f1/control_filtered.va -c $f1/F2_control_comparison.va -mode $my_operation_mode -primary 1  2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during first execution of variants-operations.py .' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': VCF operations finished.' >> $my_log_file

# (3) Generate files for variant density analysis
{
	python2 $location/scripts_snp/variants-filter.py -a $f1/F2_control_comparison.va -b $f1/F2_filtered_EMS.va -step 3 -fasta $f1/$my_gs -dp_min $dp_min -dp_max $dp_max -qual_min $problemSample_snpQualityTheshold -mut_type EMS  2>> $my_log_file
	python2 $location/scripts_snp/variants-filter.py -a $f1/F2_filtered_EMS.va -b $f1/F2_filtered_EMS_hz.va -step 3 -fasta $f1/$my_gs -af_min 0.98 -dp_min $dp_min -dp_max $dp_max -qual_min $problemSample_snpQualityTheshold -mut_type EMS  2>> $my_log_file
	python2 $location/scripts_snp/variants-filter.py -a $f1/F2_control_comparison.va -b $f1/F2_hz.va -step 3 -fasta $f1/$my_gs -af_min 0.98 -dp_min $dp_min -dp_max $dp_max -qual_min $problemSample_snpQualityTheshold   2>> $my_log_file

} || {
	echo 'Error during execution of variants-filter.py with F2 data.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': First VCF filtering step of F2 data finished.' >> $my_log_file

# (4) Run variant density analysis on the files
{
	python2 $location/scripts_snp/dens-counter.py -in $f1/F2_filtered.va  -out $f1/F2_filtered_dens.va -width $dens_width -step $dens_step -f_input $f1/$my_gs 2>> $my_log_file
	python2 $location/scripts_snp/dens-counter.py -in $f1/F2_control_comparison.va  -out $f1/F2_control_comparison_dens.va -width $dens_width -step $dens_step -f_input $f1/$my_gs 2>> $my_log_file
	python2 $location/scripts_snp/dens-counter.py -in $f1/F2_filtered_EMS.va  -out $f1/F2_filtered_EMS_dens.va -width $dens_width -step $dens_step -f_input $f1/$my_gs 2>> $my_log_file
	python2 $location/scripts_snp/dens-counter.py -in $f1/F2_filtered_EMS_hz.va  -out $f1/F2_filtered_EMS_hz_dens.va -width $dens_width -step $dens_step -f_input $f1/$my_gs 2>> $my_log_file
	python2 $location/scripts_snp/dens-counter.py -in $f1/F2_hz.va  -out $f1/F2_hz_dens.va -width $dens_width -step $dens_step -f_input $f1/$my_gs 2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during execution of dens-counter.py .' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': Mutation density counter module finished.' >> $my_log_file

# (5) Run variant density mapping 
{
	# Step 1
	python2 $location/scripts_snp/dens-mapping.py -in $f1/F2_filtered_dens.va  -out $f1/mapping_max.va -step 1							2>> $my_log_file
	python2 $location/scripts_snp/dens-mapping.py -in $f1/F2_control_comparison_dens.va  -out $f1/mapping_max.va -step 1				2>> $my_log_file
	python2 $location/scripts_snp/dens-mapping.py -in $f1/F2_filtered_EMS_dens.va  -out $f1/mapping_max.va -step 1						2>> $my_log_file
	python2 $location/scripts_snp/dens-mapping.py -in $f1/F2_filtered_EMS_hz_dens.va  -out $f1/mapping_max.va -step 1 					2>> $my_log_file
	python2 $location/scripts_snp/dens-mapping.py -in $f1/F2_hz_dens.va  -out $f1/mapping_max.va -step 1  								2>> $my_log_file
	# Step 2
	python2 $location/scripts_snp/dens-mapping.py -in $f1/mapping_max.va  -out $f1/candidate_region.txt -step 2 -mut_type $exp_mut_type -pname $project_name	2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during execution of dens-mapping.py .' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': Mutation density mapping module finished.' >> $my_log_file

# (6) Candidate region filtering
{
	python2 $location/scripts_snp/variants-filter.py -a $f1/F2_control_comparison.va -b $f1/final_variants.va -step 2 -cand_reg_file $f1/candidate_region.txt -mut_type $exp_mut_type  2>> $my_log_file
} || {
	echo $(date "+%F > %T")': Error during determinarion of candidate region .' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': Candidate region filter finished.' >> $my_log_file

# Create input for varanalyzer and run varanalyzer.py (one file for the candidate region and one for the whole genome)
{
	python2 $location/scripts_snp/snp-to-varanalyzer.py -a $f1/final_variants.va -b $f1/snp-to-varanalyzer.txt  2>> $my_log_file
	python2 $location/scripts_snp/snp-to-varanalyzer.py -a $f1/F2_control_comparison.va -b $f1/snp-to-varanalyzer-total.txt  2>> $my_log_file		

} || {
	echo $(date "+%F > %T")': Error during execution of snp-to-varanalyzer.py .' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': Input for varanalyzer finished.' >> $my_log_file

# Varanalyzer
{
	python2 $location/varanalyzer/varanalyzer.py -itp snp -con $f1/$my_gs -gff $f0/$my_gff -var $f1/snp-to-varanalyzer.txt -rrl $my_rrl -pname $project_name -ann $f0/$my_ann -out $f1/varanalyzer_output_snp.txt 2>> $my_log_file
	python2 $location/varanalyzer/varanalyzer.py -itp snp -con $f1/$my_gs -gff $f0/$my_gff -var $f1/snp-to-varanalyzer-total.txt -rrl $my_rrl -pname $project_name -ann $f0/$my_ann -out $f1/varanalyzer_output_total_temp.txt  2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during execution of varanalyzer.py .' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': Varanalyzer finished.' >> $my_log_file

# (7) Generate graphic output
{
	python2 $location/graphic_output/graphic-output.py -my_mut dens -pname $project_name  -bsnp $f1/$my_gs  -cr_file $f1/candidate_region.txt -1 $f1/F2_filtered_dens.va -2 $f1/F2_control_comparison_dens.va -3 $f1/F2_filtered_EMS_dens.va -4 $f1/F2_filtered_EMS_hz_dens.va -5 $f1/F2_hz_dens.va -mut_type $exp_mut_type 2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during generation of graphic output.' >> $my_log_file
	exit_code=1
	echo $exit_code
	exit
}
echo $(date "+%F > %T")': Graphic output generated.' >> $my_log_file

# (7) Generate report
{ 
	cp $f1/varanalyzer_output_snp.txt $f3/candidate_variants.txt
	cp $f1/varanalyzer_output_total_temp.txt $f3/candidate_variants_total.txt
	python2 $location/graphic_output/report.py  -log $my_log_file -output_html $f3/report.html -project $project_name -mut_type dens -files_dir $f3 -cand_reg_file $f1/candidate_region.txt  2>> $my_log_file

} || {
	echo $(date "+%F > %T")': Error during generation of report file.' >> $my_log_file
	#exit_code=1									# REMOVED FOR TESTING
	#echo $exit_code
	#exit
}
echo $(date "+%F > %T")': Report file generated.' >> $my_log_file


#Intermediate files cleanup
rm -f $f1/*.bam
rm -f $f1/*.bai

echo $exit_code
