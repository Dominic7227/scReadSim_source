import pandas as pd
import numpy as np
import random
from tqdm import tqdm
import sys


def ErrorBase(base, prop, base_call_ref):
	err_base_call_id = np.random.choice(a=[0, 1, 2], size=1, p=prop)[0]      
	err_base_call = base_call_ref[base][err_base_call_id]
	return err_base_call


def ErroneousRead(real_error_rate_read, read_df, output_fq_file):
	n_read = int(np.shape(read_df)[0]/4)
	## Prepare Error rate
	real_error_rate_read_A = real_error_rate_read[['a_to_c_error_rate', 'a_to_g_error_rate', 'a_to_t_error_rate']].to_numpy() 
	real_error_rate_read_A_prop = real_error_rate_read_A/real_error_rate_read_A.sum(axis=1,keepdims=1)
	real_error_rate_read_A_prop[np.isnan(real_error_rate_read_A_prop).any(axis=1),:] = 1/3
	real_error_rate_read_C = real_error_rate_read[['c_to_a_error_rate', 'c_to_g_error_rate', 'c_to_t_error_rate']].to_numpy() 
	real_error_rate_read_C_prop = real_error_rate_read_C/real_error_rate_read_C.sum(axis=1,keepdims=1)
	real_error_rate_read_C_prop[np.isnan(real_error_rate_read_C_prop).any(axis=1),:] = 1/3
	real_error_rate_read_G = real_error_rate_read[['g_to_a_error_rate', 'g_to_c_error_rate', 'g_to_t_error_rate']].to_numpy() 
	real_error_rate_read_G_prop = real_error_rate_read_G/real_error_rate_read_G.sum(axis=1,keepdims=1)
	real_error_rate_read_G_prop[np.isnan(real_error_rate_read_G_prop).any(axis=1),:] = 1/3
	real_error_rate_read_T = real_error_rate_read[['t_to_a_error_rate', 't_to_c_error_rate', 't_to_g_error_rate']].to_numpy() 
	real_error_rate_read_T_prop = real_error_rate_read_T/real_error_rate_read_T.sum(axis=1,keepdims=1)
	real_error_rate_read_T_prop[np.isnan(real_error_rate_read_T_prop).any(axis=1),:] = 1/3
	# Base decision matrix
	real_error_rate_read_prop_dict = {'A': real_error_rate_read_A_prop, 'C': real_error_rate_read_C_prop, 'G': real_error_rate_read_G_prop, 'T': real_error_rate_read_T_prop}
	# Error decision vector
	real_error_rate_read_perbase = real_error_rate_read['error_rate'].to_numpy()
	## Decide whether error occurs for each read
	read_length = real_error_rate_read.shape[0]
	error_read_perbase_indicator = np.zeros((n_read, read_length), dtype=int)
	random.seed(1)
	for base_id in range(read_length):
		error_read_perbase_indicator[:,base_id] = np.random.binomial(n=1, p=real_error_rate_read_perbase[base_id], size=n_read)
	erroneous_read_id = np.where(np.sum(error_read_perbase_indicator, axis=1) > 0)[0]
	## For erroneous reads, generate erroneous base based on the probability matrix
	base_call_ref = {'A': ['C', 'G', 'T'], 'C': ['A', 'G', 'T'], 'G': ['A', 'C', 'T'], 'T': ['A', 'C', 'G']}
	random.seed(2021)
	read_df_witherror = read_df
	# Test. See if conditions reduce the errors in the output reads
	synthetic_error_rate_read_perbase = np.sum(error_read_perbase_indicator, axis=0) / n_read
	base_error_failure_count_vec = np.zeros(90, dtype=int)
	for read_id_tqdm in tqdm(range(len(erroneous_read_id))):
		read_id = erroneous_read_id[read_id_tqdm]
		read_cur = read_df[(read_id*4) : (read_id*4 + 4)]
		bases = list(read_cur[1][0].upper())
		Qscores = list(read_cur[3][0])
		for errorneous_base_id in np.where(error_read_perbase_indicator[read_id,:] > 0)[0]:
			if errorneous_base_id < len(bases):
				base_cur = bases[errorneous_base_id]
				if base_cur in real_error_rate_read_prop_dict:
					prop = real_error_rate_read_prop_dict[base_cur][errorneous_base_id]
					# Decide error base
					err_base_call = ErrorBase(base_cur, prop, base_call_ref)
					bases[errorneous_base_id] = err_base_call
					# Decide Q score
					# Use 9 (Phred score 24) for erroneous base
					Qscores[errorneous_base_id] = '9'
				else:
					base_error_failure_count_vec[errorneous_base_id,] =+ 1
					print("Read ID: " + str(read_id) + " Base Position " + str(errorneous_base_id) + " Base Call " + str(base_cur))
		read_df_witherror[read_id*4+1] = ''.join(bases)
		read_df_witherror[read_id*4+3] = ''.join(Qscores)
	## Write out
	np.savetxt(output_fq_file, read_df_witherror, fmt='%s')

# outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20220407_e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1_NONINPUT"
# real_error_rate_file = outdirectory + "/" + "VerifyQuality/fgbio/Real.error_rate_by_read_position.txt"
# synthetic_fastq_prename = "e18_mouse_brain_fresh_5k_gex_possorted_bam_chr1.syntheticBAM.combined.CBincluded"

def main():
	user_args = sys.argv[1:]
	real_error_rate_file, outdirectory, synthetic_fastq_prename, output_prename = user_args
	# Read in real error rates
	real_error_rate_dir = real_error_rate_file
	real_error_rate = pd.read_csv(real_error_rate_dir, header=0, delimiter="\t")
	# Read in perfect reads
	read2_fq = outdirectory  + "/" + synthetic_fastq_prename + ".read2.bed2fa.fq"
	read2_df = pd.read_csv(read2_fq, header=None).to_numpy()
	# Real data quality score
	# Error rate to Qscore
	# err_rate_qscore_read1 = -10 * np.log10(real_error_rate_read1['error_rate'])
	# fastqc_result = "/home/gayan/Projects/scATAC_Simulator/results/20220408_e18_mouse_brain_fresh_5k_atac_possorted_bam_chr1_1_5000000_NONINPUT/VerifyQuality/fastqc_Real/10X_ATAC_chr1_1_5000000_fastqc/fastqc_data.txt"
	# with open(fastqc_result) as f:
	#     fastqc_result_contents = f.read().splitlines()
	# fastqc_result_contents = np.array(fastqc_result_contents)
	# real_quality_range = np.where(fastqc_result_contents == '>>END_MODULE')[0][0:2]
	# real_quality = list(fastqc_result_contents[(real_quality_range[0]+2) : real_quality_range[1]])
	# real_quality_array = np.array([i.split('\t') for i in real_quality])
	ErroneousRead(real_error_rate, read2_df, outdirectory + "/" + output_prename + ".ErrorIncluded.read2.bed2fa.fq") 

if __name__ == '__main__':
  main()     


