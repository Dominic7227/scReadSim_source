import pandas as pd
import numpy as np
import sys

def GeneratePseudoTruePeak(outdirectory, input_transcript_peakfile, output_true_peakfile):
	mm10_transcript_coordinates = pd.read_csv("%s/%s" % (outdirectory, input_transcript_peakfile), delimiter="\t",  names=['chr', 'start', 'end'])
	peak_len = 400
	jitter_value_vec = np.random.random_integers(-150,150,size=np.shape(mm10_transcript_coordinates)[0])  # nrow(reads_cur) should equal to nfrag_cur
	output_df = mm10_transcript_coordinates[['chr', 'start']]
	output_df['end'] = mm10_transcript_coordinates['start'] + peak_len + jitter_value_vec
	output_df.to_csv("%s/%s" % (outdirectory, output_true_peakfile), header=None, index=None, sep='\t')

def main():
	user_args = sys.argv[1:]
	outdirectory, input_transcript_peakfile, output_true_peakfile = user_args
	GeneratePseudoTruePeak(outdirectory, input_transcript_peakfile, output_true_peakfile)

if __name__ == '__main__':
  main()     