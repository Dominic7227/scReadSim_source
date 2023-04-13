import pandas as pd
import pickle
import numpy as np
import random
import time
import matplotlib.pyplot as plt 
import matplotlib

outdirectory = "/home/gayan/Projects/scATAC_Simulator/results/20230304_Minnow/coverage_comparison"

# Read Coverage File
real_coverage = pd.read_csv("%s/real.coverage" % (outdirectory), header=None, delimiter="\t")
scReadSim_coverage = pd.read_csv("%s/scReadSim.coverage" % (outdirectory), header=None, delimiter="\t")
minnow_coverage = pd.read_csv("%s/minnow.coverage" % (outdirectory), header=None, delimiter="\t")
minnow_coverage = minnow_coverage.loc[minnow_coverage[0]=="chr1"]
scReadSim_coverage = scReadSim_coverage.loc[scReadSim_coverage[0]=="chr1"]

# Cut ref genome and plot
chr1_length = 195471971 # https://hgdownload-test.gi.ucsc.edu/goldenPath/mm10/bigZips/mm10.chrom.sizes
window_size = 1000
bin_break = [*range(0, chr1_length, window_size)]
def bin_readcoveragedf(df, bin_break):
    print("Binning genome...")
    tic = time.time()
    df['bin'] = pd.cut(df[1], bin_break) # 54 sec
    toc = time.time()
    tic = time.time()
    print("Count each bin...")
    bin_count = df.groupby('bin')[2].sum() # 5.15 sec
    toc = time.time()
    return bin_count

# real
real_bin_count = bin_readcoveragedf(real_coverage, bin_break)
# scReadSim
scReadSim_bin_count = bin_readcoveragedf(scReadSim_coverage, bin_break)
# minnow
minnow_bin_count = bin_readcoveragedf(minnow_coverage, bin_break)

# Draw line plots
from matplotlib.pyplot import figure

figure(figsize=(13, 5.8), dpi=300)
# Individual
plt.plot(bin_break[1:], real_bin_count, rasterized=True, label="real")
plt.plot(bin_break[1:], scReadSim_bin_count, rasterized=True, label="scReadSim")
plt.plot(bin_break[1:], minnow_bin_count, rasterized=True, label="minnow")
plt.legend(loc="upper left")
plt.savefig(outdirectory + '/' + 'Fig_Minnow.20230317.pdf', format="pdf", bbox_inches="tight")
# Stacked
fig, axs = plt.subplots(3, sharex=True, sharey=True, figsize=(16, 6.5), dpi=300)
axs[0].plot(bin_break[1:], real_bin_count, rasterized=True, color='#F8766D')
axs[1].plot(bin_break[1:], scReadSim_bin_count, rasterized=True, color='#00BFC4')
axs[2].plot(bin_break[1:], minnow_bin_count, rasterized=True, color='#7CAE00')
# plt.tick_params(left= False, labelleft=False, labelbottom = False, bottom = False)
plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[])
plt.xticks([86425000], [])
# # y label
# fig.text(0.08, 0.5, 'read coverage', va='center', rotation='vertical', fontsize=25)
# # x label
# fig.text(0.5, 0.04, "reference genome (chr1)", ha='center')
plt.savefig(outdirectory + '/' + 'Fig_Minnow.stacked.20230317.withhightwindow.pdf', format="pdf", bbox_inches="tight")

# Calculate pearson correlation with real
import scipy.stats
pearson_scReadSim = scipy.stats.pearsonr(real_bin_count, scReadSim_bin_count)  
spearmanr_scReadSim = scipy.stats.spearmanr(real_bin_count, scReadSim_bin_count)
kendalltau_scReadSim = scipy.stats.kendalltau(real_bin_count, scReadSim_bin_count)  

pearson_minnow = scipy.stats.pearsonr(real_bin_count, minnow_bin_count)  
spearmanr_minnow = scipy.stats.spearmanr(real_bin_count, minnow_bin_count) 
kendalltau_minnow = scipy.stats.kendalltau(real_bin_count, minnow_bin_count) 
