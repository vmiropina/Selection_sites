import pathlib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys


input_dir = "/nfs/users/dweghorn/projects/Selection_sites/"
dir = "/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/"

def do_plots(name, cov, effect_considered, effect_simulated, bin_type, suffix, reps):
    if 'align' in bin_type: 
        realdatafile = input_dir + 'data/data_all_covariates_' + effect_simulated + '_align_sites/data_all_covariates_' + effect_simulated + '_align_sites_%s.bed.gz'%name
    else:
        realdatafile = input_dir + 'data/data_all_covariates_' + effect_simulated + '_sites/data_all_covariates_' + effect_simulated + '_sites_%s.bed.gz'%name
    columns = ['obs']
    df =  pd.read_csv(realdatafile, sep = '\t', usecols = columns, compression='gzip')
    m = np.min(df['obs'])
    M = np.max(df['obs'])
    counts , bins = np.histogram(df['obs'], bins = np.linspace(m,M+2,int(M-m+3)))
    del df


    # Read simulated data
    simulatedfolder =  dir + effect_simulated + '/' + effect_considered + '/' + name + "/simulated_data/"
    simulatedfile = simulatedfolder + bin_type + suffix + "_" + cov + ".csv"

    df =  pd.read_csv(simulatedfile, sep = '\t')

    counts_0, bins_0 = np.histogram(df['0'] , bins = np.linspace(m,M+2,int(M-m+3)))
    for i in range(1,reps):
        counts_0 += np.histogram(df[str(i)] , bins = np.linspace(m,M+2,int(M-m+3)))[0]
    del df

    # Read real data
    log_counts = np.nan_to_num((np.log(counts+1)), posinf=0.0, neginf=0.0, nan=0.0)
    log_counts_0 = np.nan_to_num((np.log(counts_0/reps+1)), posinf=0.0, neginf=0.0, nan=0.0)

    output = dir + effect_simulated + '/' + effect_considered + '/' + name + '/plots/'
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    plt.bar(np.array(bins[0:-1])-0.1, log_counts, width=0.2, label='Observed')
    plt.bar(np.array(bins_0[0:-1])+0.1, log_counts_0, width=0.2, label='Simulated')
    plt.savefig(output + '/hist_' + bin_type + suffix + "_" + cov + '.png')
    plt.legend()
    plt.close()
    return log_counts, log_counts_0

def compute_error(log_counts, log_counts_0, name, bin_type, effect_considered, cov, suffix):
    mse = sum((log_counts_0-log_counts)**2)/len(log_counts)
    output = dir + effect_simulated + '/' + effect_considered + '/' + name + '/mse_tmp/'
    pathlib.Path(output).mkdir(parents=True, exist_ok=True)
    fout = open(output + '/mse_' + bin_type + suffix + "_" + cov + '.txt', 'w')
    fout.write(str(mse))
    fout.close()

name = sys.argv[1]
context = name[-5:]
bin_type =  str(sys.argv[2]) 
effect_considered = str(sys.argv[3]) # effect considerted to compute the mutation matrix : synonymous, nonsynonymous or all

effect_simulated = 'syn'
reps = 20

if 'genes' in bin_type:
    cov = 'cbase'
    if 'cpgi' in bin_type:
        log_counts, log_counts_0 = do_plots(name, cov, effect_considered, effect_simulated, bin_type, '', reps)
        compute_error(log_counts, log_counts_0, name, bin_type, effect_considered, cov, '')
        log_counts, log_counts_0 = do_plots(name, cov, effect_considered, effect_simulated, bin_type, '_01', reps)
        compute_error(log_counts, log_counts_0, name, bin_type, effect_considered, cov, '_01')
    else:
        log_counts, log_counts_0 = do_plots(name, cov, effect_considered, effect_simulated, bin_type, '', reps)
        compute_error(log_counts, log_counts_0, name, bin_type, effect_considered, cov, '')

cov = 'cov'
if 'cpgi' in bin_type:
    log_counts, log_counts_0 = do_plots(name, cov, effect_considered, effect_simulated, bin_type, '', reps)
    compute_error(log_counts, log_counts_0, name, bin_type, effect_considered, cov, '')
    log_counts, log_counts_0 = do_plots(name, cov, effect_considered, effect_simulated, bin_type, '_01', reps)
    compute_error(log_counts, log_counts_0, name, bin_type, effect_considered, cov, '_01')
else:
    log_counts, log_counts_0 = do_plots(name, cov, effect_considered, effect_simulated, bin_type, '', reps)
    compute_error(log_counts, log_counts_0, name, bin_type, effect_considered, cov, '')
