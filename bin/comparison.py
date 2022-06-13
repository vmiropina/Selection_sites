import pandas as pd
import sys
import glob

name = sys.argv[1]
context = name[-5:]
effect_considered = str(sys.argv[2]) # effect considerted to compute the mutation matrix : synonymous, nonsynonymous or all

effect_simulated = 'syn'

dir = "/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/"
input_folder = dir + effect_simulated + '/' + effect_considered + '/' + name + "/mse_tmp/"

list_of_files =glob.glob(input_folder + '/mse_*.csv')
list_of_mse = []
list_of_bins = []
for file in list_of_files:
    f = open(file, 'r')
    mse = float(f.readlines()[0].strip('\n'))
    list_of_mse.append(mse)
    list_of_bins.append(file[4:-4])

# Select bin type with 
min_bin = list_of_bins[list_of_mse == min(list_of_mse)]

# Write result
fout = open(dir + effect_simulated + '/' + effect_considered +'/best/chosen_model_' + name + '.csv', 'w')
fout.write(min_bin)
fout.close()