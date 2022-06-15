import pandas as pd
import sys
import glob
import pathlib

name = sys.argv[1]
context = name[-5:]
effect_considered = str(sys.argv[2]) # effect considerted to compute the mutation matrix : synonymous, nonsynonymous or all

effect_simulated = 'syn'

dir = "/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/"
input_folder = dir + effect_simulated + '/' + effect_considered + '/' + name + "/mse_tmp/"

list_of_files =glob.glob(input_folder + '/mse_*.txt')
list_of_files_align = [file for file in list_of_files if 'align' in file]
list_of_files_noalign = [file for file in list_of_files if 'align' not in file]

list_of_mse_align = []
list_of_bins_align = []
for file in list_of_files_align:
    f = open(file, 'r')
    mse = float(f.readlines()[0].strip('\n'))
    list_of_mse_align.append(mse)
    file = file.split('/')[-1]
    list_of_bins_align.append(file[4:-4])

list_of_mse_noalign = []
list_of_bins_noalign = []
for file in list_of_files_noalign:
    f = open(file, 'r')
    mse = float(f.readlines()[0].strip('\n'))
    list_of_mse_noalign.append(mse)
    file = file.split('/')[-1]
    list_of_bins_noalign.append(file[4:-4])

# Select bin type with 
min_bin_align = list_of_bins_align[list_of_mse_align == min(list_of_mse_align)]
min_bin_noalign = list_of_bins_align[list_of_mse_noalign == min(list_of_mse_noalign)]

# Write result
output = dir + effect_simulated + '/' + effect_considered + '/best/'
pathlib.Path(output).mkdir(parents=True, exist_ok=True)
fout = open(dir + effect_simulated + '/' + effect_considered +'/best/chosen_model_align_' + name + '.csv', 'w')
fout.write(min_bin_align)
fout.close()

fout = open(dir + effect_simulated + '/' + effect_considered +'/best/chosen_model_' + name + '.csv', 'w')
fout.write(min_bin_noalign)
fout.close()