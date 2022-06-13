import numpy as np
import pandas as pd
import sys
import pathlib

from generate_muts import *


#########################################################################################
#########################################################################################
#Define functions 

def process_mutation_matrix(mut_matrix, context):
    if context == "3mers":
        number_muts = 4**3*3
        values = [str(i) for i in range(4**3) for k in range(3)]
        muts = ['1','2','3']*16 + ['0','2','3']*16 + ['0','1','3']*16 + ['0','1','2']*16 
        values = [values[i]+muts[i] for i in range(len(muts))]
    elif context == "5mers":
        number_muts = 4**5*3
        values = [str(i) for i in range(4**5) for k in range(3)]
        muts = ['1','2','3']*16 + ['0','2','3']*16 + ['0','1','3']*16 + ['0','1','2']*16 
        muts = muts*16
        values = [values[i]+muts[i] for i in range(len(muts))]
    elif context == "7mers":
        number_muts = 4**7*3
        values = [str(i) for i in range(4**7) for k in range(3)]
        muts = ['1','2','3']*4**3 + ['0','2','3']*4**3 + ['0','1','3']*4**3 + ['0','1','2']*4**3 
        muts = muts*4**3
        values = [values[i]+muts[i] for i in range(len(muts))]
    elif context == "9mers":
        number_muts = 4**9*3
        values = [str(i) for i in range(4**9) for k in range(3)]
        muts = ['1','2','3']*4**4 + ['0','2','3']*4**4 + ['0','1','3']*4**4 + ['0','1','2']*4**4 
        muts = muts*4**4
        values = [values[i]+muts[i] for i in range(len(muts))]
    mut_mat_vec = np.zeros(number_muts)
    for i in range(number_muts):
        x = int(values[i][0:-1])
        y = int(values[i][-1])
        mut_mat_vec[i] = mut_matrix.loc[x][y]
    return mut_mat_vec, number_muts, values


def loop_exons_or_genes(P, mistarget, mut_mat_vec, outfile, bin, flag, effect, theta, m, number_muts, values, reps):
    if bin=='gene':
        COLUMN_NAMES = ['gene', 'context', 'pos', 'mi', 'l_s' ,'lambda_s' ]
    elif bin=='exon':
        COLUMN_NAMES = ['exon', 'context', 'pos', 'mi', 'l_s' ,'lambda_s' ]
    #for i in range(len(P)):
    for i in range(10):
        genes = []
        contexts = []
        pos = [] 
        nobs = [] 
        mi = [] 
        ls = []
        lambdas = []
        gene = P.iloc[i][bin]
        s_obs = P.iloc[i]['obs']
        lambda_s = P.iloc[i]['prediction'] 
        mis_counts = np.array(mistarget.loc[gene])[0:-1]
        l_s = sum(mis_counts*mut_mat_vec)
        for context in range(number_muts):
            m_i = mut_mat_vec[context]
            c = int(mis_counts[context])
            try :  
                if effect == "nonsyn":
                    r = m_i/l_s     
                    n_obs = generate_n(r, s_obs, theta, m, c, 200)
                elif effect == "syn":
                    print(str(lambda_s*m_i/l_s))
                    n_obs = generate_s(lambda_s, m_i, l_s, (c, reps))
                    print(np.max(n_obs))
            except :
                n_obs = np.zeros(c)
            for j in range(len(n_obs)): 
                nobs.append(n_obs[j])
                genes.append(gene)
                contexts.append(values[context])
                mi.append(m_i)
                lambdas.append(lambda_s)
                ls.append(l_s)
                pos.append(j)
        df = pd.DataFrame(list(zip(genes, contexts, pos,  mi, ls, lambdas)), columns =COLUMN_NAMES)
        mat = pd.DataFrame(np.array(nobs))
        df2 = pd.concat([df, mat], axis = 1)
        if i == 0 and flag == 1:
            df2.to_csv(outfile , sep ="\t", index = False, header=True)
        else:
            df2.to_csv(outfile , sep ="\t", index = False, mode='a', header=False)
    print("All genes done!")   

def simulate_selection_sites(name, context, suffix, bin_type, cov, output_folder, reps, list_of_genes, effect):
    out_folder = dir + effect_simulated + '/' + effect_considered + '/' + name + "/simulated_data/"
    if suffix =='_01':
        out_file = out_folder + bin_type + suffix + "_" + cov + ".csv"
    else: 
        out_file = out_folder + bin_type + "_" + cov + ".csv"
        
    pathlib.Path(out_folder).mkdir(parents=True, exist_ok=True)
    if cov =='cbase' and suffix != '_01':
        if 'align' in bin_type:
            P = pd.read_csv(input_dir + 'New_CBaSE_covariates/' + output_folder + '/output_data_preparation_' + name + '_align' +  suffix + '.txt' , sep="\t", header=None)
        else:
            P = pd.read_csv(input_dir + 'New_CBaSE_covariates/' + output_folder + '/output_data_preparation_' + name +   suffix + '.txt' , sep="\t", header=None)
        P.columns = ['gene', 'lm', 'lk', 'ls', 'mobs', 'kobs', 'obs', 'len', 'fcov', 'prediction']
    else: 
        P = pd.read_csv(dir  + effect_simulated + '/' + effect_considered + '/'  + name + '/glm' + '/fcovariates/fcovariates_%s%s.csv.gz'%(bin_type, suffix), sep=',', compression='gzip')

    P = P[P['gene'].isin(list_of_genes)]
    P = P.reset_index(drop =True)

    if effect == 'nonsyn':
        #THIS HAS TO BE CHECKED
        with open(input_dir + 'New_CBaSE_covariates/' + output_folder + '/used_params_and_model_' + name + '.txt' ) as file:
            for line in file:
                line = line.split(',')
                m = float(line[-1])
                theta = [float(i) for i in line[0:-1]]
    else:
        theta = None
        m = None

    if 'cpgi' in bin_type:
        bin_type_target = bin_type[:-5]
    else:
        bin_type_target = bin_type

    if 'sites' in bin_type: 
        for i in range(reps):
            P[i] = np.random.poisson(P['prediction'])
        P.to_csv(out_file , sep ="\t", index = False, header=True)
    else: 
        if suffix != '_01': 
            mis_target =  pd.read_csv(input_dir + 'Target_sizes/synonymous_' + bin_type_target +'_'+ context  + suffix + '.csv.gz', sep="\t", header=None, index_col = 0, compression='gzip')
            mut_matrix = pd.read_csv(input_dir + 'New_CBaSE_covariates/' + output_folder + '/mutation_matrix_' + name + suffix + '.txt', sep = '\t', header = None).drop(columns = [0]).rename(columns = {1:0, 2:1, 3:2, 4:3})
            mut_mat_vec, number_muts, values = process_mutation_matrix(mut_matrix, context)
            if suffix != '_1':
                if 'genes' in bin_type:
                    loop_exons_or_genes(P, mis_target, mut_mat_vec, out_file, 'gene', 1, effect, theta, m, number_muts, values, reps)
                elif 'exons' in bin_type :
                    loop_exons_or_genes(P, mis_target, mut_mat_vec, out_file, 'exon', 1, effect, theta, m, number_muts, values, reps)  
            else: 
                if 'genes' in bin_type:
                    loop_exons_or_genes(P, mis_target, mut_mat_vec, out_file, 'gene', 0, effect, theta, m, number_muts, values, reps)
                elif 'exons' in bin_type :
                    loop_exons_or_genes(P, mis_target, mut_mat_vec, out_file, 'exon', 0, effect, theta, m, number_muts, values, reps)      
        elif suffix =='_01': 
            mis_target1 =  pd.read_csv(input_dir + 'Target_sizes/synonymous_' + bin_type_target +'_%s_1.csv.gz'%context , sep="\t", header=None, index_col = 0, compression='gzip')
            mis_target0 =  pd.read_csv(input_dir + 'Target_sizes/synonymous_' + bin_type_target +'_%s_0.csv.gz'%context , sep="\t", header=None, index_col = 0, compression='gzip')
            mut_matrix1 = pd.read_csv(input_dir + 'New_CBaSE_covariates/' + output_folder + '/mutation_matrix_' + name  + '_1.txt', sep = '\t', header = None).drop(columns = [0]).rename(columns = {1:0, 2:1, 3:2, 4:3})
            mut_mat_vec1, number_muts1, values1 = process_mutation_matrix(mut_matrix1, context) 
            mut_matrix0 = pd.read_csv(input_dir + 'New_CBaSE_covariates/' + output_folder + '/mutation_matrix_' + name  + '_0.txt', sep = '\t', header = None).drop(columns = [0]).rename(columns = {1:0, 2:1, 3:2, 4:3})
            mut_mat_vec0, number_muts0, values0 = process_mutation_matrix(mut_matrix0, context)
            if 'genes' in bin_type :
                loop_exons_or_genes(P[P['CpGI']==0], mis_target1, mut_mat_vec0, out_file, 'gene', 1, effect, theta, m, number_muts0, values0, reps)
                loop_exons_or_genes(P[P['CpGI']==1], mis_target0, mut_mat_vec1, out_file, 'gene', 0, effect, theta, m, number_muts1, values1, reps) 
            elif 'exons' in bin_type :
                loop_exons_or_genes(P[P['CpGI']==0], mis_target1, mut_mat_vec0, out_file, 'exon', 1, effect, theta, m, number_muts0, values0, reps)
                loop_exons_or_genes(P[P['CpGI']==1], mis_target0, mut_mat_vec1, out_file, 'exon', 0, effect, theta, m, number_muts1, values1, reps)         
        



#########################################################################################
#########################################################################################

# Input parameters
name = sys.argv[1] # e.g. ACC_3mers
context = name[-5:]
bin_type = str(sys.argv[2])
effect_considered = str(sys.argv[3]) # effect considered to compute the mutation matrix : synonymous, nonsynonymous or all

reps = 20
input_dir = '/nfs/users/dweghorn/projects/Selection_sites/'
dir = '/nfs/users/dweghorn/projects/Selection_sites/scripts/nextflow/results/' 

if effect_considered == 'synonymous':
    output_folder = 'Output_syn'
elif effect_considered == 'nonsynonymous':
    output_folder = 'Output_nonsyn'
else:
    output_folder = 'Output'

effect_simulated = "syn"

Q = pd.read_csv(dir + effect_simulated + '/'+ effect_considered + '/' + name + "/glm/fcovariates/fcovariates_genes.csv.gz", sep=',', compression='gzip')

list_of_genes = list(Q['gene'])

cov = 'cov'

if 'cpgi' in bin_type:
    simulate_selection_sites(name, context, '_0', bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)
    simulate_selection_sites(name, context, '_1',   bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)
    simulate_selection_sites(name, context, '_01', bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)
else: 
    simulate_selection_sites(name, context, '', bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)

if 'genes' in bin_type: 
    cov = 'cbase'
    if 'cpgi' in bin_type:
        simulate_selection_sites(name, context, '_0', bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)
        simulate_selection_sites(name, context, '_1',   bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)
        simulate_selection_sites(name, context, '_01', bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)
    else: 
        simulate_selection_sites(name, context, '', bin_type, cov, output_folder, reps, list_of_genes, effect_simulated)
