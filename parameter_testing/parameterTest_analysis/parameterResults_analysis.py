# run on ipython
import pandas as pd
import os

tool = 'kaiju'
run_nums = [34, 35, 36, 37, 38, 39]
result_filename = 'Algae_m3-15results.txt'

wdir = '/home/ae42909/Scratch/parameter_test/' + tool + '/'
runs = ['run' + str(nums) for nums in run_nums]


sample_names = []
parameter_names = []
parameter_values = []
find_tag = []
for run in runs:
    run_dir = os.path.join(wdir, run)
    for root, subdirs, files in os.walk(run_dir):
        for filename in files:
            if 'log_file_' in filename:
                tag = filename.split('log_file_',1)[1]
                find_tag.append(tag)
                with open(os.path.join(root, filename), 'r') as in_file:
                    for lines in in_file:
                        if 'input1 =' in lines:
                            sample_names.append(lines.rpartition('/')[2].rstrip())
                        if 'expected virus score =' in lines:
                                parameter_values.append(round(float(lines.partition('expected virus score = ')[2].rstrip()),2))
                        if tool == 'kraken':
                            if 'kraken database =' in lines:
                                parameter_names.append('k-mer_' + lines.partition('krakenDB_k')[2].rpartition('_m')[0])
                        elif tool == 'kaiju':
                            if 'kaiju_minlen =' in lines:
                                kaiju_parameter = 'mL_' + lines.partition('kaiju_minlen = ')[2].rstrip()
                            if 'kaiju_mismatch =' in lines:
                                kaiju_parameter += '/mM_' + lines.partition('kaiju_mismatch =')[2].rstrip()
                            if 'kaiju_score =' in lines:
                                kaiju_parameter += '/s_' + lines.partition('kaiju_score =')[2].rstrip()
                                parameter_names.append(kaiju_parameter)

                        
                            
    

results = pd.DataFrame(columns=set(sample_names), index=set(parameter_names))
for position, values in enumerate(parameter_values):
    col = list(results.columns).index(sample_names[position])
    row = list(results.index).index(parameter_names[position])
    results.ix[row,col] = values
results.to_csv(wdir + result_filename, sep='\t')

which_tag = pd.DataFrame(columns=set(sample_names), index=set(parameter_names))
for position, values in enumerate(find_tag):
    col = list(which_tag.columns).index(sample_names[position])
    row = list(which_tag.index).index(parameter_names[position])
    which_tag.ix[row,col] = values
which_tag.to_csv(wdir + result_filename[:-4] + '_tag', sep='\t')

