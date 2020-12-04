#!/usr/bin/env python
# coding: utf-8

# In[1]:


import multiprocessing
import os
import sys
import time

import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

ecpy_path = '../../../ecpy/'
sys.path.insert(0, os.path.abspath(ecpy_path))
import pickle

# from curate_kcat import *

# from multiprocessing import Process,cpu_count,Manager


# In[4]:
# %%

from curate_kcat_light_goslim_asc_step2 import *

results02 = pickle.load(open('../Results/smc_abc_light_step_02.pkl', 'rb'))
dfpheno = pd.read_csv('../proteomics/phynotype.csv', index_col=0, comment='#')
dffrac_pro = pd.read_csv('../Results/protein_abundance_go_slim_level_uniq_asc.csv', index_col=0)
dffrac_mRNA = pd.read_excel('../proteomics/mRNA_abundance_go_slim-20201124.xlsx', sheet_name='mRNA', index_col=0,
                            header=0)
mass_columns = [i for i in dffrac_mRNA.columns if 'Mass' in i]
dffrac_mRNA[mass_columns] = dffrac_mRNA[mass_columns]*1.651567

dftot = pd.read_csv('../proteomics/total_protein_abandance_mean.csv', index_col=0)

# %%

model_go_pools_file = '../Results/template_ecModel_goslim_pools.pkl'
model_go_pools = pickle.load(open(model_go_pools_file, 'rb'))
model_go_pools.objective = 'Biomass_v1'
model_go_pools.objective_direction = 'max'

# %%
test_conds = ['NACL60', 'NACL20', 'NACL100', 'HN', 'Fermentation-9h', 'Fermentation-19h', 'Fermentation-30h']
for cond_id in dftot.index:
    if cond_id not in dfpheno.index:
        dfpheno.loc[cond_id, 'P3HB'] = 0
print()


# %%

def fva_special(model_go_pools, results02, test_conds, dftot, dffrac, output_file):
    result_df = pd.DataFrame()
    result_df_growth = pd.DataFrame(columns=['growth'])
    kcat_number = -1
    for kcat in results02.population:
        # eModel = pickle.load(open(eModel_file,'rb'))
        model = model_go_pools.copy()
        kcat_number += 1

        updateKcats(model, kcat)

        with model:

            for condition_id in test_conds:

                print('\nkcat index:', kcat_number, '\nCondition:', condition_id, )
                Ptot = dftot.loc[condition_id, 'Ptot']
                one_pool, go_pools = buildPoolDict(dffrac, condition_id, Ptot)
                pools = go_pools
                updateProteinConstraints(model, pools)

                PHA = dfpheno.loc[condition_id, 'P3HB']
                if ~np.isnan(PHA):
                    model.reactions.PHA_secretion.lower_bound = PHA

                try:
                    s = model.optimize()
                    r1 = s.objective_value if s.status == 'optimal' else -1

                    # model.reactions.PHA_secretion.lower_bound = 0
                    # model.reactions.NGAM.lower_bound = 0
                    #
                    # s = model.optimize()
                    # r2 = s.objective_value if s.status == 'optimal' else -1

                except:
                    r1 = -1
                print('with pha ngam: %f \n' % (r1))
                df = pd.DataFrame(data={'growth': [r1]}, index=[str(kcat_number) + '_' + condition_id])
                result_df_growth = result_df_growth.append(df)

                for rea_i in model.reactions:
                    if rea_i.upper_bound > 9999:
                        rea_i.upper_bound = 9999
                # fva_realist = [i for i in model.reactions if 'prot_' not in i.id]
                try:
                    if r1 != -1:
                        result_df_i_pro = flux_variability_analysis(model, model.reactions,
                                                                    fraction_of_optimum=0.95)
                        result_df_i_pro.columns = str(kcat_number) + '_' + condition_id + '_' + result_df_i_pro.columns
                        result_df = result_df.merge(result_df_i_pro, how='outer', left_index=True, right_index=True, )
                        # print(result_df)
                except:
                    continue
            result_df.to_csv(output_file, sep='\t')
            result_df_growth.to_csv(output_file.replace('.tsv', '_growth.tsv'), sep='\t')


# %%
start = time.time()

p1 = multiprocessing.Process(target=fva_special, args=(
    model_go_pools, results02, test_conds, dftot, dffrac_pro, '../Results/FVA_result_df_pro_1204.tsv'))

p2 = multiprocessing.Process(target=fva_special, args=(
    model_go_pools, results02, test_conds, dftot, dffrac_mRNA, '../Results/FVA_result_df_mRNA_1204.tsv'))
p1.start()
p2.start()

p1.join()
p2.join()

# fva_special(model_go_pools, results02, test_conds, dftot, dffrac=dffrac_pro,
#             output_file='../Results/FVA_result_df_pro_1204.tsv')
#
#
# fva_special(model_go_pools, results02, test_conds, dftot, dffrac=dffrac_mRNA,
#             output_file='../Results/FVA_result_df_mRNA_1204.tsv')


print("--- %s seconds ---" % (time.time() - start))
