#!/usr/bin/env python
# coding: utf-8

# In[1]:


import sys

import numpy as np
import pandas as pd
from cobra.flux_analysis import flux_variability_analysis

sys.path.insert(0, '/Users/lhao/Documents/Git/py_Halo-GEM/ecpy')
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
dffrac = dffrac_pro
result_df_pro = pd.DataFrame()
for kcat in results02.population:
    # eModel = pickle.load(open(eModel_file,'rb'))
    model = model_go_pools.copy()

    updateKcats(model, kcat)

    with model:

        for condition_id in test_conds:

            print('Condition:', condition_id)
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

                model.reactions.PHA_secretion.lower_bound = 0
                model.reactions.NGAM.lower_bound = 0

                s = model.optimize()
                r2 = s.objective_value if s.status == 'optimal' else -1

            except:
                r1, r2 = -1, -1
            print(r1, r2)

            for rea_i in model.reactions:
                if rea_i.upper_bound > 9999:
                    rea_i.upper_bound = 9999
            # fva_realist = [i for i in model.reactions if 'prot_' not in i.id]
            try:
                result_df_i_pro = flux_variability_analysis(model, model.reactions[:3512], fraction_of_optimum=0.95)
                df = pd.DataFrame(data={'minimum': [r1], 'maximum': [r2]}, index=['growth'])
                result_df_i_pro = result_df_i_pro.append(df)

                result_df_i_pro.columns = condition_id + '_' + result_df_i_pro.columns
                result_df_pro = result_df_pro.merge(result_df_i_pro, how='outer', left_index=True, right_index=True, )
                print(result_df_pro)
            except:
                continue
    result_df_pro.to_csv('../Results/FVA_result_df_pro.tsv', sep='\t')

# %%
dffrac = dffrac_mRNA
result_df_mRNA = pd.DataFrame()
for kcat in results02.population:
    # eModel = pickle.load(open(eModel_file,'rb'))
    model = model_go_pools.copy()

    updateKcats(model, kcat)

    with model:

        for condition_id in test_conds:

            print('Condition:', condition_id)
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

                model.reactions.PHA_secretion.lower_bound = 0
                model.reactions.NGAM.lower_bound = 0

                s = model.optimize()
                r2 = s.objective_value if s.status == 'optimal' else -1

            except:
                r1, r2 = -1, -1
            print(r1, r2)

            for rea_i in model.reactions:
                if rea_i.upper_bound > 9999:
                    rea_i.upper_bound = 9999
            # fva_realist = [i for i in model.reactions if 'prot_' not in i.id]
            try:
                result_df_i_pro = flux_variability_analysis(model, model.reactions[:3512], fraction_of_optimum=0.95)
                result_df_i_pro.columns = condition_id + '_' + result_df_i_pro.columns

                df = pd.DataFrame(data={'minimum':[r1], 'maximum':[r2]},index=['growth'])
                result_df_pro = result_df_pro.append(df)
                result_df_mRNA = result_df_mRNA.merge(result_df_i_pro, how='outer', left_index=True, right_index=True, )
                print(result_df_mRNA)
            except:
                raise
    result_df_mRNA.to_csv('../Results/FVA_result_df_mRNA.tsv', sep='\t')
