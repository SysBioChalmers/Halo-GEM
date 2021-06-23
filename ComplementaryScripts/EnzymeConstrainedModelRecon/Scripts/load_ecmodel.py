#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 6/23/21

"""load_ecmocel.py
:description : script
:param : 
:returns: 
:rtype: 
"""
import os
import sys
import pickle
import pandas as pd
import cobra
import itertools
import multiprocessing

ecpy_path = '../../../ecpy/'
sys.path.insert(0, os.path.abspath(ecpy_path))
from curate_kcat_light_goslim_asc_step2 import *

os.chdir('.')


# %%

def load_files():
    # ecGEM required files
    global posterior, dffrac, dftot, tmp_model
    posterior = pickle.load(open('../Results/smc_abc_light_go_and_one_asc_step_2.pkl', 'rb'))

    dffrac = pd.read_csv('../Results/protein_abundance_go_slim_level_uniq_asc.csv', index_col=0)
    dftot = pd.read_csv('../proteomics/total_protein_abandance_mean.csv', index_col=0)

    tmp_model_file = '../Results/template_ecModel_goslim_pools.pkl'
    tmp_model = pickle.load(open(tmp_model_file, 'rb'))
    tmp_model.solver = 'cplex'


def get_ecGEM(condition_id, random_index):
    global posterior, dffrac, dftot, tmp_model
    _, pools = buildPoolDict(dffrac, condition_id, dftot.loc[condition_id, 'Ptot'], sigma=0.5)
    kcat_dct = posterior.population[random_index]

    model = tmp_model.copy()
    updateKcats(model, kcat_dct)
    updateProteinConstraints(model, pools)
    return model


def test_growth(model):
    model.objective = 'Biomass_v1'
    model.objective_direction = 'max'
    s = model.optimize()
    r = s.objective_value if s.status == 'optimal' else -1
    print('Results:')
    print('  Status        :', s.status)
    print('  growth_rate   :', r)
    print('  glucose uptake:', s.fluxes['Exchange_Glucopyranose'])
    print('  PHA           :', s.fluxes['PHA_secretion'])
    print('  NGAM lb:', model.reactions.NGAM.lower_bound)
    print('  PHA lb:', model.reactions.PHA_secretion.lower_bound)


def pool_save_modelFiles(condition_id, random_index, outdir='../modelFiles/'):
    model = get_ecGEM(condition_id, random_index)
    file_name = '{}Halo_ecGEM_{}_{}'.format(outdir, condition_id, random_index)
    # cobra.io.write_sbml_model(model, file_name+'.xml')
    cobra.io.save_json_model(model, file_name + '.json')
    pickle.dump(model, open(file_name + '.pkl', 'wb'))


# %%
if __name__ == '__main__':
    global posterior, dffrac, dftot, tmp_model
    load_files()

    # all_conds = dftot.index.to_list()
    all_conds = ['NACL60', 'NACL20', 'NACL100', 'HN',
                 'Fermentation-9h', 'Fermentation-19h', 'Fermentation-30h',
                 'MU']
    condition_id = all_conds[np.random.randint(len(all_conds))]

    # len(posterior.population) == 100
    random_index = np.random.randint(100)

    # test get_ecGEM
    model = get_ecGEM(condition_id, random_index)
    test_growth(model)

    # test pool_save_modelFiles
    outdir = '../modelFiles/'
    os.system('mkdir {}'.format(outdir))
    pool_save_modelFiles(condition_id, random_index)

    # save model files
    cc = list(itertools.product(all_conds[0:1], [*range(0, 3)]))  # test,use all cores of CUP
    # cc = list(itertools.product(all_conds, [*range(0, 100)]))
    pool = multiprocessing.Pool()
    pool.starmap(pool_save_modelFiles, cc)
    pool.close()

    # test load model
    model_ = cobra.io.load_json_model(outdir + 'Halo_ecGEM_NACL60_3.json')
    test_growth(model_)

    model = get_ecGEM('NACL60', 3)
    test_growth(model)
