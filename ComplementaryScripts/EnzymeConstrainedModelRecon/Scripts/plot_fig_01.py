#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Created by Hao Luo at 1/7/21

"""plot_fig1.py
:description : script
:param : 
:returns: 
:rtype: 
"""

import itertools
import os
import pickle
import re
import sys

import matplotlib
import numpy as np
from matplotlib import pyplot as plt

ecpy_path = '../../../ecpy/'
sys.path.insert(0, os.path.abspath(ecpy_path))
from curate_kcat_light_goslim_asc_step2 import *

# %%md
### Fig1a from PPT

# %%
# Fig1b
# smc_abc_fitting_curve_0203
# ComplementaryScripts/EnzymeConstrainedModelRecon/Scripts/analyze_smc_abc.ipynb

results02 = pickle.load(open('../Results/smc_abc_light_step_02.pkl', 'rb'))
results03 = pickle.load(open('../Results/smc_abc_light_step_03.pkl', 'rb'))

# %%

# pickle.dump([results02.epsilons,results03.epsilons],open('../Results/validation_growth_rate_plot_data.pkl','wb'))

# [results02_epsilons,results03_epsilons] = pickle.load(open('../Results/validation_growth_rate_plot_data.pkl', 'rb'))

matplotlib.rc('font', family="Arial")

matplotlib.rcParams["font.family"] = 'Arial'  # 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial']


def plot_fitting_curve(results, ig=1, log=True, label=None, color=None):
    ax = plt.axes()
    ax.plot(np.arange(len(results.epsilons) - ig), results.epsilons[ig:], marker='o',
            markersize=3, label=label, color=color)
    # plt.scatter(np.arange(len(results.epsilons)-ig),results.epsilons[ig:])
    ax.set_xlabel('Iterations', fontname="Arial", fontsize=9)
    ax.set_ylabel('MSE', fontname="Arial", fontsize=9)
    if log: ax.set_yscale('log')
    # ax.tick_params(axis="y", direction="in")
    ax.tick_params(which="minor", axis="y", length=0, direction="in")
    # ax.legend(['Individual','GoSlim'],fontsize = 20)
    plt.setp(ax.texts, family='Arial')


plt.figure(figsize=(2.6, 2.4))
plot_fitting_curve(results03, ig=1, log=True, label='Individual', color='#91cf60')
plot_fitting_curve(results02, ig=1, log=True, label='GoSlim', color='#fc8d59')

plt.legend(prop={'family': 'Arial', 'size': 8})
# plt.tight_layout()
plt.savefig('../figures/fig_b_smc_abc_fitting_curve_0203.pdf', bbox_inches='tight')
plt.show()

# %%
# Fig1c

[val02, val03] = pickle.load(open('../Results/validation_growth_on_unused_conditions.pkl', 'rb'))
[val_prior02, val_prior03] = pickle.load(open('../Results/validation_growth_on_unused_conditions_prior.pkl', 'rb'))


# %%
def collect_data(res):
    res_with = {}
    res_without = {}

    for p in res:
        for k, lst in p.items():
            lst = lst.copy()
            for i in range(len(lst)):
                if lst[i] < 0: lst[i] = 0
            res_with[k] = res_with.get(k, []) + [lst[0]]
            try:
                res_without[k] = res_without.get(k, []) + [lst[1]]
            except:
                None

    return res_with, res_without


def plot_dot_bar(data, positions, label=None, color=None):
    data = np.array(data)
    data[data == -1] = 0
    means = np.mean(data, axis=1)

    plt.bar(positions, means, label=label, color=color)
    rep_x = []
    for item in positions: rep_x += [item for _ in range(data.shape[1])]
    p_rep_x = rep_x + np.random.normal(scale=0.1, size=len(rep_x))

    plt.scatter(p_rep_x, data.flatten(), s=2, zorder=3, color='gray', alpha=0.5)


# training data
res_train_02, _ = collect_data(results02.simulated_data)
res_train_03, _ = collect_data(results03.simulated_data)

res_train_prior_02, _ = collect_data(results02.simulated_data_t0)
res_train_prior_03, _ = collect_data(results03.simulated_data_t0)

train_ids = list(res_train_02.keys())
val_ids = list(val02.keys())
train_ids, val_ids

prior_data02 = [res_train_prior_02[item] for item in train_ids] + [val_prior02[item] for item in val_ids]
prior_data03 = [res_train_prior_03[item] for item in train_ids] + [val_prior03[item] for item in val_ids]
post_data02 = [res_train_02[item] for item in train_ids] + [val02[item] for item in val_ids]
post_data03 = [res_train_03[item] for item in train_ids] + [val03[item] for item in val_ids]

pos = np.arange(len(prior_data02)) * 5
plt.figure(figsize=(5.5, 2.7))
# plot_dot_bar(prior_data02,positions=pos,label='Prior02')
# plot_dot_bar(prior_data03,positions=pos+1,label='Prior03')
exp_data = [results02.datasets['dfpheno'].loc[item, 'SpecificGrowthRate'] for item in train_ids]
plt.bar(pos[:len(exp_data)], exp_data, label='Experiment')
plot_dot_bar(post_data03[:-1], positions=pos[:-1] + 1, label='Individual', color='#91cf60')
plot_dot_bar(post_data02[:-1], positions=pos[:-1] + 2, label='GoSlim', color='#fc8d59')

ylim = (0, 0.3)
plt.plot([pos[4] + 4] * 2, ylim, 'k--')

plt.ylim(ylim)
# plt.xlim((-1,7))
plt.xticks(pos[:-1] + 1, train_ids + val_ids[:-1], rotation=60, fontname="Arial", fontsize=9)
# plt.legend(bbox_to_anchor=(1.3, 1), ncol=1, loc='upper right', prop={'family': 'Arial', 'size': 8})
plt.legend(bbox_to_anchor=(1.0, -0.5), ncol=2, loc='lower right', prop={'family': 'Arial', 'size': 8})

plt.ylabel('Specific growth rate (h$^{-1}$)', fontname="Arial", fontsize=9)

plt.text(pos[2] + 0.5, 0.27, 'Training', ha='center', va='center')
plt.text(pos[6] - 0.5, 0.27, 'Validation', ha='center', va='center')
plt.savefig('../figures/fig_c_train_validation_02_03.pdf', bbox_inches='tight')
# plt.tight_layout()
plt.show()

# %% md

### 2. Compare Prior and Posterior


# %%

from scipy import stats


def f_test(x1, x2):
    F = np.var(x1) / np.var(x2)
    p = stats.f.cdf(F, len(x1) - 1, len(x2) - 1)
    return 1 - p


def extract_mean_var_inRV(results):
    def _format_key(k):
        return '{0}__{1}'.format(k[0], k[1])

    # kcat_model02 = [item.loc for item in results02.posterior.values()]  # those values are log10 transformed 1/h
    # kcat_model02 = [10 ** item / 3600 for item in kcat_model02]  # convert it to 1/s

    df = pd.DataFrame()
    for k, rv in results.prior.items():
        k = _format_key(k)
        df.loc[k, 'prior_mean'] = np.log10(10 ** rv.loc / 3600)  # rv.loc
        df.loc[k, 'prior_std'] = np.log10(10 ** rv.scale / 3600)  # rv.scale
        # loc = np.log10(10**x*3600)
        # scale = np.log10(10**x*3600)

    for k, rv in results.posterior.items():
        k = _format_key(k)
        df.loc[k, 'post_mean'] = np.log10(10 ** rv.loc / 3600)  # rv.loc
        df.loc[k, 'post_std'] = np.log10(10 ** rv.scale / 3600)

    print(df.shape)
    print(df.head(n=5))

    return df  # in the unit of 1/s log10-transformed


def plot_mean_scatter_comp(df):
    def _scatter_plot(prior, post, lim, title):
        # plt.figure(figsize=(3, 2.5))
        ax = plt.figure(figsize=(2.6, 2.4)).gca()
        # ax.yaxis.set_major_locator(MaxNLocator(integer=True))
        plt.scatter(prior, post, zorder=1, edgecolor='#5ab4ac',
                    facecolor='none',
                    s=20,
                    label=title,
                    alpha=0.8)

        plt.xlabel('Prior log10($k_{cat}$) ($s^{- 1}$)', fontname="Arial", fontsize=9)
        plt.ylabel('Posterior log10($k_{cat}$) ($s^{- 1}$)', fontname="Arial", fontsize=9)
        # plt.yticks([-6,-4,-2,0,2,4,6])
        # plt.xticks([-6,-4,-2,0,2,4,6])

        # lim = [-6,7]
        plt.plot(lim, lim, zorder=0, color='k')

        plt.xlim(lim)
        plt.ylim(lim)
        # plt.title(title)
        # plt.legend(title,prop={'family': 'Arial', 'size': 8})
        plt.legend(prop={'family': 'Arial', 'size': 8})
        # plt.tight_layout()
        plt.savefig('../figures/fig_de_Compare_Prior_and_Posterior_' + title + '_.pdf', bbox_inches='tight')
        plt.show()

    _scatter_plot(df['prior_mean'], df['post_mean'], [-8, 8], 'mean')
    _scatter_plot(df['prior_std'], df['post_std'], [-4, 4], 'std')


def visualize(results):
    dfmstd = extract_mean_var_inRV(results)
    # plot_mean_dist_comp(dfmstd)
    plot_mean_scatter_comp(dfmstd)


visualize(results02)


# %%
# Fig1f

def load_post_and_ecoli(results02, results03):
    kcat_model02 = [item.loc for item in results02.posterior.values()]  # those values are log10 transformed 1/h
    kcat_model02 = [10 ** item / 3600 for item in kcat_model02]  # convert it to 1/s

    kcat_model03 = [item.loc for item in results03.posterior.values()]  # those values are log10 transformed 1/h
    kcat_model03 = [10 ** item / 3600 for item in kcat_model03]  # convert it to 1/s

    # load prior, do not use results.prior, as it was not the real prior. forgot to use .copy in the very begining
    model_one_pool_file = '../Results/template_ecModel_one_pool.pkl'
    tmp_ecmodel = pickle.load(open(model_one_pool_file, 'rb'))
    df_enz_kcat = pd.read_csv('../Results/mapped_kcats_updated_with_ko.csv', index_col=0)

    priors = build_priors(tmp_ecmodel, df_enz_kcat)
    kcat_model_prior = [item.loc for item in priors.values()]  # those values are log10 transformed 1/h
    kcat_model_prior = [10 ** item / 3600 for item in kcat_model_prior]  # convert it to 1/s

    # load ecoli kcats
    kcatDB = pd.read_csv('../../../ComplementaryData/max_KCAT.txt', sep='\t', header=None)
    ecoli_kcat = kcatDB[kcatDB[2].str.contains('escherichia coli')][3].values

    X02 = np.sort(kcat_model02)
    Y02 = np.array(range(len(X02))) / float(len(X02))

    X03 = np.sort(kcat_model03)
    Y03 = np.array(range(len(X03))) / float(len(X03))

    X_model_prior = np.sort(kcat_model_prior)
    Y_model_prior = np.array(range(len(X_model_prior))) / float(len(X_model_prior))

    X_ecoli = np.sort(ecoli_kcat)
    Y_ecoli = np.array(range(len(X_ecoli))) / float(len(X_ecoli))

    # X_bac = np.sort(bac_kcat)
    # Y_bac = np.array(range(len(X_bac)))/float(len(X_bac))

    # X_all = np.sort(kcatDB[3].values)
    # Y_all = np.array(range(len(X_all)))/float(len(X_all))

    return X02, Y02, X03, Y03, X_model_prior, Y_model_prior, X_ecoli, Y_ecoli, kcat_model_prior


def do_plot_compare_to_ecoli(X02, Y02, X03, Y03, X_model_prior, Y_model_prior):
    plt.figure(figsize=(2.7, 2.5))
    # ax = plt.figure(figsize=(3, 2.7)).gca()

    plt.plot(X_model_prior, Y_model_prior, label='$Prior$ $H$. TD01')
    plt.plot(X02, Y02, label='$Posterior$ (GoSlim)')  #
    plt.plot(X03, Y03, label='$Posterior$ (Individual)')  #
    plt.plot(X_ecoli, Y_ecoli, label='$E. coli$')

    # plt.plot(X_all  , Y_all  ,label = 'All')

    plt.xscale('log')
    plt.legend(loc=0, prop={'family': 'Arial', 'size': 6.7})
    plt.xlabel('$k_{cat}$ ($s^{- 1}$)', fontname="Arial", fontsize=9)
    plt.ylabel('Cumulative distribution', fontname="Arial", fontsize=9)
    # plt.tight_layout()
    plt.savefig('../figures/fig_f_prior_post_kcat_aculative_ecoli.pdf', bbox_inches='tight')
    plt.show()
    plt.close()


X02, Y02, X03, Y03, X_model_prior, Y_model_prior, X_ecoli, Y_ecoli, kcat_model_prior = load_post_and_ecoli(results02,
                                                                                                           results03)

do_plot_compare_to_ecoli(X02, Y02, X03, Y03, X_model_prior, Y_model_prior)

results = results02.population.copy()  # TODO: populations are kcat?
kcat_model = np.array(list(results[0].values())) / 3600  # TODO: there are 100 sets kcat?

# loda table and get 'Genes' column
go_slim_db = pd.read_excel('../proteomics/go-slim-in model.xlsx', index_col=None, header=0)
go_slim_db['Genes'] = [list(set(i[4:]) - {np.nan}) for i in go_slim_db.values.tolist()]  # remove np.nan
go_slim_db = go_slim_db[['Description', 'GoTerm', 'GeneNumber', 'Classification', 'Genes']]
go_slim_db.head

# get 'Kcats' column for 'Genes' column
Gens_id_list = np.array([re.sub('(prot_)|(_(c|p|e))', '', key[1]) for key in results02.posterior.keys()])
kcat_model = np.array(kcat_model)
kcat_model_prior = np.array(kcat_model_prior)

kcats_column = []
kcats_prior_column = []
for index in go_slim_db.index:
    genes_row = go_slim_db.iloc[index]['Genes']
    kcats_row = []
    kcats_prior_row = []
    for gene_i in genes_row:
        kcats_i = kcat_model[Gens_id_list == gene_i]
        kcats_row.append(list(kcats_i))
        kcats_prior_i = kcat_model_prior[Gens_id_list == gene_i]
        kcats_prior_row.append(list(kcats_prior_i))
    kcats_column.append(kcats_row)
    kcats_prior_column.append(kcats_prior_row)
go_slim_db['Kcats'] = kcats_column
go_slim_db['Kcats_prior'] = kcats_prior_column
go_slim_db.head()
# plot for each classification
go_slim_db['Kcats_all'] = go_slim_db['Kcats'].apply(itertools.chain.from_iterable)  # convert to One-dimensional list
go_slim_db['Kcats_all'] = go_slim_db['Kcats_all'].apply(list)

go_slim_db['Kcats_prior_all'] = go_slim_db['Kcats_prior'].apply(itertools.chain.from_iterable)
go_slim_db['Kcats_prior_all'] = go_slim_db['Kcats_prior_all'].apply(list)

descriptions = ['cellular respiration', 'ion transport', 'ion transmembrane transport',
                'nucleic acid metabolic process', 'transcription, DNA-dependent', 'ribosome assembly',
                'alpha-amino acid biosynthetic process', 'translation', 'response to osmotic stress',
                'ATP biosynthetic process', 'oxidation-reduction process', 'tricarboxylic acid cycle']
# %%
plt.close()

plt.figure(figsize=(2.7, 2.5))
# figure = plt.gcf()
# figure.set_size_inches(5, 2.7)
# ax = plt.figure(figsize=(5, 2.7)).gca()
for description_i in descriptions:
    kcat_i = go_slim_db[go_slim_db['Description'] == description_i]['Kcats_all'].values
    kcat_i = list(itertools.chain.from_iterable(kcat_i))
    # print(len(kcat_i),description_i)
    if len(kcat_i) == 0: continue
    X_i = np.sort(kcat_i)
    Y_i = np.array(range(len(kcat_i))) / float(len(kcat_i))
    plt.plot(X_i, Y_i, label=description_i + ' ({0})'.format(len(kcat_i)))
    plt.xscale('log')
plt.legend(bbox_to_anchor=(1, 0.9), prop={'family': 'Arial', 'size': 6.7})
plt.xlabel('$k_{cat}$ ($s^{- 1}$)', fontname="Arial", fontsize=9)
plt.ylabel('Cumulative distribution', fontname="Arial", fontsize=9)

plt.savefig('../figures/fig_g_post_kcat_aculative_go_terms_3.pdf', bbox_inches='tight')

plt.show()

plt.close()
