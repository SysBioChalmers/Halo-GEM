from matplotlib import pyplot as plt
import numpy as np
import cobra
import os

def report_model_rxns(model):
    for rxn in model.reactions: 
        print(rxn.id, rxn.name, rxn.gene_reaction_rule)
        print(rxn.reaction,'\n')
        
def report_model(model):
    print('Reactions:',len(model.reactions))
    print('Metabolites:',len(model.metabolites))
    
def print_rxn(rxn):
    eq = rxn.reaction
    for met in rxn.metabolites:
        eq = eq.replace(met.id,met.name)
    print(eq)
    
def correct_genename(model,gene_id,gene_name):
    gene = model.genes.get_by_id(gene_id)
    gene.name = gene_name
    
def print_reactions_of_gene(model,gene_id):
    gene = model.genes.get_by_id(gene_id)
    for rxn in gene.reactions: 
        print(rxn.gene_reaction_rule)
        print_rxn(rxn)
        print()
        
        
def plot_kcat_dist(df,col='kcat',title=None,logform=False):
    plt.figure(figsize=(4,3))
    if logform: plt.hist(df[col])
    else: plt.hist(np.log10(df[col]))
    plt.xlabel('log$_{10}(k_{cat}$ (s$^{-1}$))')
    plt.ylabel('Count')
    if title is not None: plt.title(title)
    plt.tight_layout()
    plt.show()
    
def plot_mw_dist(df):
    plt.figure(figsize=(4,3))
    plt.hist([df['mw']],20)
    plt.xlabel('MW (kDa)')
    plt.ylabel('Count')
    plt.tight_layout()
    plt.show()
    
    
def save_model_into_formats(model,model_dir,model_name):
    '''
    # save the model into three formats: json, mat and sbml
    
    Usage: save_model_into_formats(model,model_dir,model_name)
    
    Gang Li, last updated 2020-03-10
    '''
    xml_model = os.path.join(model_dir,'xml',model_name+'.xml')
    cobra.io.write_sbml_model(model, xml_model)
    print('Saved',xml_model)
    
    mat_model = os.path.join(model_dir,'mat',model_name+'.mat')
    cobra.io.save_matlab_model(model, mat_model)
    print('Saved',mat_model)
    
    json_model = os.path.join(model_dir,'json',model_name+'.json')
    cobra.io.save_json_model(model, json_model)
    print('Saved',json_model)
    print()
    
def report_match_kcats(case_count):
    keys = list(case_count.keys())
    keys.sort()
    
    operations = {'0': 'matching EC number and substrate',
                  '1': 'matching only EC number' }
    
    for key in keys:
        lst = key.split('_')
        print('With {0} wildcards, kcats were found for {1} reactions by {2}'.format(lst[0],case_count[key],operations[lst[1]]))

def test_biomass_production(model,biomass_id):
    model.objective = biomass_id
    model.objective_direction = 'max'
    print(model.optimize())