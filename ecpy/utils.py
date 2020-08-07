from matplotlib import pyplot as plt
import numpy as np
import cobra
import os
import ecpy

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

def test_biomass_production(model,show=False):
    with model:
        set_yeast_extraction(model,lb=0,ub=0)
        model.objective = 'Biomass_v1'
        model.objective_direction = 'max'
        s1 = model.optimize()
        if show:print('Without YE:',s1)

    with model:
        set_yeast_extraction(model,lb=0,ub=1000)
        model.objective = 'Biomass_v1'
        model.objective_direction = 'max'
        s2 = model.optimize()
        if show:print('With unlimited YE:',s2)
    
    if not show:return s1,s2

def test_atp_production(model):
    atp_sink = cobra.Reaction('atp_sink')
    atp_sink.add_metabolites({model.metabolites.ATP_c:-1})
    
    with model:
        model.add_reactions([atp_sink])
        set_yeast_extraction(model,lb=0,ub=0)
        model.objective = 'atp_sink'
        model.objective_direction = 'max'
        s1 = model.optimize()

    with model:
        model.add_reactions([atp_sink])
        set_yeast_extraction(model,lb=0,ub=1000)
        model.objective = 'atp_sink'
        model.objective_direction = 'max'
        s2 = model.optimize()
    return s1,s2

def test_ngam_flux(model):
    with model:
        set_yeast_extraction(model,lb=0,ub=0)
        set_bound(model,'NGAM',lb=-1000,ub=1000)
        model.objective = 'NGAM'
        model.objective_direction = 'max'
        s1 = model.optimize()

    with model:
        set_yeast_extraction(model,lb=0,ub=1000)
        set_bound(model,'NGAM',lb=-1000,ub=1000)
        
        model.objective = 'NGAM'
        model.objective_direction = 'max'
        s2 = model.optimize()
    return s1,s2
    


def set_yeast_extraction(model,ub=1000,lb=-1000):
    # given the evidence that TD01 strain can growth wihout supplyment of amino acids, 
    # the strain should be able to produce 20 amino acids by itself.
    # Block the uptake of 20 amino acids
    amino_acids = ['CYS_c','L-ASPARTATE_c','GLT_c','PHE_c','GLY_c',
                  'HIS_c','ILE_c','LYS_c','LEU_c','MET_c','ASN_c',
                  'PRO_c','GLN_c','ARG_c','SER_c','THR_c','VAL_c',
                  'TRP_c','L-ALPHA-ALANINE_c','TYR_c','L-LACTATE','CPD-15990', 'ADENINE']
    for met_id in amino_acids:
        exg_rxn = model.reactions.get_by_id('Exchange_'+met_id.replace('_c',''))
        exg_rxn.upper_bound = ub
        exg_rxn.lower_bound = lb
     

def report_gene_status(model):
    print('Genes')
    print('  Number of genes:',len(model.genes))
    
    # number of protein complexes
    complexes = []
    genes_in_complex = []
    for rxn in model.reactions:
        cmplx = ecpy.parse_gr_rule(rxn.gene_reaction_rule)
        for item in cmplx:
            if 'and' in item: 
                complexes.append(item)
                genes_in_complex += [g.strip() for g in item.split('and')]
    print('  Number of Multigene complexes:',len(set(complexes)))
    print('  Number of Genes involved in complexes:',len(set(genes_in_complex)))

def is_metabolic_rxn(rxn):
    # with different mets in two sides
    pass
    
def report_rxn_status(model):
    print('Reactions')
    print('  Number of reactions:',len(model.reactions))
    #met_rxns = [rxn for rxn in model.reactions if ]
def report_gem_status(model):
    report_gene_status(model)
    
def set_bound(model,rxn_id,lb=None,ub=None,eq=None):
    rxn = model.reactions.get_by_id(rxn_id)
    if lb is not None: rxn.lower_bound = lb
    if ub is not None: rxn.upper_bound = ub
    if eq is not None: 
        rxn.lower_bound = eq
        rxn.upper_bound = eq
def test_met_production(model,met_id):
    with model:
        sink = cobra.Reaction('tmp_sink')
        sink.add_metabolites({model.metabolites.get_by_id(met_id):-1})
        model.add_reactions([sink])
        model.objective = 'tmp_sink'
        s = model.optimize()
        print(s)
    return s
        
def test_glucose(model):
    rgs = []
    glcs = np.arange(0,20)
    for glc in glcs:
        set_bound(model,'Exchange_Glucopyranose',ub=glc)
        s1,s2 = test_biomass_production(model,show=False)
        rgs.append([s1.objective_value,s2.objective_value])
    
    plt.figure(figsize=(7,3))
    rgs = np.array(rgs)
    titles = ['Without YE','With YE']
    for i in range(2):
        plt.subplot(1,2,i+1)
        #plt.scatter(glcs,rgs[:,0],label='Without YE')
        plt.scatter(glcs,rgs[:,i])
        plt.ylabel('Specific growth rate (h$^{-1}$)')
        plt.xlabel('Glucose uptake (mmol/dDW/h)')
        plt.title(titles[i])
        
    plt.tight_layout()
    plt.show()

def test_NGAM(model):
    rgs = []
    ngams = np.arange(0,20)
    for ngam in ngams:
        set_bound(model,'NGAM',eq=ngam)
        s1,s2 = test_biomass_production(model)
        rgs.append([s1.objective_value,s2.objective_value])
    
    plt.figure(figsize=(7,3))
    rgs = np.array(rgs)
    titles = ['Without YE','With YE']
    for i in range(2):
        plt.subplot(1,2,i+1)
        #plt.scatter(glcs,rgs[:,0],label='Without YE')
        plt.scatter(ngams,rgs[:,i])
        plt.title(titles[i])
        plt.ylabel('Specific growth rate (h$^{-1}$)')
        plt.xlabel('NGAM (mmol/dDW/h)')
        
    #plt.ylim(0,5)
    plt.tight_layout()
    plt.show()
    
def test_PHA(model,ys=np.arange(0,1000,20)):
    rgs = []
    for y in ys:
        set_bound(model,'PHA_secretion',lb=y,ub=1000)
        s1,s2 = test_biomass_production(model)
        rgs.append([s1.objective_value,s2.objective_value])
        
    plt.figure(figsize=(7,3))
    rgs = np.array(rgs)
    titles = ['Without YE','With YE']
    for i in range(2):
        plt.subplot(1,2,i+1)
        #plt.scatter(glcs,rgs[:,0],label='Without YE')
        plt.scatter(ys,rgs[:,i])
        plt.ylabel('Specific growth rate (h$^{-1}$)')
        plt.xlabel('PHA secretion (mmol/dDW/h)')
        plt.title(titles[i])
        
   
    #plt.ylim(0,5)
    plt.tight_layout()
    plt.show()
    
def test_Glc_to_ATP(model):
    rgs = []
    glcs = np.arange(0,20)
    for glc in glcs:
        set_bound(model,'Exchange_Glucopyranose',ub=glc)
        s1,s2 = test_atp_production(model)
        rgs.append([s1.objective_value,s2.objective_value])
        
    plt.figure(figsize=(7,3))
    rgs = np.array(rgs)
    titles = ['Without YE','With YE']
    for i in range(2):
        plt.subplot(1,2,i+1)
        #plt.scatter(glcs,rgs[:,0],label='Without YE')
        plt.scatter(glcs,rgs[:,i])
        plt.title(titles[i])
        plt.ylabel('ATP_c (mmol/dDW/h)')
        plt.xlabel('Glucose uptake (mmol/dDW/h)')
    
    plt.tight_layout()
    plt.show()
    
    
def change_rxn_coeff(rxn,met,new_coeff):
    '''
    # This is based on the rxn.add_metabolites function. If there the metabolite is already in the reaction,
    # new and old coefficients will be added. For example, if the old coeff of metA is 1, use
    # rxn.add_metabolites({metA:2}), After adding, the coeff of metA is 1+2 = 3
    #
    '''

    diff_coeff = new_coeff-rxn.metabolites[met]
    rxn.add_metabolites({met:diff_coeff})
    
def test_GAM_on_Growth(model):
    rgs = []
    
    gams = np.arange(40,80)
    atp_c = model.metabolites.ATP_c
    biomass = model.reactions.Biomass_v1
    for gam in gams:
        change_rxn_coeff(biomass,atp_c,-gam)
        s1,s2 = test_atp_production(model)
        rgs.append([s1.objective_value,s2.objective_value])
        
    plt.figure(figsize=(7,3))
    rgs = np.array(rgs)
    titles = ['Without YE','With YE']
    for i in range(2):
        plt.subplot(1,2,i+1)
        #plt.scatter(glcs,rgs[:,0],label='Without YE')
        plt.scatter(gams,rgs[:,i])
        plt.title(titles[i])
        plt.ylabel('Specific growth rate (h$^{-1}$)')
        plt.xlabel('GAM (mmol/gDW)')
    
    plt.tight_layout()
    plt.show()

def test_Glc_to_ATP_ecoli(eco_model):
    glcs = np.arange(0,20)
    rgs = []
    atp_sink = cobra.Reaction('atp_sink')
    atp_sink.add_metabolites({eco_model.metabolites.atp_c:-1})
    
    with eco_model as model:
        model.add_reactions([atp_sink])
        for glc in glcs:
            with model as m:
                m.reactions.EX_glc__D_e.lower_bound = -glc
                m.reactions.ATPM.lower_bound = 0
                m.objective = 'atp_sink'

                m.objective_direction = 'max'
                rgs.append(m.optimize().objective_value)
                
    plt.figure(figsize=(4,3))
    plt.scatter(glcs,rgs)
    plt.ylabel('ATP (mmol/dDW/h)')
    plt.xlabel('Glucose uptake (mmol/dDW/h)')
    print(rgs)
    plt.tight_layout()
    plt.show()
    
    
def test_Glc_to_NGAM(model):
    rgs = []
    glcs = np.arange(0,20)
    for glc in glcs:
        set_bound(model,'Exchange_Glucopyranose',ub=glc)
        s1,s2 = test_ngam_flux(model)
        rgs.append([s1.objective_value,s2.objective_value])
    print(rgs)
    plt.figure(figsize=(7,3))
    rgs = np.array(rgs)
    titles = ['Without YE','With YE']
    for i in range(2):
        plt.subplot(1,2,i+1)
        #plt.scatter(glcs,rgs[:,0],label='Without YE')
        plt.scatter(glcs,rgs[:,i])
        plt.title(titles[i])
        plt.ylabel('NGAM (mmol/dDW/h)')
        plt.xlabel('Glucose uptake (mmol/dDW/h)')
    
    plt.tight_layout()
    plt.show()

def test_Glc_to_NGAM_ecoli(eco_model):
    glcs = np.arange(0,20)
    rgs = []
    with eco_model as model:
        for glc in glcs:
            with model as m:
                m.reactions.EX_glc__D_e.lower_bound = -glc
                m.reactions.ATPM.lower_bound = 0
                m.objective = 'ATPM'

                m.objective_direction = 'max'
                rgs.append(m.optimize().objective_value)
    plt.figure(figsize=(4,3))
    plt.scatter(glcs,rgs)
    plt.ylabel('NGAM (mmol/dDW/h)')
    plt.xlabel('Glucose uptake (mmol/dDW/h)')
    print(rgs)
    plt.tight_layout()
    plt.show()
    
    
def metacyc_id_indexed_ecoli_rxns(ecoli_model):
    met_rxns = dict()
    for rxn in ecoli_model.reactions:
        # ['META:3.6.3.10-RXN', 'META:3.6.3.12-RXN', 'META:TRANS-RXN-2']
        metcyc = rxn.annotation.get('biocyc') 
        if metcyc is None: continue
        elif type(metcyc) is list:
            for item in metcyc: met_rxns[item.replace('META:','')] = rxn
        else: met_rxns[metcyc.replace('META:','')] = rxn
    return met_rxns