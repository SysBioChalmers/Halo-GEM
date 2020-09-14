#!/usr/bin/env python
# coding: utf-8

# ### Constrain based on go-term

# In[1]:


import cobra
import os 
import sys
import numpy as np
import pandas as pd
import time
from sklearn.metrics import mean_squared_error
import scipy.stats as ss
from multiprocessing import Process,cpu_count,Manager
import pickle


# In[2]:


ecpy_path = '../../../ecpy/'
sys.path.append(os.path.abspath(ecpy_path))
import utils
import ecpy


# import importlib
# importlib.reload(utils)
# importlib.reload(ecpy)

# In[3]:


class RV:
    def __init__(self,loc,scale):
        '''
        Normal distribution assumed
        loc, scale: the same as used in scipy.stats. kcat in log10(1/h)
        
        '''
       
        
        self.loc = loc
        self.scale = scale
        self.lb = np.log10(1e-5*3600) # kcat in range 1e-5, 1e6 1/s
        self.ub = np.log10(1e6*3600)
        
        
        
    def sample(self):
        '''
        Generate a random sample from the given prior distribution
        '''
        kcat = np.random.normal(self.loc,self.scale)
        if kcat < self.lb: kcat = self.lb
        if kcat > self.ub: kcat = self.ub
        return kcat


# In[4]:


def initialize_kcats_dict(tmp_ecmodel):
    # Initialized kcats_dict, h-1, values do not matter since only keys will be used 
    initial_kcats = {}
    for met in tmp_ecmodel.metabolites:
        # look for those metabolites: prot_protid
        if not met.id.startswith('prot_'): continue

        # ingore metabolite: prot_pool
        if met.id.startswith('prot_pool'): continue

        prot_id = met.id.split('_')[1]
        for rxn in met.reactions:
            if rxn.id.startswith('draw_prot'): continue
            initial_kcats[(rxn.id,met.id)] = -1/rxn.metabolites[met]
    return initial_kcats


# In[5]:


def build_priors(tmp_ecmodel,df_enz_kcat):
    # tmp_ecmodel, an ecModel with all enzymes and pools
    # in df_enz_kcat, value are in the unit of 1/s, log1o-transformed
    # in the resulting prior, kcat is in log10-transformed, 1/h
    
    initial_kcats = initialize_kcats_dict(tmp_ecmodel)
            
    # Build priors
    priors = dict()
    for rxn_id,prot_met_id in initial_kcats.keys():
        rxn_short_id = rxn_id.split('No')[0]
        loc = np.log10(10**df_enz_kcat.loc[rxn_short_id,'log10_kcat_mean']*3600)
        std_in_log_seconds = df_enz_kcat.loc[rxn_short_id,'log10_kcat_std']
 
        scale = np.log10(10**std_in_log_seconds*3600)
        priors[(rxn_id,prot_met_id)] = RV(loc,scale)
    return priors


# In[6]:


def updateProteinConstraints(tmp_ecmodel,pools):
    # pools, a dictionary, {'prot_pool_exchange_GO0006457': ub, ...}
    for rxn_id, ub in pools.items():
        rxn = tmp_ecmodel.reactions.get_by_id(rxn_id)
        rxn.upper_bound = ub


# In[7]:


def buildPoolDict(dffrac,condition_id,Ptot,sigma=0.5):
    # usego = True  : build go based pools
    # usego = False : build only one pool
    one_pool = dict()
    go_pools = dict()
    
    one_pool['prot_pool_exchange'] = np.sum(dffrac['MassFrac_{0}'.format(condition_id)]) * Ptot * sigma
    
    for go in dffrac.index:
        rxn_id = 'prot_pool_exchange_{0}'.format(go.replace(':',''))
        go_pools[rxn_id] = dffrac.loc[go,'MassFrac_{0}'.format(condition_id)] * Ptot * sigma
        
    return one_pool,go_pools


# In[8]:


def updateKcats(ecmodel,kcats):
    # kcats: dict() {(rxn_id, prot_met_id): kcat in h-1}
    for (rxn_id, prot_met_id),kcat in kcats.items():
        rxn = ecmodel.reactions.get_by_id(rxn_id)
        prot_met = ecmodel.metabolites.get_by_id(prot_met_id)
        
        new_coeff = -1./kcat
        
        ecpy.change_rxn_coeff(rxn,prot_met,new_coeff)


# In[9]:


def optimize_model(model,pools,tag=''):
    model.objective = 'Biomass_v1'
    model.objective_direction = 'max'
    s = model.optimize()
    r = s.objective_value if s.status == 'optimal' else -1
    print('Results {0}:'.format(tag))
    print('  Status        :',s.status)
    print('  growth_rate   :',r)
    print('  glucose uptake:',s.fluxes['Exchange_Glucopyranose'])
    print('  PHA           :',s.fluxes['PHA_secretion'])
    print('  NGAM:',model.reactions.NGAM.lower_bound)
    print('   PHA:',model.reactions.PHA_secretion.lower_bound)
    for pool_rxn in pools.keys():
        print('  Protein pool  :',model.reactions.get_by_id(pool_rxn).upper_bound)
    print()
            
    return r


# In[10]:


def simulate_growth(model,pools,dfpheno,kcats,condition_id,tag=''):
    # kcats: dict() {(rxn_id, prot_met_id): kcat in h-1}
    # return two growth rates: one for model with NGAM and PHA constraints and one for without
    
    with model:
        updateProteinConstraints(model,pools)
        updateKcats(model,kcats)
        
        PHA = dfpheno.loc[condition_id,'P3HB']
        if ~np.isnan(PHA): 
            model.reactions.PHA_secretion.lower_bound = PHA
        
        try: r = optimize_model(model,pools)    
        except: r = -1
    return [r]


# In[11]:


def test_model(#model_one_pool,      # json format
               model_go_pools,
               kcats,          # kcats: dict() {(rxn_id, prot_met_id): kcat in h-1}
               dftot,          # contains the total protein abandance, in the unit of gram protein/gram CDW
               dfpheno,        # contains the P3HB specific synthetic rate
               dffrac,
               condition_id,   # condition id
                ):
    
    
    # constrain protein pool
    print('Condition:',condition_id)
    Ptot = dftot.loc[condition_id,'Ptot']
    one_pool,go_pools = buildPoolDict(dffrac,condition_id,Ptot)
    
    #r1 = simulate_growth(model_one_pool,one_pool,dfpheno,kcats,condition_id,tag='one pool')
    r2 = simulate_growth(model_go_pools,go_pools,dfpheno,kcats,condition_id,tag='go pools')
    
    print(''.join(['\n']*2))
    return r2


# In[12]:


def simulate_a_single_kcat_set_on_all_datasets(datasets,kcats):
    start_time = time.time()
    results = {}
    #model_one_pool = pickle.load(open(datasets['model_one_pool_file'],'rb'))
    model_go_pools = pickle.load(open(datasets['model_go_pools_file'],'rb'))
    
    datasets_copy = datasets.copy()
    #datasets_copy.pop('model_one_pool_file')
    datasets_copy.pop('model_go_pools_file')
    
    for condition_id in datasets['dfpheno'].index:
        datasets_copy['condition_id'] = condition_id
        res = test_model(model_go_pools,kcats=kcats,**datasets_copy)
        results[condition_id] = res
    print('Time:',time.time() - start_time)
    return results


# In[13]:


class smc_abc:
    def __init__(self,datasets,priors,epsilon,outfile,cores=cpu_count(),population_size=100):
        self.datasets   = datasets
        self.population = []
        self.epsilons   = [np.inf] # store the distance after each generation
        self.prior      = priors.copy()
        self.posterior  = priors.copy() # a dictionary with rxn_id as key, RV as value
        self.cores      = cores
        
        self.outfile    = outfile
        self.population = []  # a list of populations [p1,p2...]
        self.distances  = []    # a list of distances for particles in population
        self.simulations = 0  # number of simulations performed 
        self.simulated_data = []
        
        self.population_size = population_size
        
        #self.all_simulated_particle = []
        #self.all_simulated_data     = []
        #self.all_simulated_distance = []
        self.epsilon = epsilon
        
    def simulate_one(self,particle,index,Q):
        '''
        particle:  parameters 
        Q:      a multiprocessing.Queue object
        index:  the index in particles list
        '''
        
        res = simulate_a_single_kcat_set_on_all_datasets(self.datasets,particle)

        Q.put((index,res))
    
    def distance(self,res):
        y_sim = []
        y_exp = []
        for condition_id,lst in res.items():
            y_sim.extend(lst)
            y_exp.extend([self.datasets['dfpheno'].loc[condition_id,'SpecificGrowthRate']])
            
        return mean_squared_error(y_exp,y_sim)
        
    def calculate_distances_parallel(self,particles):
        Q = Manager().Queue()
        jobs = [Process(target=self.simulate_one,args=(particle,index,Q)) 
                               for index,particle in enumerate(particles)]
        
        for p in jobs: p.start()
        for p in jobs: p.join()
        
        distances = [None for _ in range(len(particles))]
        simulated_data = [None for _ in range(len(particles))]

        for index,res in [Q.get(timeout=1) for p in jobs]: 
            distances[index] = self.distance(res)
            simulated_data[index] = res
        
        # save all simulated results
        #self.all_simulated_data.extend(simulated_data)
        #self.all_simulated_distance.extend(distances)
        #self.all_simulated_particle.extend(particles)
        
        del Q, jobs
        return distances,simulated_data
        
    def simulate_a_generation(self):
        particles_t, simulated_data_t, distances_t = [], [], []
        while len(particles_t) < self.population_size:
            self.simulations += self.cores
            particles = [{idp: 10**rv.sample() for idp,rv in self.posterior.items()} for i in range(self.cores)]

            distances,simulated_data = self.calculate_distances_parallel(particles)
            
            particles_t.extend(particles)
            simulated_data_t.extend(simulated_data)
            distances_t.extend(distances)
        
        return particles_t, simulated_data_t, distances_t
    
    def update_population(self,particles_t, simulated_data_t, distances_t):
        print ('updating population')
        # save first generation
        if len(self.population) == 0:
            self.population_t0 = particles_t.copy()
            self.distances_t0 = distances_t.copy()
            self.simulated_data_t0 = simulated_data_t.copy()
        
        
        combined_particles = np.array(self.population + particles_t)
        combined_distances = np.array(self.distances + distances_t)
        combined_simulated = np.array(self.simulated_data + simulated_data_t)
        
        sort_index = np.argsort(combined_distances)
        self.population = list(combined_particles[sort_index][:self.population_size])
        self.distances = list(combined_distances[sort_index][:self.population_size])
        self.simulated_data = list(combined_simulated[sort_index][:self.population_size])
        self.epsilons.append(np.max(self.distances))
        
        print('Model: epsilon=',str(self.epsilons[-1]))
        
    def update_posterior(self):
        print ('Updating prior')
        parameters = dict()   # {'':[]}
        for particle in self.population:
            for p,val in particle.items(): 
                lst = parameters.get(p,[])
                lst.append(np.log10(val))
                parameters[p] =lst
        
        for p, lst in parameters.items():
            self.posterior[p] = RV(loc = np.mean(lst), scale = np.std(lst))
        
        
    def run_simulation(self):
        while self.epsilons[-1] > self.epsilon:
            particles_t, simulated_data_t, distances_t = self.simulate_a_generation()
            self.update_population(particles_t, simulated_data_t, distances_t)
            self.update_posterior()
            pickle.dump(self,open(self.outfile,'wb'))
        


# In[ ]:


if __name__ == '__main__':
    #model_one_pool_file = '../Results/template_ecModel_one_pool.pkl'
    model_go_pools_file = '../Results/template_ecModel_goslim_pools.pkl'
    
    dffrac = pd.read_csv('../Results/protein_abundance_go_slim_level_uniq_asc.csv',index_col=0)
    dftot = pd.read_csv('../proteomics/total_protein_abandance_mean.csv',index_col=0)
    dfpheno = pd.read_csv('../proteomics/phynotype.csv',index_col=0,comment='#')
    df_enz_kcat = pd.read_csv('../Results/mapped_kcats_updated_with_ko.csv',index_col=0)
    
    datasets = {
        #'model_one_pool_file': model_one_pool_file,
        'model_go_pools_file': model_go_pools_file,
        'dffrac'      : dffrac,
        'dfpheno'     : dfpheno,
        'dftot'       : dftot,
    }
    
    step1_outfile = '../Results/smc_abc_light_go_and_one_asc_step_1.pkl'
    step1_res = pickle.load(open(step1_outfile,'rb'))
    priors = step1_res.posterior
    #tmp_ecmodel = pickle.load(open(model_one_pool_file,'rb'))
    #priors = build_priors(tmp_ecmodel,df_enz_kcat)
    
    
    epsilon = 0
    outfile = '../Results/smc_abc_light_go_and_one_asc_step_2.pkl'
    if os.path.isfile(outfile):
        experiment = pickle.load(open(outfile,'rb'))
    else:
        experiment = smc_abc(datasets,priors,epsilon,outfile)
    experiment.run_simulation()


# def run_test():
#     #model_one_pool_file = '../Results/template_ecModel_one_pool.pkl'
#     model_go_pools_file = '../Results/template_ecModel_goslim_pools.pkl'
#     
#     dffrac = pd.read_csv('../Results/protein_abundance_go_slim_level_uniq_asc.csv',index_col=0)
#     dftot = pd.read_csv('../proteomics/total_protein_abandance_mean.csv',index_col=0)
#     dfpheno = pd.read_csv('../proteomics/phynotype.csv',index_col=0,comment='#')
#     df_enz_kcat = pd.read_csv('../Results/mapped_kcats_updated_with_ko.csv',index_col=0)
#     
#     datasets = {
#         #'model_one_pool_file': model_one_pool_file,
#         'model_go_pools_file': model_go_pools_file,
#         'dffrac'      : dffrac,
#         'dfpheno'     : dfpheno,
#         'dftot'       : dftot,
#     }
#     
#     step1_outfile = '../Results/smc_abc_light_go_and_one_asc_step_1.pkl'
#     step1_res = pickle.load(open(step1_outfile,'rb'))
#     priors = step1_res.posterior
#     #tmp_ecmodel = pickle.load(open(model_one_pool_file,'rb'))
#     #priors = build_priors(tmp_ecmodel,df_enz_kcat)
#     
#     
#     epsilon = 0
#     outfile = '../Results/smc_abc_light_go_and_one_asc_step_2_test.pkl'
#     experiment = smc_abc(datasets,priors,epsilon,outfile,cores=4,population_size=4)
#     experiment.run_simulation()

# run_test()

# def load_kcats_test(tmp_ecmodel,df_enz_kcat):
#     # tmp_ecmodel, an ecModel with all enzymes and pools
#     # in df_enz_kcat, value are in the unit of 1/s, log1o-transformed
#     # in the resulting prior, kcat is in log10-transformed, 1/h
#     
#     initial_kcats = initialize_kcats_dict(tmp_ecmodel)
#     
#     # Build priors
#     kcats = dict()
#     for rxn_id,prot_met_id in initial_kcats.keys():
#         rxn_short_id = rxn_id.split('No')[0]
#         kcats[(rxn_id,prot_met_id)] = 10**df_enz_kcat.loc[rxn_short_id,'log10_kcat_mean']*3600
#     return kcats

# tmp_ecmodel.reactions.get_by_id('RXN-12565No1').reaction

# if __name__ == '__main__':
#     model_one_pool_file = '../Results/template_ecModel_one_pool.pkl'
#     #model_go_pools_file = '../Results/template_ecModel_goslim_pools.pkl'
#     
#     dffrac = pd.read_csv('../Results/protein_abundance_go_slim_level_uniq_asc.csv',index_col=0)
#     dftot = pd.read_csv('../proteomics/total_protein_abandance_mean.csv',index_col=0)
#     dfpheno = pd.read_csv('../proteomics/phynotype.csv',index_col=0,comment='#')
#     df_enz_kcat = pd.read_csv('../Results/mapped_kcats_updated_with_ko.csv',index_col=0)
#     
#     datasets = {
#         'model_one_pool_file': model_one_pool_file,
#         #'model_go_pools_file': model_go_pools_file,
#         'dfpheno'     : dfpheno,
#         'dftot'       : dftot,
#     }
#     
#     
#     tmp_ecmodel = pickle.load(open(model_one_pool_file,'rb'))
#     kcats=load_kcats_test(tmp_ecmodel,df_enz_kcat)
#     

# In[ ]:





# In[ ]:




