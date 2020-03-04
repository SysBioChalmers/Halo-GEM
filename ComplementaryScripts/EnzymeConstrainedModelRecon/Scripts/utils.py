from cobra import Model, Reaction, Metabolite

def parse_gr_rule(gr):
    '''
    Parse gr rule into a list of components. 
    gr: gene reaction rule, defined in cobrapy.
    
    For example: 
    
    Input         : Output
    A or B        : ["A", "B"]
    (A and B)     : ["A and B"]
    (A and B) or C: ["A and B","C"]
    
    Usage: complexes = parse_gr_rule(gr)
    
    Gang Li, last updated 2020-03-04
    
    '''
    complexes = [item.strip().replace('(','').replace(')','') for item in gr.split('or')]
    if len(complexes) < 2 and len(complexes[0]) < 1: complexes = []
    
    return complexes

def constrainPool(model,MWs, non_measured,UB):
    '''
    Adapted from gekecomat, consstrainPool.m
    model       : eModel from convertToEnzymeModel()
    MWs         : a dictionary with molecular weight of enzymes
    non_measured: a list of enzymes without proteomics data
    UB          : upper bound for the combined pool of those non_measured enzymes
    
    Define new rxns: For each enzyme, add a new rxn that draws enzyme from the
    enzyme pool (a new metabolite), and remove previous exchange rxn. The new
    rxns have the following stoichiometry (T is the enzyme pool):
     MW[i]*P[T] -> P[i]
     
    Usage: model = constrainPool(model,MWs, non_measured,UB)
    
    Gang Li, last updated 2020-03-04
    '''
    model = model.copy()
    # create prot_pool metabolite 
    prot_pool = Metabolite('prot_pool')
    prot_pool.name = prot_pool.id
    
    rxns_to_add  = list()
    rxns_to_drop = list()
    for prot in non_measured:
        prot_exchange_rxn = model.reactions.get_by_id('prot_{0}_exchange'.format(prot))
        
        draw_rxn = Reaction('draw_prot_{0}'.format(prot))
        draw_rxn.name = draw_rxn.id
        draw_rxn.add_metabolites({prot_pool:-MWs[prot],list(prot_exchange_rxn.metabolites)[0]:1})
        
        rxns_to_add.append(draw_rxn)
        rxns_to_drop.append(prot_exchange_rxn)
        
    # add draw reaction into model
    model.add_reactions(rxns_to_add)
    model.remove_reactions(rxns_to_drop)
    
    # add prot_pool_exchange rxn
    rxn_prot_pool_exg = Reaction('prot_pool_exchange')
    rxn_prot_pool_exg.name = rxn_prot_pool_exg.id
    rxn_prot_pool_exg.add_metabolites({prot_pool:1})
    rxn_prot_pool_exg.lower_bound = 0
    rxn_prot_pool_exg.upper_bound = UB
    
    model.add_reaction(rxn_prot_pool_exg)
    
    return model
        
    
def convertToEnzymeModel(model,kcats):
    '''
    model .   : irrevModel
    kcats     : a dictionary with kcat values {('protein_id',rxn_id):100,...}
    
    Usage: eModel = convertToEnzymeModel(model,kcats)
    
    Gang Li, last updated 2020-03-04
    '''
    converted_reaction_list = []
    protein_exchange_rxns = {}
    for rxn in model.reactions:
        complexes = parse_gr_rule(rxn.gene_reaction_rule)
        
        # 1. for those reactions without genes 
        if len(complexes) <1: 
            converted_reaction_list.append(rxn)
            continue
        
        # 2. for those reactions with genes, but no kcat
        first_gene = [gene.id for gene in rxn.genes][0]
        if kcats.get((first_gene,rxn.id),None) is None:
            converted_reaction_list.append(rxn)
            continue
        
        # 3. for those reactions with isoenzymes, add arm reaction
        if len(complexes) >1:
            rxn_new, arm_rxn = getArmReaction(rxn)
            converted_reaction_list.append(arm_rxn)
            
            for i,complx in enumerate(complexes):
                prots = [item.strip() for item in complx.split('and')]
                kcat = kcats[(prots[0],rxn.id)]
                e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn_new, kcat, complx,rxn_index=i+1)
                
                converted_reaction_list.append(e_rxn)
                for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rxn.id] = prot_exchange_rxn
                
            continue
         
        if len(complexes) == 1:
            complx = complexes[0]
            prots = [item.strip() for item in complx.split('and')]
            kcat = kcats[(prots[0],rxn.id)]
            e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn, kcat, complx)
            converted_reaction_list.append(e_rxn)
            for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rxn.id] = prot_exchange_rxn
    
    eModel = Model()
    eModel.add_reactions(converted_reaction_list)
    eModel.add_reactions(protein_exchange_rxns.values())
    
    return eModel


def addEnzymesToRxn(rxn, kcat, protIDs, rxn_index=None):
    '''
    Add each enzyme as one of metabolites in the model, Current version does not support stoichiometric cofficients of subunits
    rxn      : the input Reaction object in cobrapy
    rxn_index: a integer like 1, for isoenzymes
    kcats    : kcat value for the reaction
    protIDs  : a single protein name, like "A", or a complex like "A and B". String
    MWs      : a dictionary with prot_id as key and molecular weight as value
    
    Usage: e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn, kcat, protIDs,MWs)
    
    Gang Li, last updated 2020-03-03
    '''
    
    e_rxn      = rxn.copy()
    if rxn_index is not None:
        e_rxn.id   = e_rxn.id + 'No{0}'.format(rxn_index)
        e_rxn.name = e_rxn.name + ' (No{0})'.format(rxn_index)
    prots = [item.strip() for item in protIDs.split('and')]
    
    
    # get compartment
    comp = None
    for met in rxn.metabolites:
        comp = met.compartment
        if rxn.get_coefficient(met)<0: comp = met.compartment
    
    # create Metabolite object for each protein and create exchange reaction
    prot_mets = []
    prot_exchange_rxns = []
    for prot in prots:
        prot_met = Metabolite('prot_{0}_{1}'.format(prot,comp))
        prot_met.compartment =  comp
        prot_mets.append(prot_met)
        
        # add excange reaction of protein
        excg_rxn = Reaction('prot_{0}_exchange'.format(prot))
        excg_rxn.lower_bound = 0
        excg_rxn.gene_reaction_rule = prot
        excg_rxn.add_metabolites({prot_met:1})
        prot_exchange_rxns.append(excg_rxn)
        
    # add enzymes into reaction
    e_rxn.add_metabolites({prot_met:-1./kcat for prot_met in prot_mets})
    e_rxn.gene_reaction_rule = protIDs

    return e_rxn, prot_exchange_rxns
    
    

def convertToIrrev(model,rxns=None):
    '''
    Adapted from convertToIrrev.m in RAVEN. Split all reversible reactions into one forward and one reverse reaction.
    model: the input model
    rxns : a list of reactions to be converted. default is model.reactions
    
    Usage: irrevModel = convertToIrrev(model,rxns)
    
    Gang Li, last updated 2020-03-03
    '''
    if rxns is None: rxns = model.reactions
    
    converted_reaction_list = []
    for rxn in rxns:
        
        # irreversible reactions
        if rxn.lower_bound>=0: converted_reaction_list.append(rxn) 
        
        # reversible_reacions
        else:
            rxn_REV      = rxn.copy()
            rxn_REV.id   = rxn.id + '_REV'
            rxn_REV.name = rxn.name + ' (reversible)'
            
            rxn_REV.add_metabolites({met:-rxn.get_coefficient(met)*2 for met in rxn.metabolites})
            rxn_REV.lower_bound = 0
            rxn_REV.upper_bound = -rxn.lower_bound
            
            rxn.lower_bound = 0
            converted_reaction_list.extend([rxn,rxn_REV])
    
    
    # build irrevModel
    irrevModel = Model()
    irrevModel.add_reactions(converted_reaction_list)
    
    return irrevModel
            

def getArmReaction(rxn):
    '''
    Adapted from addArmReaction.m from geckomat. Add an arm reaction for the selected reaction in the model.
    
    rxn: the reaction Object in cobrapy
    
    Original reaction: A + B --> C + D
    
    Arm reaction    : A + B --> pmet   (no gr rule)
    Change the orginial reaction to:  pmet --> C + D (use old gr rule)
    
    The arm reaction has a id format of "arm_rxnID" and a name format of "rxnName (arm)"
    
    The intermediate metabilite has a name format of "pmet_rxnID"
    
    The arm reaction shares the same lb, ub, gr rules, subsystems with original reaction.
    
    Compartment: fistly try to use the same compartment as substrates, then products', otherwise None.
    
    
    Usage: rxn_new, arm_rxn = addArmReaction(model,rxn_id).
    
    Gang Li, Last update: 2020-03-03
    '''
    
    # 1. create intermediate metabilite
    rxnID = rxn.id
    comp = None
    for met in rxn.metabolites:
        comp = met.compartment
        if rxn.get_coefficient(met)<0: comp = met.compartment
    
    pmet = Metabolite('pmet_{0}'.format(rxnID),compartment=comp)
    
    # 2. create arm reaction: 
    arm_rxn                    = Reaction('arm_{0}'.format(rxnID))
    arm_rxn.name               = rxn.name + ' (arm)'
    arm_rxn.subsystem          = rxn.subsystem
    arm_rxn.lower_bound        = rxn.lower_bound  
    arm_rxn.upper_bound        = rxn.upper_bound  
    arm_rxn.gene_reaction_rule = ''
    
    mets = {met:rxn.get_coefficient(met) for met in rxn.metabolites if rxn.get_coefficient(met)<0}
    mets[pmet] = 1
    
    arm_rxn.add_metabolites(mets)
    
    # 3. change orignal reaction to pmet --> C + D 
    rxn_new = rxn.copy()
    rxn_new.subtract_metabolites({met:rxn_new.get_coefficient(met) for met in rxn_new.metabolites if rxn_new.get_coefficient(met)<0})
    rxn_new.add_metabolites({pmet:-1})
    
    return rxn_new, arm_rxn

    