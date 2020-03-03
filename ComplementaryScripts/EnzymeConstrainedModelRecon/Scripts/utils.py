from cobra import Model, Reaction, Metabolite


def addEnzymesToRxn(rxn, kcat, protIDs,MWs):
    '''
    Add each enzyme as one of metabolites in the model
    rxn    : the input Reaction object in cobrapy
    kcats  : kcat value for the reaction
    protIDs: a single protein name, like "A", or a complex like "A and B". String
    MWs    : a dictionary with prot_id as key and molecular weight as value
    
    Usage: e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn, kcat, protIDs,MWs)
    
    Gang Li, last updated 2020-03-03
    '''
    
    e_rxn = rxn.copy()
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
        prot_met.formula_weight = MWs[prot]
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
    irrevModel.objective = model.objective
    
    return irrevModel
            

def getArmReaction(model,rxnID):
    '''
    Adapted from addArmReaction.m from geckomat. Add an arm reaction for the selected reaction in the model.
    
    model: input cobra model
    rxnID: the reaction id
    
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
    rxn = model.reactions.get_by_id(rxnID)
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
    rxn_new.subtract_metabolites([met for met in rxn_new.metabolites if rxn_new.get_coefficient(met)<0])
    rxn_new.add_metabolites({pmet:-1})
    
    return rxn_new, arm_rxn

    