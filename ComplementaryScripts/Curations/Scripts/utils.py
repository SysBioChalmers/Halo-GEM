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