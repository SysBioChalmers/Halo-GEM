`halo-network.tsv` is a tab-sperated file that has two columns: Source and Target. It's a petri net representation of the metabolic network. The directed edges connect reactants to reactions and reactions to products. Thereby for a single connection listed in this file, there are two possible connections: "Reaction ID to Product ID" or "Reactant ID to Reaction ID".  
  
`node_type.tsv` is a tab-sperated file that tells you the type of IDs in `halo-network.tsv`: either "RXN" or "MET".
   
`rxns.tsv` is a tab-sperated file that has two columns: reaction ID and reaction name. 
  
  
`mets.tsv` is a tab-sperated file that has two columns: metabolite ID and metabolite name. 


