## Renconstruction of the first draft Halo-GEM: Halo_GEM_v1

#### Based on MetaCyc and KEGG
The RAVEN 2.0 toolbox were applied. For KEGG based reconstruction, an evalue cutoff of 1e-50 was used. For MetaCyc based reonstruction, an bit score cutoff of 100 was used. 
**Halo_meta**: `halo_metacycmodel_100.mat`
```
Number of reactions: 1271
Number of metabolits: 1520
Number of compartments: 1
Number of genes: 1230
Number of missing genes: 0
Number of reactions with missing genes: 0
```

**Halo_kegg**: `halo_keggmodel_50.mat`
```
Number of reactions: 1462
Number of metabolits: 1639
Number of compartments: 1
Number of genes: 852
Number of missing genes: 0
Number of reactions with missing genes: 0
```

Script for this step: `buildModelFromKEGG_Metacyc.m`

#### Based on template iML1515
The orgholog pairs between H TD01 and E.coli were identified with a bi-directional best hit (BBH) apporach. The cutoff for evalue is 1e-20 and the length coverage is 45%. With this approach, 718 BBH pairs were identified. Those genes are involved in 1235 reactions in iMM1515. After refining gene-reaction rules, there are 136 reactions contain extra 144 required genes that are failed  in mapping to H. TD01 proteome with BBH approach. Then the best homolog of those genes in H TD01 were identified (evalue<1e-20 and length coverage >45%). After this step, there are 116 reactions with extra 96 missing ensential genes to catalyze those reactions. Those reactions were then removed from the model. 

**Halo_eco:**`../Results/halo_iML1515_template_without_missing_genes.json`
```
Number of reactions: 1119
Number of metabolits: 1099
Number of compartments: 3
Number of genes: 662
Number of missing genes: 0
Number of reactions with missing genes: 0
```

Script for this step: `gap_filling-core_eco_biomass_add_template_rxns_without_missing_genes.ipynb`

#### Curate metabolite ids in the Metacyc based reconstruction
(1) We identified that some metabolites in Metacyc based reconstruction are in the form of -GLU-N(poly-glutamine). After carefully check the MetaCyc website and cross-referenced to KEGG database, N in those metabolites should be 0. Then we replaced those ones in -GLU-N format with non-GLU-N format. (2)  In the case that one metabolite has more than one metacyc id existing in the model, remove one of them. (3) Since there is no compartment information in the transport reactions, the metabolites transported are missing in the reaction because they exit in both side of the reaction. We manually checked the metacyc website for those transporters and then added them back (`../../../ComplementaryData/transport_rxns_mannual.tsv`). 
Script for this step: `gap_filling-core_eco_biomass_add_template_rxns_without_missing_genes.ipynb`

#### Combination of different reconstructions
Since Metacyc is a database containing reactions with experiemntal evidences, the reconstruction based on Metacyc is considered with higher quality than KEGG based. Thereby, only Metacyc based and iM1515 based reonstructions were used for combination.

1. Compartments  
Only cytosol and extracellular space were considered. For those metabolites in **Halo_eco**, consider all metabolites present in periplasm as in cytosol.   

2. Metabolite and Reaction IDs  
Use metacyc id for metabolites and reactions if available. 

3. Combine Gr rules.
If the rxn exist in both models, combine their gr rule as follows: If 'and' not in both gr rules. Combine all genes with 'or'. In all other cases, use the one from iMM1515 based reonstruction.

4. Add exchange reactions for all metabolites present in the medium, as defined in `medium` sheet of `../Results/core_biomass_medium_summary_iML1515.xlsx`

5. Add exchange reactions for all secretion metabolites, as defined in `secretion` sheet of `../Results/core_biomass_medium_summary_iML1515.xlsx`  
Now the model contains
```
Number of reactions: 1889
Number of metabolits: 2005
Number of compartments: 2
Number of genes: 1270
Number of missing genes: 0
Number of reactions with missing genes: 0
```

6. Add protein, DNA, RNA, ions, lipid synthsis reactions by using the coefficients in the core biomass equation of iMM1515 model.


#### Gap-filling with iMM1515 as a universal model
Since the H. TD01 can grow in the medium without amino acids, all uptake reactions for amino acids were blocked before gap-filling. Then step by step gap-filling were performed to make the model be able to produce the components in the biomass equation: DNA, RNA, Protein, lipids and others. As defined in `../Results/core_biomass_medium_summary_iML1515.xlsx`.

35 reactions were added by gap-filling. The resulting model contains 
```
Number of reactions: 1929
Number of metabolits: 2014
Number of compartments: 2
Number of genes: 1298
Number of missing genes: 28
Number of reactions with missing genes: 26
```
Those genes that were added together with gap-filled reactions were mapped back to the H. TD01 proteome. With a evalue cutoff of 1e-20 and lenght coverage of 45%, 25/28 genes have a homolog in H TD01 proteome. The resulting model contains
```
Number of reactions: 1929
Number of metabolits: 2014
Number of compartments: 2
Number of genes: 1273
Number of missing genes: 3
Number of reactions with missing genes: 1
```

For those reactions still with missing genes, remove their gr rules and put a note as ```{'source':'gapfilled from iJO1366, genes were not found in H.TD01'}```  

Script for this step: `gap_filling-core_eco_biomass_add_template_rxns_without_missing_genes.ipynb`
#### Add Biomass equation
```
0.000223 10-FORMYL-THF_c + 0.000223 2-OCTAPRENYL-6-HYDROXYPHENOL_c + 2.6e-05 2fe2s_c + 75.37723 ATP_c + 2e-06 BIOTIN_c + 0.000576 CO-A_c + 0.00026 CPD-7_c + 5.5e-05 CPD-9649_c + 0.013894 CPD0-2278_c + 0.000223 FAD_c + 0.000223 METHYLENE-THF_c + 0.000447 NADP_c + 0.001831 NAD_c + 0.000223 PROTOHEME_c + 0.000223 PYRIDOXAL_PHOSPHATE_c + 0.000223 RIBOFLAVIN_c + 0.000223 S-ADENOSYLMETHIONINE_c + 0.000223 SIROHEME_c + 9.8e-05 SUC-COA_c + 0.000223 THF_c + 0.000223 THIAMINE-PYROPHOSPHATE_c + 70.028756 WATER_c + dna_c + ions_c + lipids_c + protein_c + rna_c --> 75.37723 ADP_c + 0.773903 PPI_c + 75.37723 PROTON_c + 75.37323 Pi_c
```
Save the resulted model in different formats, including `.mat`, `.xml` and `.json`
```
Halo-GEM/ModelFiles/mat/combine_halo_meta100_iML1515_gap_fill_with_core_biomass_without_missing_genes.mat
Halo-GEM/ModelFiles/xml/combine_halo_meta100_iML1515_gap_fill_with_core_biomass_without_missing_genes.xml
Halo-GEM/ModelFiles/json/combine_halo_meta100_iML1515_gap_fill_with_core_biomass_without_missing_genes.json
```
Script for this step: `gap_filling-core_eco_biomass_add_template_rxns_without_missing_genes.ipynb`

#### Rename model
Rename those models as 
```
Halo-GEM/ModelFiles/mat/Halo_GEM_v1.mat
Halo-GEM/ModelFiles/xml/Halo_GEM_v1.xml
Halo-GEM/ModelFiles/json/Halo_GEM_v1.json
```



















