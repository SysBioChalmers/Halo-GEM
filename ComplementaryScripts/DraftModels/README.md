## Renconstruction of the first draft Halo-GEM: Halo_GEM_v1

Run those script under `Source/` sequentially:
* `1st_buildModelFromKEGG_Metacyc.m`
* `2nd_buildModelFromTemplate.ipynb`
* `3rd_curate_metacyc_ids.ipynb`
* `4th_add_compartments_and_update_transport_rxns.ipynb`
* `5th_combine_models_from_ecoli_metacyc_kegg.ipynb`
* `6th_add_medium_biomass.ipynb`
* `7th_gap_filling.ipynb`
* `8th_add_PHA_NGAM.ipynb`
* `9th_Simulate_Growth.ipynb`


### 1. Draft models with RAVEN 2.0 toolbox
**Script:** `1st_buildModelFromKEGG_Metacyc.m`
The RAVEN 2.0 toolbox were applied. For KEGG based reconstruction, an evalue cutoff of 1e-50 was used. For MetaCyc based reonstruction, an bit score cutoff of 100 was used.   
In the end, only the metacyc based reconstruction was used. 

**Output:**`Results/halo_metacycmodel_100.mat`
```
Number of reactions: 1321
Number of metabolits: 1553
Number of compartments: 1 {'s': ''}
Number of genes: 1252
Number of missing genes: 0
Number of reactions with missing genes: 0
```


### 2. Draft model Based on template iML1515 with a BBH search approach
**Script**: `2nd_buildModelFromTemplate.ipynb`
The orgholog pairs between H TD01 and E.coli were identified with a bi-directional best hit (BBH) apporach. The cutoff for evalue is 1e-3 and the length coverage is 45%. 

**Output:**`Results/halo_iML1515_template_without_missing_genes.json`
```
Number of reactions: 1215
Number of metabolits: 1173
Number of compartments: 3
Number of genes: 751
Number of missing genes: 0
Number of reactions with missing genes: 0
```


### 3. Curate metabolite ids in the Metacyc based reconstruction
**Script:** `3rd_curate_metacyc_ids.ipynb`
(1) We identified that some metabolites in Metacyc based reconstruction are in the form of -GLU-N(poly-glutamine). After carefully check the MetaCyc website and cross-referenced to KEGG database, N in those metabolites should be 0. Then we replaced those ones in -GLU-N format with non-GLU-N format. 
(2) In the case that one metabolite has more than one metacyc id existing in the model, remove one of them.  
(3) Split the reactions with `NAD-P-OR-NOP` or `NADH-P-OR-NOP` into two reactions, for instance 
```
HYDROXYPYRUVATE-REDUCTASE-RXN
GLYCERATE + NAD-P-OR-NOP <=> NADH-P-OR-NOP + OH-PYR + PROTON

HYDROXYPYRUVATE-REDUCTASE-RXN_NAD
GLYCERATE + NAD <=> NADH + OH-PYR + PROTON

HYDROXYPYRUVATE-REDUCTASE-RXN_NADP
GLYCERATE + NADP <=> NADPH + OH-PYR + PROTON
```
Manualy curations on the reversibility and ec numbers of those reactions were performed based on MetaCyc website.
**Output:**`Results/halo_metacycmodel_100_curated.pkl`


### 4. Assign compartments to metabolites in MetaCyc based reconstruction. 
**Script:** `4th_add_compartments_and_update_transport_rxns.ipynb`
For all non-transport reactions, consider all of them occurs in the cytosal.  
For transport reactions, since there is no compartment information in the transport reactions, the metabolites transported are missing in the reaction because they exit in both side of the reaction. We manually checked the metacyc website for those transporters and then added them back together with compartment information.  
**Output:**`Results/halo_metacycmodel_100_curated_compart.pkl`
```
Number of reactions: 1321
Number of metabolits: 1679
Number of compartments: 3 {'c': '', 'p': '', 'e': ''}
Number of genes: 1252
Number of missing genes: 0
Number of reactions with missing genes: 0
```

### 5.Combination of different reconstructions
**Script:** `5th_combine_models_from_ecoli_metacyc_kegg.ipynb`
At this step, reconstructions from MetaCyc and iML1515 were considered.  
Step 1: update the reaction and metabolite ids with priority metacyc_id > kegg_id > bigg_id
Step 2: With iML1515 based reconstruction, add new reactions from MetaCyc based reconstruction. If a reaction present in both reconstructions and if two gene reaction rules are only "or" relationship, then combine all the genes. Otherwise use the gene reaction rule from iML1515 based reconstruction.
**Output:**`Results/Results/halo_metacycmodel_100_curated_compart_with_eco_without_missing.pkl`
```
Number of reactions: 2139
Number of metabolits: 2279
Number of compartments: 3 {'c': '', 'p': '', 'e': ''}
Number of genes: 1345
Number of missing genes: 0
Number of reactions with missing genes: 0
```

### 6. Add medium uptake and biomass reaction
**Script:**`6th_add_medium_biomass.ipynb`
* Add exchange reactions for all metabolites present in the medium, as defined in `medium` sheet of `Results/core_biomass_medium_summary_iML1515.xlsx`

* Add exchange reactions for all secretion metabolites, as defined in `secretion` sheet of `Results/core_biomass_medium_summary_iML1515.xlsx`  
* Add general reactions for synthesis of protein, DNA, RNA, lips and ions for biomass. The core-biomass equation in IML1515 were used at this step.

```
0.295792 ARG_c + 0.241055 ASN_c + 0.09158 CYS_c + 0.26316 GLN_c + 0.26316 GLT_c + 0.612638 GLY_c + 0.094738 HIS_c + 0.290529 ILE_c + 0.513689 L-ALPHA-ALANINE_c + 0.241055 L-ASPARTATE_c + 0.450531 LEU_c + 0.343161 LYS_c + 0.153686 MET_c + 0.185265 PHE_c + 0.221055 PRO_c + 0.215792 SER_c + 0.253687 THR_c + 0.056843 TRP_c + 0.137896 TYR_c + 0.423162 VAL_c --> protein_c

0.026166 DATP_c + 0.027017 DCTP_c + 0.027017 DGTP_c + 0.026166 TTP_c --> dna_c

0.175 ATP_c + 0.133508 CTP_c + 0.215096 GTP_c + 0.144104 UTP_c --> rna_c

0.013013 AMMONIA_c + 0.005205 CA+2_c + 0.005205 CL-_c + 2.5e-05 CO+2_c + 7e-06 CPD-3_c + 0.000709 CU+2_c + 0.006715 FE+2_c + 0.007808 FE+3_c + 0.004338 HSO4_c + 0.195193 K+_c + 0.008675 MG+2_c + 0.000691 MN+2_c + 0.000323 NI+2_c + 0.000341 ZN+2_c --> ions_c

0.063814 CPD-12819_p + 0.075214 CPD-17086_p + 0.019456 KDO2-LIPID-IVA_e --> lipids_c
```

* Add biomass equation
```
0.000223 10-FORMYL-THF_c + 0.000223 2-OCTAPRENYL-6-HYDROXYPHENOL_c + 2.6e-05 2fe2s_c + 75.37723 ATP_c + 2e-06 BIOTIN_c + 0.000576 CO-A_c + 70.028756 CPD-15815_c + 0.00026 CPD-7_c + 5.5e-05 CPD-9649_c + 0.013894 CPD0-2278_p + 0.000223 FAD_c + 0.000223 METHYLENE-THF_c + 0.000447 NADP_c + 0.001831 NAD_c + 0.000223 PROTOHEME_c + 0.000223 PYRIDOXAL_PHOSPHATE_c + 0.000223 RIBOFLAVIN_c + 0.000223 S-ADENOSYLMETHIONINE_c + 0.000223 SIROHEME_c + 9.8e-05 SUC-COA_c + 0.000223 THF_c + 0.000223 THIAMINE-PYROPHOSPHATE_c + dna_c + ions_c + lipids_c + protein_c + rna_c --> 75.37723 ADP_c + 75.37323 CPD-16459_c + 0.773903 PPI_c + 75.37723 PROTON_c
```
**Output:**`Results/halo_metacycmodel_100_curated_compart_with_eco_without_missing_medium_biomass.pkl`
```
Number of reactions: 2236
Number of metabolits: 2337
Number of compartments: 3 {'c': '', 'p': '', 'e': ''}
Number of genes: 1345
Number of missing genes: 0
Number of reactions with missing genes: 0
```


### 7. Gap-filling with iML1515 as a universal model
**Script:**`7th_gap_filling.ipynb`
Since the H. TD01 can grow in the medium without amino acids, all uptake reactions for amino acids were blocked before gap-filling. Then step by step gap-filling were performed to make the model be able to produce the each of components in the biomass equation: DNA, RNA, Protein, lipids and others.  

Those genes that were added together with gap-filled reactions were mapped back to the H. TD01 proteome. With a evalue cutoff of 1e-3 and lenght coverage of 45%. 
**Output:**`esults/halo_metacycmodel_100_curated_compart_with_eco_without_missing_medium_biomass_gapfilled.pkl'
```
Number of reactions: 2266
Number of metabolits: 2338
Number of compartments: 3 {'c': '', 'p': '', 'e': ''}
Number of genes: 1356
Number of missing genes: 5
Number of reactions with missing genes: 3
```


Script for this step: `Source/reconstruction_pipeline.ipynb`
### 8. Add PHA synthesis reaction and Non-growth Associated ATP Maintanance.
**Script:**`8th_add_PHA_NGAM.ipynb`
```
CPHA_synthetase: PD-650_c --> CO-A_c + PHA_c
PHA_secretion: PHA_c --> 
NGAM: ATP_c + CPD-15815_c --> ADP_c + PROTON_c + Pi_c, lower_bound = 6.86
```
**Output:**`../Results/halo_metacycmodel_100_curated_compart_with_eco_without_missing_medium_biomass_gapfilled_PHA_NGAM.pkl`
```
Number of reactions: 2269
Number of metabolits: 2339
Number of compartments: 3 {'c': '', 'p': '', 'e': ''}
Number of genes: 1358
Number of missing genes: 5
Number of reactions with missing genes: 3
```


### 9. Some simulations to test the effect of glucose uptake on PHA, ATP, specific growth rate
**Script:**`9th_Simulate_Growth.ipynb`
There are 5 genes involving in 3 reactions wihout any homologs in H. TD01. 
```
SULFITE-REDUCT-RXN
H2SO3_c + 3.0 NADPH_c + 5.0 PROTON_c --> 3.0 CPD-15815_c + CPD-7046_c + 3.0 NADP_c
old: TD01GL001837 and b2764_missing
new: TD01GL001837

K2L4Aabctex
ATP_c + CPD-15815_c + KDO2-LIPID-IVA_p --> ADP_c + CPD-16459_c + KDO2-LIPID-IVA_e + PROTON_c
old: TD01GL000592 and TD01GL000593 and TD01GL002340 and TD01GL002341 and TD01GL002916 and b0641_missing and b3199_missing
new: TD01GL000592 and TD01GL000593 and TD01GL002340 and TD01GL002341 and TD01GL002916

THZPSN3
ATP_c + CPD-12279_c + DEOXYXYLULOSE-5P_c + NADPH_c + PROTON_c + iscssh_c --> AMP_c + CARBON-DIOXIDE_c + 2.0 CPD-15815_c + NADP_c + PPI_c + THZ-P_c + iscs_c
old: TD01GL000663 and TD01GL001192 and TD01GL003669 and b0423_missing and b3990_missing
new: TD01GL000663 and TD01GL001192 and TD01GL003669
```

For those reactions still with missing genes, remove missing genes and add a note in 
```python
rxn.notes['old_gr']= old_gr
```
The model contains
```
Number of reactions: 2269
Number of metabolits: 2339
Number of compartments: 3 {'c': '', 'p': '', 'e': ''}
Number of genes: 1353
Number of missing genes: 0
Number of reactions with missing genes: 0
```
 
Save the resulted model as 
```
Halo-GEM/ModelFiles/mat/Halo_GEM_v1.mat
Halo-GEM/ModelFiles/xml/Halo_GEM_v1.xml
Halo-GEM/ModelFiles/json/Halo_GEM_v1.json
```
