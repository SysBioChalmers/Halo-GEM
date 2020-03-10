### Pipeline

* Step 1: Match kcat values to enzymes and reactions: only those enzymes with EC number. This step is done by 
```
match_kcat_values_and_MWs.ipynb
```
This script outputs two csv files: 
```
Results/mapped_kcats.csv
Results/enzyme_MWs.csv 
```
MW isin the unit of g/mmo(kDa).  
kcats and associated uncertainties were given by log10-transfromed 1/s, when integrating into the model, it will be converted to 1/h.

* Step 2: Convert the model  to irreversible model 
```
irrModel = utils.convertToIrrev(model)
```

* Step 3: Convert the irreversible model to enzyme model according kcats
```
eModel = utils.convertToEnzymeModel(irrModel,kcats)
```

* Step 4: Constrain protein pool
```
ecModel = utils.constrainPool(eModel,MWs, ['E1','E2','E4'],10)
```

An small scale toy model was provided in `test_with_a_toy_model.ipynb`


### TO DO
* Check free upper bound to +inf
* Specify protein pool.
* Integration of proteomics data