### Pipeline

* Step 1: Match kcat values to enzymes and reactions: only those enzymes with EC number. Not transport reaction

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

### Match kcat values

### TO DO
MW should be in the unit of g/mol
kcat should be in 1/h

* Matach kcat values.

* Check free upper bound to +inf
* Try to test it on convert Halo-GEM