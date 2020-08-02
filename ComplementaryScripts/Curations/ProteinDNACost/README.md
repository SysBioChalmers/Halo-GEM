### The simulation for cost esimation of single amino acid or dNTP:  
(1) block the free uptake of 20 amino acids  
(2) set the glucose uptake as 1 mmol/gDW/h  
(3) add reaction "--> met", met is one of amino acids or dNTP   
(4) maximize the flux in (3), the objective value is the maximal amount of aa or dNTP produced from 1 mmol/gDW/h glucose  
(5) the cost is defined as the amount of glucose required to produce 1 mmol/gDW/h aa or dNTP. The reverse of the value from (4)   

Two files were obtained: `aa_cost.csv` and `dNTP_cost.csv`.  

### Protein cost
Protein and DNA cost: weighted summation of values obtained in the above, based on the amino acid composition or dNTP composition.
  
Result file: `protein_cost.csv`