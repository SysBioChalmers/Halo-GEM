# Halomonas-TD01-GEM

- Brief Repository Description

This repository contains a genome-scale metabolic model **Halo** for a salt-tolerant bacteria _Halomonas bluephagenesis_ TD01.

- Abstract:
The halophile _Halomonas bluephagenesis_ TD01 and its derivatives have been successfully developed as a low-cost platform for the unsterile and continuous production of chemicals. Here we present the de novo reconstruction of its GEM. Two draft models from E. coli GEM iJO1366 (by cobrapy) and MetaCyc (by RAVEN 2.0) were firstly generated and then combined. The core-biomass equation from iJO1366
was used. Gap-filling was performed with iJO1366 as universal model.

- Model KeyWords:

**GEM Category:** Species; **Utilisation:** Predictive simulation; **Field:** Metabolic-network reconstruction; **Type of Model:** Reconstruction; **Taxonomy:** _Halomonas bluephagenesis_; **Metabolic System:** General Metabolism; **Strain:** TD01; **Condition:** Defined medium;

- Reference:
> None

- Pubmed ID: TBA

- Last update: 2019-03-21

- The model contains:

| Taxonomy | Template Model | Reactions | Metabolites| Genes |
| ------------- |:-------------:|:-------------:|:-------------:|-----:|
| _Halomonas bluephagenesis_ TD01 | None | 1,909 | 1,989 | 1,282 |

This repository is administered by Gang Li ([@SysBioChalmers](https://github.com/SysBioChalmers)), Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology

## Installation

### Recommended Software:
* A functional Matlab installation (MATLAB 7.3 or higher).
* [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) for MATLAB (required for contributing to development).
* libSBML MATLAB API ([version 5.16.0](https://sourceforge.net/projects/sbml/files/libsbml/5.13.0/stable/MATLAB%20interface/)  is recommended).
* [Gurobi Optimizer for MATLAB](http://www.gurobi.com/registration/download-reg).
* For contributing to development: a [git wrapper](https://github.com/manur/MATLAB-git) added to the search path.
* cobrapy


### Contribute To Development
1. Fork the repository to your own Github account
2. Create a new branch from [`devel`](https://github.com/SysBioChalmers/Streptomyces_coelicolor-GEM/tree/devel).
3. Make changes to the model
    + [RAVEN Toolbox 2](https://github.com/SysBioChalmers/RAVEN) for MATLAB is highly recommended for making changes
    + Before each commit, run in Matlab the `newCommit(model)` function from the `ComplementaryScripts` folder
    + Make a Pull Request to the `devel` folder, including changed `txt`, `yml` and `xml` files

## Contributors
* [Gang Li](https://www.chalmers.se/en/staff/Pages/gangl.aspx) ([@Gangl2016](https://github.com/Gangl2016)), Chalmers University of Technology, Göteborg, Sweden
* [Hao Wang](https://www.chalmers.se/en/staff/Pages/hao-wang.aspx) ([@Hao-Chalmers](https://github.com/Hao-Chalmers)), Chalmers University of Technology, Göteborg, Sweden
