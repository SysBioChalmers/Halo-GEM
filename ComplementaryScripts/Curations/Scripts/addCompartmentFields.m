%
% FILE NAME:    addCompartmentFields.m
%
% PURPOSE: This script is to modify the Halo-GEM by adding compartment-
%          associated fields, including 'comps', 'compNames' and 'metComps'
%


%% Load model
load('../../../ModelFiles/mat/Halo_GEM_v1.mat');
model = model_v1;


%% add compartment related fields to model

% There are 3 compartments (cytosol, peroxisome and extracellular) are
% currently modeled in Halo-GEM. Here is to get the index of metabolites
% that are from the corresponding compartments.
indC = find(endsWith(model.mets,'_c'));  % cytosol
indP = find(endsWith(model.mets,'_p'));  % peroxisome
indE = find(endsWith(model.mets,'_e'));  % extracellular

% add comps and compNames fields
model.comps = {'c';'p';'e'};
model.compNames = {'Cytosol';'Peroxisome';'Extracellular'};

% add metComps field
model.metComps = model.b;
model.metComps(indC) = 1;
model.metComps(indP) = 2;
model.metComps(indE) = 3;


%% Save updated model file

% write to SBML
exportModel(model, '../../../ModelFiles/xml/Halo_GEM_v1.xml');

% save to Matlab format
model_v1 = model;
save('../../../ModelFiles/mat/Halo_GEM_v1.mat','model_v1');


