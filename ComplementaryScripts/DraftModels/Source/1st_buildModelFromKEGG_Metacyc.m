%% BUild the draft model with RAVEN2 tool box.
% Try different blast settings to generate several models.
% 
% Gang Li
% 2018-10-28
clc
clear

%% 1. get model from KEGG
dataDir = 'prok90_kegg87';
fastaFile = '../../../ComplementaryData/protein_sequences_v2.0.fasta';
% scores = [-50];
% for i = 1:length(scores)
%     score = 10^scores(i);
%     keggmodel = getKEGGModelForOrganism('halo',fastaFile,dataDir,'',1,0,0,0,score);
%     outname = strcat('../Results/halo_keggmodel_',num2str(-scores(i)),'.mat');
%     save(outname,'keggmodel');
% end
%% 2. get model from Metacyc
scores = [100];
for i = 1:length(scores)
    metacycmodel=getMetaCycModelForOrganism('halo',fastaFile,1,0,0,scores(i));
    outname = strcat('../Results/halo_metacycmodel_',num2str(scores(i)),'.mat')
    save(outname,'metacycmodel');
end

%% 3. combine models
clear
load('../Results/halo_metacycmodel_100.mat')
load('../Results/halo_keggmodel_50.mat')
model = combineMetaCycKEGGModels(metacycmodel,keggmodel);
save('../Results/halo_metacyc_kegg.mat','model')