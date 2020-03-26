import cobra
import pandas as pd
import xlrd
import time
import os
from Bio import SeqIO
from scipy.io import loadmat
import re
import matplotlib.pyplot as plt
import numpy as np
import urllib
import pickle

def load_ecoli_model(mat_file,sbml_file):
    # iJO1366: load mat file: contains subsystem information. while sbml doesn't
    # mat file doesn't contain annotation field for reaction and metabolites. smbl does
    # 
    # add subsystem information into sbml model
    # 
    mat_model = cobra.io.load_matlab_model(mat_file)
    sbml_model = cobra.io.read_sbml_model(sbml_file)
    
    for rxn in mat_model.reactions:
        rxn_xml = sbml_model.reactions.get_by_id(rxn.id)
        rxn_xml.subsystem = rxn.subsystem
    
    return sbml_model

def retrive_from_uniprot(uniprot_id):
    # get fasta content for given uniprot_id, return as string format
    url = 'https://www.uniprot.org/uniprot/{0}.fasta'.format(uniprot_id)
    response = urllib.request.urlopen(url)
    rec = response.read()
    return rec

def download_genes(model,outname):
    #### Download genes in template model based on uniprot ids
    # 1. load ids that have been done
    done_ids = [line.split('|')[1] for line in open(outname) if line.startswith('>')]
    
    fhand = open(outname,'a')
    for g in model.genes:
        uniprot_id = g.annotation.get('uniprot')
        if uniprot_id is None: 
            print(g.id,g.name,g.annotation)
            continue
        if uniprot_id in done_ids: continue
            
        try:
            rec = retrive_from_uniprot(uniprot_id)
            fhand.write(str(rec)[2:-1].replace('\\n','\n'))
        except: print('Failed to retrive sequence for',uniprot_id)

    fhand.close()
    
def update_fasta_format(model,fastafile):
    # Replace the uniprot gene ids in fastafile with gene ids in the model
    
    if '.fasta' in fastafile: outname = fastafile.replace('.fasta','_formated.fasta')
    elif '.fa' in fastafile: outname = fastafile.replace('.fa','_formated.fasta')
    else: outname = fastafile + '_formated.fasta'
    fhand = open(outname,'w')
    
    # load uniprotid-geneid
    uni2gene = dict()
    for gene in model.genes:
        uni = gene.annotation.get('uniprot')
        if uni is not None: uni2gene[uni] = gene.id
    
    for rec in SeqIO.parse(fastafile,'fasta'):
        uni = rec.id.split('|')[1]
        fhand.write('>{0}\n{1}\n'.format(uni2gene[uni],rec.seq))
    fhand.close()
    
def add_coverage(combined_fasta,blast_results):
    # combined_fasta: combined fasta file, which contains all sequnces from both query and target proteome.
    # blast_results: tsv file produced by blastp
    # 
    # Output: add a new column: coverage at the end of the blast_results.
    # Coverage is caculated as: length of aligned query sequence/length of query sequence * 100
    # 
    
    # Load sequence lenght
    Lenghs = dict()
    for rec in SeqIO.parse(combined_fasta,'fasta'): Lenghs[rec.id] = len(rec.seq)
    
    
    # calculate coverage
    outfile = blast_results.replace('.tab','_cov.tab')
    if outfile == blast_results: outfile += '_cov.tab'
    
    fhand = open(outfile,'w')
    for line in open(blast_results):
        cont = line.split()
        aln_lengh = float(cont[7])-float(cont[6])
        cov = aln_lengh/Lenghs[cont[0]]*100
        fhand.write(line.strip()+'\t{0}\n'.format(cov))
    fhand.close()
    
    # replace the orginal blast results file with the new one.
    os.system('mv {0} {1}'.format(outfile,blast_results))
    
def do_blast(query_fasta,target_fasta):
    # query_fasta: fasta file for the organsims to be modelled
    # target_fasta: fasata file that contains enzyme sequences in the template model
    
    # combine two fasta files, do blast. Use default evalue. The blast results will be filetered based on evalue
    # in the BBH extraction step
    
    # build combined fasta file for blast
    fhand = open('input.fa','w')
    for rec in SeqIO.parse(query_fasta,'fasta'):
        fhand.write('>query|{0}\n{1}\n'.format(rec.id,rec.seq))
    
    for rec in SeqIO.parse(target_fasta,'fasta'):
        fhand.write('>target|{0}\n{1}\n'.format(rec.id,rec.seq))
    
    fhand.close()
    
    blast_cmd = '''mkdir blast_tmps
    mv input.fa blast_tmps
    makeblastdb -dbtype prot -in blast_tmps/input.fa  -out blast_tmps/DB
    blastp -query blast_tmps/input.fa -db blast_tmps/DB -outfmt 6  -max_target_seqs 500 -max_hsps 1  -num_threads 1 -out ../Results/bidirectional_blast.tab
    '''
    os.system(blast_cmd)
    add_coverage('blast_tmps/input.fa','../Results/bidirectional_blast.tab')
    os.system('rm -r blast_tmps')
    
    
    
def extract_BBHs(blastfile,BBH_evalue,report=False,coverage=45):
    # BBH_evalue: like 1e-10
    # report: if True, print out the nubmer of BBHs
    # coverage: the lower bound of query coverage to determine a true hit. from 0-100
    # 
    
    target_query = dict()
    query_target = dict()
    for line in open(blastfile):
        cont = line.split()
        evalue = float(cont[10])
        cov = float(cont[12])
        # filter based on evalue and coverage
        if evalue > BBH_evalue: continue
        if cov < coverage: continue
        
        if cont[0].startswith('target|') and cont[1].startswith('query|'):
            target_gene = cont[0].split('|')[1]
            query_gene = cont[1].split('|')[1]
            if target_query.get(target_gene,(None,1000))[1]>evalue: target_query[target_gene] = (query_gene,evalue)
        elif cont[0].startswith('query|') and cont[1].startswith('target|'):
            target_gene = cont[1].split('|')[1]
            query_gene = cont[0].split('|')[1]
            if query_target.get(query_gene,(None,1000))[1]>evalue: query_target[query_gene] = (target_gene,evalue)
        else: None

    # 3. find BBH, the overlaped gene pairs
    BBH = list()
    for target_gene, (query_gene,evalue) in target_query.items():
        if query_target.get(query_gene,(None,1000))[0] == target_gene: BBH.append((query_gene,target_gene))
    
    if report:
        print('Number of homologs for query genes (evalue = {0}):'.format(BBH_evalue),len(query_target))
        print('Number of homologs for target genes (evalue = {0}):'.format(BBH_evalue),len(target_query))
        print('Number of BBHs (evalue = {0}):'.format(BBH_evalue),len(BBH))
    return BBH

def get_all_rxns_in_BBH(template_model, BBHs):
    # model file: template model file
    # BBH, best bidirectional hit

    # get all reactions
    candidate_rxn_ids = list()
    candidate_rxns = list()
    for query_gene,target_gene in BBHs:
        gene = template_model.genes.get_by_id(target_gene)
        for rxn in gene.reactions:
            if rxn.id not in candidate_rxn_ids: 
                candidate_rxn_ids.append(rxn.id)
                candidate_rxns.append(rxn)
    print('Number of candiate reactions:',len(candidate_rxn_ids))
    
    return candidate_rxns

def parse_gr(gr):
    # return a list of genes in gr rule
    genes = list(set([item for item in gr.replace('(','').replace(')','').split() if item not in ['and','or']]))
    genes.sort()
    return genes

def parse_gr_parts(gr):
    if ' and ' in gr and ' or ' not in gr: gene_parts = ['( {0} )'.format(gr)]
    else:
        complexes = re.findall('(\(.*?\))',gr)
        gr1 = gr
        for i in range(len(complexes)):
            gr1 = gr1.replace(complexes[i],'comp_{0}'.format(i))
        gene_parts = parse_gr(gr1)

        for i in range(len(gene_parts)):
            if 'comp_' in gene_parts[i]: gene_parts[i] = complexes[int(gene_parts[i].split('_')[1])]
    return list(set(gene_parts))

def missing_type(gene_part):
    # determin if a gene/coplex is totally, partly missing or complete
    assert ' or ' not in gene_part
    genes = parse_gr(gene_part)
    missing_parts = re.findall('_missing',gene_part)
    
    if len(missing_parts) == 0: tpe = 'complete' # complete
    elif len(missing_parts) == len(genes): tpe = 'missing'
    else: tpe = 'partial'
    return tpe

def remove_repeated_gene_in_and(gene_part):
    # for case: A and A
    assert ' or ' not in gene_part
    genes = parse_gr(gene_part)
    
    new_part = ''
    for gene in genes: new_part += gene + ' and '
    new_part = new_part[:-5]
    
    # in case the final gr rule is (A and B), remove parathesis
    if len(genes) > 1: new_part = '( {0} )'.format(new_part) 
    return new_part

def remove_repeated_gene_parts(gr):
    gene_parts = parse_gr_parts(gr)
    
    # remove repeated genes in each part 
    new_gene_parts = list()
    for gene_part in gene_parts:
        new_gene_part = remove_repeated_gene_in_and(gene_part)
        if new_gene_part not in new_gene_parts: new_gene_parts.append(new_gene_part)
    
    new_gr = ''
    for gene_part in new_gene_parts: 
        new_gr += gene_part + ' or '
    new_gr = new_gr[:-4]
    
    # in case the final gr rule is (A and B), remove parathesis
    if len(new_gene_parts) == 1 and '(' in new_gr: new_gr = new_gr[1:-1].strip()
    return new_gr

def update_and_or_case(gr):
    genes1 = parse_gr_parts(gr)
    tpes = list()
    for g in genes1: tpes.append(missing_type(g))

    # case 1: if one of parts is complete, remove all parts with are missing or partial
    if 'complete' in tpes:
        new_gr = ''
        for i in range(len(genes1)):
            if tpes[i] == 'complete': new_gr += genes1[i] + ' or '
    # case 2: if none of parts is complete, only remove those ones that are missing
    else:
        new_gr = ''
        for i in range(len(genes1)):
            if tpes[i] != 'missing': new_gr += genes1[i] + ' or '
    new_gr = new_gr[:-4]
    
    # in case the final gr rule is (A and B), remove parathesis
    if 'or' not in new_gr and '(' in new_gr: new_gr = new_gr[1:-1].strip()
    
    new_gr = remove_repeated_gene_parts(new_gr)
    return new_gr

def update_gr(gr,BBHs):
    # replace the targetli gene ids in gr rule with the ones in querymonas based on BBH. For those targetli genes
    # wihout BBH in querymonas, add a tag '_missing' at the end of the gene id
    
    target2query_map = dict()
    for query_gene,target_gene in BBHs: target2query_map[target_gene] = query_gene
    
    genes = parse_gr(gr)
    new_gr = gr
    for gene in genes:
        query_gene = target2query_map.get(gene)
        if query_gene is None: query_gene = gene+'_missing'
        new_gr = new_gr.replace(gene,query_gene)
    
    # refine new gr rule. 
     
    # case 1: if 'missing' not in gr rule, do not need to update
    if '_missing' not in new_gr: refined_gr = new_gr
    
    # case 2: only 'or' relationships are in gr rule, remove all genes with '_missing' tag
    elif ' and ' not in new_gr: 
        new_genes = [item for item in parse_gr(new_gr) if '_missing' not in item]
        refined_gr = ''
        for g in new_genes: refined_gr += g + ' or '
        refined_gr = refined_gr[:-4]
    
    # case 3: if 'or' not in gr rule, keep missing genes
    elif 'or' not in new_gr: refined_gr = new_gr
    
    # case 4: if (A and B) or (C and D_missing) or E. In such or structure, if one of the gene/gene complex is 
    #         complete, then remove all others with missing gene(s). if all parts are incompelte, only removes 
    #         those part which are totally missing.
    elif 'and' in new_gr and 'or' in new_gr: refined_gr = update_and_or_case(new_gr)

    # Ohter cases: export to a txt file and manully edit
    else:  refined_gr = new_gr
    return refined_gr

def report_model_status(model):
    print('Number of reactions:',len(model.reactions))
    print('Number of metabolits:',len(model.metabolites))
    print('Number of compartments:',len(model.compartments))
    k = 0
    for gene in model.genes:
        if 'missing' in gene.id: k += 1
    print('Number of genes:',len(model.genes))
    print('Number of missing genes:',k)
    
    k = 0
    for rxn in model.reactions:
        if '_missing' in rxn.gene_reaction_rule: k += 1
    print('Number of reactions with missing genes:',k)
    
def build_model_from_template(candidate_rxns,BBHs,report=False):
    rxns_to_add = list()
    for rxn in candidate_rxns:
        new_gr = update_gr(rxn.gene_reaction_rule,BBHs)
        new_rxn = rxn.copy()
        new_rxn.gene_reaction_rule = new_gr
        rxns_to_add.append(new_rxn)
    model_from_template = cobra.Model('model_from_template')
    model_from_template.add_reactions(rxns_to_add)
    if report: report_model_status(model_from_template)
    return model_from_template

def update_missing_genes(model,homolog_evalue,report=False,coverage=45):
    # homolog_evalue: evalue used in blast to infer homology relationships
    # model: tmeplate model from BBH, which in most cases contains missing genes in gene-reaction rules
    
    # load best hit based on homolog_evalue
    best_hits = dict()
    # best_hits = {eco_gene_id:(halo_gene,evalue)}
    for line in open('../Results/bidirectional_blast.tab'):
        cont = line.split()
        evalue = float(cont[10])
        cov = float(cont[12])
        
        # filter hits based on evalue and coverage
        if evalue > homolog_evalue: continue
        if cov < coverage: continue
        
        if line.startswith('target|') and cont[1].startswith('query'):
            target = cont[0].split('|')[1]
            query = cont[1].split('|')[1]
            best_hit, best_evalue = best_hits.get(target,(None,0))
            if best_hit is None or best_evalue > evalue: best_hits[target] = (query,evalue)
    
    # update missing genes in the model. Replace the missing genes with homologs found in above step
    
    updated_rxns = list()
    for rxn in model.reactions:
        rxn_copy = rxn.copy()
        gr = rxn_copy.gene_reaction_rule

        if '_missing' in gr:
            genes = parse_gr(gr)
            for gene in genes:
                if '_missing' in gene: 
                    homolog = best_hits.get(gene.replace('_missing',''),(None,None))[0]
                    if homolog is not None: gr = gr.replace(gene,homolog)
        
        gr = remove_repeated_gene_parts(gr)
        gr = update_and_or_case(gr)
        rxn_copy.gene_reaction_rule = gr

        updated_rxns.append(rxn_copy)
    
    # create a new model
    new_model = cobra.Model('model_from_template')
    new_model.add_reactions(updated_rxns)
    if report: report_model_status(new_model)
    return new_model

def remove_rxns_with_missing_genes(model,report=True):
    rxns = list()
    for rxn in model.reactions:
        if '_missing' in rxn.gene_reaction_rule: continue
        rxns.append(rxn)
    new_model = cobra.Model('model_from_template')
    new_model.add_reactions(rxns)
    if report: report_model_status(new_model)
    return new_model

def load_metacyc_mat(model_file):
    # model_file: a mat file of the model. This model is the model that is reconstructed from MetaCyc with RAVEN.
    # return a dictionary, with filed name as key, field values as value
    # like data['id'] = 'metacycmodel'
    #
    
    ###### 
    # load mat model with cobrapy
    model = cobra.io.load_matlab_model(model_file)
    
    
    ###### 
    # load mnx of metabolites
    meta2mnx = dict()
    dfmnx = pd.read_csv('../../../ComplementaryData/MNX2metacyc.tsv',index_col=1,sep='\t')

    for ind in dfmnx.index:meta2mnx[ind] = meta2mnx.get(ind,[]) + [dfmnx.loc[ind,'MNX']]
    
    
    ######
    # load annotation fields
    data = dict()
    
    mat = mat = loadmat(model_file)
    meta_vars = {"__globals__", "__header__", "__version__"}
    for item in mat.keys(): 
        if item not in meta_vars: model_id = item
    
    # field names
    names = mat[model_id].dtype.names
    
    # save data into a dictionary
    for i,item in enumerate(mat[model_id][0][0]):
        if names[i] in ['rxns','rxnNames','eccodes','rxnReferences','mets','metNames','metFormulas','inchis']:
            data[names[i]] = [i[0][0] if len(i[0])!=0 else '' for i in item]
        
        elif names[i] in ['subSystems']:
            lst = []
            for sub_item in item:
                try: lst.append([k[0] for k in sub_item[0][0]])
                except: lst.append([])
            data[names[i]] = lst
        elif names[i] in ['rxnMiriams','metMiriams']:
            lst = []
            for sub_item in item:
                anno = dict()
                try:
                    sub_item = sub_item[0][0][0]
                    for j in range(len(sub_item['name'])):
                        name = sub_item['name'][j][0][0]
                        value = sub_item['value'][j][0][0]
                        anno[name] = str(value)
                except:pass
                lst.append(anno)
            data[names[i]] = lst
        elif names[i] in ['rxnConfidenceScores','metCharges']:
            data[names[i]] = [k[0] for k in item]
        else: pass
    
    print('Loaded fields:',data.keys())
    
    ######
    # there are some metabolites formula cannot be correctly parsed by cobrapy. like ZN1 from RAVEN, which will
    # be treated as Z1N1. Most of such metabolites are ions
    for i,item in enumerate(data['metFormulas']):
        for element in ['ZN','MN','CU','MG','FE','CL','BR','AS','CA','CD','MO','TE','SE','NA']: 
            if element in item: item = item.replace(element,element.lower().title())
            data['metFormulas'][i] = item
    
    ######
    # assign annotation to the cobrapy model
    ## 1. assign rxn annotation and reaction subsystem
    
    for i in range(len(data['rxns'])):
        rxn_id = data['rxns'][i]
        anno = data['rxnMiriams'][i]
        
        # ec number
        ec = data['eccodes'][i]
        ec = ec.replace('ec-code/','')

        if len(ec) == 0: pass
        elif len(ec.split(';')) == 1: anno['ec-code'] = ec
        else: anno['ec-code'] = ec.split(';')
        
        # metacyc id
        if anno.get('kegg.reaction') != rxn_id:
            anno['biocyc'] = 'META:'+rxn_id
            anno['id_source'] = 'biocyc'
        else: anno['id_source'] = 'kegg'
        
        rxn = model.reactions.get_by_id(rxn_id)
        rxn.annotation = anno
        
        # assign reaction subsystem
        subss = ''
        for item in data['subSystems'][i]: 
            subss += item +';'
        if subss.endswith(';'): subss = subss[:-1]
        rxn.subsystem = subss
    
    ## assign metabolites annotation
    for i in range(len(data['mets'])):
        met_id = data['mets'][i]
        anno = data['metMiriams'][i]
        
        # assign metacyc id
        if anno.get('kegg.compound') != met_id:
            anno['biocyc'] = 'META:'+met_id
            anno['id_source'] = 'biocyc'
        else: anno['id_source'] = 'kegg'
        
        # assign mnx id
        mnx_id = meta2mnx.get(met_id)
        if mnx_id is not None: anno['metanetx.chemical'] = mnx_id
        
        met = model.metabolites.get_by_id(met_id)
        met.annotation = anno
        
        # assign charge
        met.charge = data['metCharges'][i]
        
        # assign formula
        met.formula = data['metFormulas'][i]

    return model

def print_keggrxn(rxn):
    eq = rxn.reaction
    for met in rxn.metabolites:
        eq = eq.replace(met.id,met.name)
    print('KEGG:',eq)
    
def update_rxn(model,old_met_id,new_met_id):
    old_met = model.metabolites.get_by_id(old_met_id)
    for old_rxn in old_met.reactions:
        
        # create a new reaction with replaced metabolite
        new_rxn = cobra.Reaction(old_rxn.id)
        new_rxn.name = old_rxn.name
        new_rxn.subsystem = old_rxn.subsystem
        new_rxn.lower_bound = old_rxn.lower_bound
        new_rxn.upper_bound = old_rxn.upper_bound
        
        new_met = cobra.Metabolite(new_met_id,name=new_met_id,compartment=old_met.compartment)
        
        coeffs = dict()
        for met in old_rxn.metabolites:
            if met.id == old_met.id: coeffs[new_met] = old_rxn.get_coefficient(met)
            else: coeffs[met] = old_rxn.get_coefficient(met)
        new_rxn.add_metabolites(coeffs)
        new_rxn.gene_reaction_rule = old_rxn.gene_reaction_rule 
        
        # remove old reactions and add new reaction
        model.remove_reactions([old_rxn])
        model.add_reaction(new_rxn)
    
    # findally, remove old metabolites
    model.remove_metabolites([old_met])
    return model


def report_model_status(model):
    print('Number of reactions:',len(model.reactions))
    print('Number of metabolits:',len(model.metabolites))
    print('Number of compartments:',len(model.compartments),model.compartments)
    k = 0
    for gene in model.genes:
        if 'missing' in gene.id: k += 1
    print('Number of genes:',len(model.genes))
    print('Number of missing genes:',k)
    
    k = 0
    for rxn in model.reactions:
        if '_missing' in rxn.gene_reaction_rule: k += 1
    print('Number of reactions with missing genes:',k)
    print()
    
    
def load_MNX(infile):
    mnx2meta_met = dict()
    mnx2kegg_met = dict()
    for line in open(infile):
        if line.startswith('#'):continue
        if line.startswith('kegg'): 
            cont = line.split()
            mnx2kegg_met[cont[1]] = mnx2kegg_met.get(cont[1],[]) + [cont[0].replace('kegg:','')]
            
        if line.startswith('metacyc'): 
            cont = line.split()
            mnx2meta_met[cont[1]] = mnx2meta_met.get(cont[1],[]) + [cont[0].replace('metacyc:','')]
    
    print('Metacyc',len(mnx2meta_met))
    print('KEGG',len(mnx2kegg_met))
    
    return mnx2meta_met,mnx2kegg_met

## Build Universal model from iML1515 based recon
def load_met_rxn_ids_in_raven_model(model):
    # model, raven model
    meta_raven_ids_rxn, kegg_raven_ids_rxn = {},{}
    meta_raven_ids_met, kegg_raven_ids_met = {},{}
    
    for rxn in model.reactions:
        if rxn.annotation.get('id_source','biocyc') == 'biocyc': meta_raven_ids_rxn[rxn.id.split('_')[0]] = True
        if rxn.annotation.get('id_source','biocyc') == 'kegg'  : kegg_raven_ids_rxn[rxn.id.split('_')[0]] = True
    
    for met in model.metabolites:
        if met.annotation.get('id_source','biocyc') == 'biocyc': meta_raven_ids_met[rxn.id.split('_')[0]] = True
        if met.annotation.get('id_source','biocyc') == 'kegg'  : kegg_raven_ids_met[rxn.id.split('_')[0]] = True
            
    return meta_raven_ids_rxn, kegg_raven_ids_rxn, meta_raven_ids_met, kegg_raven_ids_met

def convert_met_id(met,mnx2meta_met, mnx2kegg_met,meta_raven_ids, kegg_raven_ids):
    '''
    met         : met in iML1515
    mnx2meta_met: {mnx_id: [metacyc ids]}
    mnx2kegg_met: {mnx_id: [kegg ids]}
    meta_ids    : {metcyc_id:True}, the metacyc ids in the RAVEN model
    kegg_ids    : {kegg_id: True}, the kegg ids in the RAVEN model
    '''
    # Firstly match existing metacyc ids and kegg ids
    meta_ids = None#met.annotation.get('biocyc')
    if meta_ids is None: meta_ids = []
    elif type(meta_ids) is list: meta_ids = [item.replace('META:','') for item in meta_ids]
    else: meta_ids = [meta_ids.replace('META:','')]
    
    kegg_ids = met.annotation.get('kegg.compound')
    if kegg_ids is None: kegg_ids = []
    elif type(kegg_ids) is list: kegg_ids = kegg_ids
    else: kegg_ids = [kegg_ids]
    
    # match based on mnx id
    mnx_id   = met.annotation.get('metanetx.chemical')
    if mnx_id is not None: 
        meta_ids = meta_ids + mnx2meta_met.get(mnx_id,[])
        kegg_ids = kegg_ids + mnx2kegg_met.get(mnx_id,[])
    
    new_id = None
    id_source = 'bigg'
    if len(meta_ids) != 0: 
        new_id = meta_ids[0]
        id_source = 'biocyc'
        for meta_id in meta_ids:
            if meta_raven_ids.get(meta_id,False): new_id = meta_id
    elif len(kegg_ids) != 0: 
        id_source = 'kegg'
        new_id = kegg_ids[0]
        for kegg_id in kegg_ids:
            if kegg_raven_ids.get(kegg_id,False): new_id = kegg_id
    else: pass
    
    new_met = met.copy()
    new_met.annotation['id_source'] = id_source
    if new_id is not None: new_met.id = new_id + '_' + met.compartment
    return new_met
        
    
def convert_rxn_id(rxn,mnx2meta_met, mnx2kegg_met,meta_raven_ids, kegg_raven_ids):
    '''
    rxn         : rxn in iML1515
    mnx2meta_met: {mnx_id: [metacyc ids]}
    mnx2kegg_met: {mnx_id: [kegg ids]}
    meta_ids    : {metcyc_id:True}, the metacyc ids in the RAVEN model
    kegg_ids    : {kegg_id: True}, the kegg ids in the RAVEN model
    '''
    # Firstly match existing metacyc ids and kegg ids
    meta_ids = rxn.annotation.get('biocyc')
    if meta_ids is None: meta_ids = []
    elif type(meta_ids) is list: meta_ids = [item.replace('META:','') for item in meta_ids]
    else: meta_ids = [meta_ids.replace('META:','')]
    
    kegg_ids = rxn.annotation.get('kegg.reaction')
    if kegg_ids is None: kegg_ids = []
    elif type(kegg_ids) is list: kegg_ids = kegg_ids
    else: kegg_ids = [kegg_ids]
    
    # match based on mnx id
    mnx_id   = rxn.annotation.get('metanetx.reaction')
    if mnx_id is not None: 
        meta_ids = meta_ids + mnx2meta_met.get(mnx_id,[])
        kegg_ids = kegg_ids + mnx2kegg_met.get(mnx_id,[])
    
    new_id = None
    id_source = 'bigg'
    if len(meta_ids) != 0: 
        new_id = meta_ids[0]
        id_source = 'biocyc'
        for meta_id in meta_ids:
            if meta_raven_ids.get(meta_id,False): new_id = meta_id
    elif len(kegg_ids) != 0: 
        new_id = kegg_ids[0]
        id_source = 'kegg'
        for kegg_id in kegg_ids:
            if kegg_raven_ids.get(kegg_id,False): new_id = kegg_id
    else: pass
    
    new_rxn = rxn.copy()
    if new_id is not None: new_rxn.id = new_id
    return new_rxn

def convert_eco_to_metacyc_kegg(eco_model,raven_model):
    mnx2meta_met, mnx2kegg_met = load_MNX('../../../ComplementaryData/chem_xref.tsv')
    mnx2meta_rxn, mnx2kegg_rxn = load_MNX('../../../ComplementaryData/reac_xref.tsv')
    
    meta_raven_ids_rxn, kegg_raven_ids_rxn, meta_raven_ids_met, kegg_raven_ids_met = load_met_rxn_ids_in_raven_model(raven_model)
    
    converted_rxn_list = []
    for rxn in eco_model.reactions:
        new_rxn = convert_rxn_id(rxn,mnx2meta_rxn, mnx2kegg_rxn,meta_raven_ids_rxn, kegg_raven_ids_rxn)
        new_mets = {}
        for met in new_rxn.metabolites:
            new_met = convert_met_id(met,mnx2meta_met, mnx2kegg_met,meta_raven_ids_met, kegg_raven_ids_met)
            new_mets[new_met] = new_rxn.metabolites[met]
        new_rxn.subtract_metabolites(new_rxn.metabolites)
        new_rxn.add_metabolites(new_mets)
        converted_rxn_list.append(new_rxn)
    
    
    # process reduandant reaction ids
    # ASPtpp
    # 5-10-METHENYL-THF_c + CPD-15815_c --> 5-FORMYL-THF_c + PROTON_c
    # 5-FORMYL-THF_c + PROTON_c --> 5-10-METHENYL-THF_c + CPD-15815_c
    rxn_ids = {}
    updated_rxn_list = []
    for rxn in converted_rxn_list: rxn_ids[rxn.id] = rxn_ids.get(rxn.id,[]) + [rxn]
    for rxn_id, rxns in rxn_ids.items():
        if len(rxns)>1: 
            #print(rxn_id)
            for i, rxn in enumerate(rxns): 
                rxn.id = rxn.id + '_copy{0}'.format(i+1)
                #print(rxn.id,rxn.reaction)
                updated_rxn_list.append(rxn)
        else: updated_rxn_list.extend(rxns)

    converted_model = cobra.Model()
    converted_model.add_reactions(updated_rxn_list)
    return converted_model

def refresh_model(old_model):
    model = cobra.Model(old_model.id)
    rxns = [rxn.copy() for rxn in old_model.reactions]
    model.add_reactions(rxns)
    report_model_status(model)
    return model

def save_pickle_model(model,outfile):
    pickle.dump(refresh_model(model),open(outfile,'wb'))

    
####GAP FILLING
def test_if_model_can_produce(met,halomodel):
    # add sink reaction
    rxn_id = 'tmp'
    rxn_sink = cobra.Reaction(rxn_id)
    rxn_sink.add_metabolites({met:-1})
    halomodel.add_reaction(rxn_sink)
    halomodel.objective = 'tmp'
    f = halomodel.optimize().objective_value
    halomodel.remove_reactions([rxn_sink])
    print('maximal flux:',met.id,f)
    
    
def check_if_model_can_produce_mets_in_biomass_one_by_one(model):
    for rxn_id in ['Biomass_v1','Protein_synthesis','DNA_synthesis','RNA_synthesis','ions_pool','lipids_synthesis']:
        print(rxn_id)
        for met in model.reactions.get_by_id(rxn_id).reactants:
            with model:
                # add sink reaction
                    rxn_id = 'tmp'
                    rxn_sink = cobra.Reaction(rxn_id)
                    rxn_sink.add_metabolites({met:-1})
                    model.add_reaction(rxn_sink)
                    model.objective = 'tmp'
                    f = model.optimize().objective_value
                    print('  ',met.id,met.name,f)
        print()
                

        
        
def do_gap_filling(met,halmodel,universal,**kwargs):
 
    rxn_id = 'tmp'
    try: 
        rxn_sink = halmodel.reactions.get_by_id(rxn_id)
        halmodel.remove_reactions([rxn_sink])
    except: None
    rxn_sink = cobra.Reaction(rxn_id)
    rxn_sink.add_metabolites({met:-1})
    halmodel.add_reaction(rxn_sink)
    halmodel.objective = 'tmp'
    f1 = halmodel.optimize().objective_value
    
    added_genes = dict()# {'rxn_id':'genes'}
    print(met.id)
    print('Before gap-filling:',f1)
    if f1 != 0: print('No need for gap-filling')
    else:
        try: 
            with universal:
                solution = cobra.flux_analysis.gapfill(halmodel, universal,**kwargs)
            all_rxns = dict()
            rxns_to_add = list()
            for rxn in solution[0]:
                try: halmodel.reactions.get_by_id(rxn.id)
                except: 
                    rxn_clean = rxn.copy()
                    added_genes[rxn_clean.id] = rxn_clean.gene_reaction_rule
                    
                    # add _missing tag for those genes
                    new_gr = rxn_clean.gene_reaction_rule
                    for gene in parse_gr(rxn_clean.gene_reaction_rule):
                        new_gr = new_gr.replace(gene,gene+'_missing')
        
                    rxn_clean.gene_reaction_rule = new_gr
                    rxns_to_add.append(rxn_clean)
            halmodel.add_reactions(rxns_to_add)
            f2 = halmodel.optimize().objective_value
            print('After gap-filling:',f2)
            print('Added reactions:')
            for rxn in rxns_to_add:print('  ',rxn.id,rxn.reaction)
        except: print('Gap-filling failed for',met.id)

    halmodel.remove_reactions([rxn_sink])
    halmodel.repair()
    print('\n\n')
    return added_genes