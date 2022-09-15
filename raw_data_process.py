from glob import glob
import pandas as pd
from collections import defaultdict
import json


def fasta_reader(fasta_file_path):

    with open(fasta_file_path, 'r') as file_open:
        file_split = file_open.read().split('\n>')

    return {each.split('\n')[0].split('|')[1]: ''.join(each.split('\n')[1:]) for each in file_split}


def protein_info_from_fasta(fasta_path):
    """
    get protein name, gene name, entry name, and description
    :param fasta_path:
    :return:
    """
    info_dict = {}
    with open(fasta_path,'r') as f:
        for line in f:
             if line.startswith('>'):
                protein_id = line.split('|')[1]
                cls = line.split('|')[0].split('>')[1]
                # print (protein_id)
                description = ' '.join(line.split('OS=')[0].split(' ')[1:])

                gene_name = line.split('GN=')[1].split(' ')[0].rstrip('\n') if 'GN=' in line else 'N/A'
                info_dict[protein_id] = (gene_name,description,cls)
    return info_dict


def txt_reader(txt_file,matrix_protein_dict):
    prot_psm_dict = defaultdict(list)
    prot_hyperscore = defaultdict(float)
    with open(txt_file,'r',newline='\r\n') as f_o:

        for line in f_o:
            line_split = line.split('\t')
            if ';' in line_split[1]:

                prot_list = [each.split('|')[1] for each in line_split[1].split(';')]
            else:

                prot_list = [line_split[1].split('|')[1]]
            psm,hyper_score = line_split[0][2:-2],float(line_split[-1])
            for prot in prot_list:
                if prot in matrix_protein_dict:  # only include matrix protein
                    prot_psm_dict[prot].append(psm)
                    prot_hyperscore[prot] += hyper_score
    return prot_psm_dict, prot_hyperscore


def psm_species_counter(text_file_list, protein_info_dict,gene_species_dict):
    species_psm_count = defaultdict(int)
    for text_file in text_file_list:
        print(text_file)
        with open(text_file, 'r', newline='\r\n') as f_o:

            for line in f_o:
                line_split = line.split('\t')
                if ';' in line_split[1]:

                    prot_list = [each.split('|')[1] for each in line_split[1].split(';')]
                else:

                    prot_list = [line_split[1].split('|')[1]]
                for prot in prot_list:
                    if prot in protein_info_dict:
                        species = gene_species_dict[protein_info_dict[prot][0]]
                        species_psm_count[species]+=1
    return species_psm_count


def nsaf(prot_psm_dict,protein_seq_dict):

    nsaf_dict = {}
    nsaf_total = sum([len(prot_psm_dict[each])/len(protein_seq_dict[each]) for each in prot_psm_dict])

    for each in prot_psm_dict:
        nsaf_single = len(prot_psm_dict[each])/len(protein_seq_dict[each])/nsaf_total
        nsaf_dict[each] = nsaf_single
    return nsaf_dict


def table_assemble(txt_file,annontation_dict,protein_seq_dict,protein_info_dict,matrix_protein_info_dict):
    data = []
    columns = ["uniprot_id","gene_name","species","tissue",
               "organ_system","sample_type","repository_id",
               "matrisome_category","matrisome_class",
               "file_name","reference_doi","protein_description",
               "note","total_psm","hyperscore_sum","NSAF","seq_cov_file"]

    raw_file = txt_file.split('\\')[-1].replace('.txt','.raw')
    if raw_file not in annontation_dict:
        raw_file = txt_file.split('\\')[-1].replace('.txt','.RAW')
    annotation_sub_dict = annontation_dict[raw_file]

    project, doi, raw_species, sys, tissue, description = annotation_sub_dict["Project"],annotation_sub_dict["DOI"],\
                                                          annotation_sub_dict["Species"], annotation_sub_dict["Sys"], \
                                                          annotation_sub_dict["Tissue"], annotation_sub_dict["Description"]

    prot_psm_dict, prot_hyperscore = txt_reader(txt_file,protein_seq_dict)
    prot_psm_count_dict = {prot:len(prot_psm_dict[prot]) for prot in prot_psm_dict}
    nsaf_dict = nsaf(prot_psm_dict,protein_seq_dict)

    for prot in prot_psm_dict:
        gene = protein_info_dict[prot][0]
        matrix_cat = matrix_protein_info_dict[gene]["Category"] if gene in matrix_protein_info_dict else ''
        matrix_sub_cat = matrix_protein_info_dict[gene]["Sub"] if gene in matrix_protein_info_dict else ''
        prot_species = matrix_protein_info_dict[gene]["Species"] if gene in matrix_protein_info_dict else ''

        protein_descript = protein_info_dict[prot][1]
        data.append([prot,gene,prot_species,tissue,sys,description,project,matrix_cat,matrix_sub_cat,raw_file,doi,
                     protein_descript,'',prot_psm_count_dict[prot],prot_hyperscore[prot],nsaf_dict[prot],description+'_'+prot+'.html'])

    df = pd.DataFrame(data,columns=columns)
    print (f"output to {txt_file.replace('.txt','_summary.csv')}")
    return df.to_csv(txt_file.replace('.txt','_summary.tsv'),sep='\t')


if __name__=='__main__':

    base_path = 'F:/matrisomedb2.0/MDB2/result/'
    files = glob(base_path+'/**/*.txt',recursive=True)

    matrix_info_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json'))
    protein_info_dict = protein_info_from_fasta('F:/matrisomedb2.0/mat.fasta')
    gene_species_dict = {gene:matrix_info_dict[gene]['Species'] for gene in matrix_info_dict}

    print (psm_species_counter(files,protein_info_dict,gene_species_dict))

    annotation_dict = json.load(open('F:/matrisomedb2.0/annotation/matdb_dict.json'))
    for f in annotation_dict:
        for each in annotation_dict[f]:

            if annotation_dict[f][each][-1] == ' ':
                print (annotation_dict[f][each])

    annotation_dict = {each.split('/')[-1]:annotation_dict[each] for each in annotation_dict}
    print ([k for k in annotation_dict.keys()])

    protein_seq_dict = fasta_reader('F:/matrisomedb2.0/mat.fasta')


    for each_f in files:
        try:
            table_assemble(each_f,annotation_dict,protein_seq_dict,protein_info_dict,matrix_info_dict)
        except KeyError:
            print (f'{each_f} not in annotation file')


