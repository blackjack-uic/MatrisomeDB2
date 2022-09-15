import json
from collections import defaultdict
import os
from raw_data_process import fasta_reader
from glob import glob


def text_file_reader(sample_type_files_dict,base_path,matrix_protein_dict):
    """
    read multiple files from one sample type
    :return:
    """
    sample_protein_psm_dict = {}
    for sample in sample_type_files_dict:
        print (sample)
        prot_psm_dict = defaultdict(list)
        file_list = []
        for f in sample_type_files_dict[sample]:
            try:
                f_path = glob(base_path+'**/'+f,recursive=True)[0]
                file_list.append(f_path)
            except IndexError:
                continue

        print (sample, f'{len(file_list)} files')

        for each_f in file_list:
            try:
                with open(each_f, 'r', newline='\r\n') as f_o:
                    for line in f_o:
                        line_split = line.split('\t')
                        if ';' in line_split[1]:

                            prot_list = [each.split('|')[1] for each in line_split[1].split(';')]
                        else:

                            prot_list = [line_split[1].split('|')[1]]
                        psm = line_split[0][2:-2]
                        for prot in prot_list:
                            if prot in matrix_protein_dict:  # only include matrix protein
                                prot_psm_dict[prot].append(psm)
                    print (f'{each_f} succeed')
            except FileNotFoundError:
                print (f"{each_f} failed")
                continue

        sample_protein_psm_dict[sample] = prot_psm_dict
    return json.dump(sample_protein_psm_dict,open('F:/matrisomedb2.0/sample_protein_psm_dict_3.json','w'))


def txt_reader_all(txt_file_list,matrix_protein_dict):
    prot_psm_dict = defaultdict(list)
    for txt_file in txt_file_list:
        with open(txt_file,'r',newline='\r\n') as f_o:

            for line in f_o:
                line_split = line.split('\t')
                if ';' in line_split[1]:

                    prot_list = [each.split('|')[1] for each in line_split[1].split(';')]
                else:

                    prot_list = [line_split[1].split('|')[1]]
                psm = line_split[0][2:-2]
                for prot in prot_list:
                    if prot in matrix_protein_dict:  # only include matrix protein
                        prot_psm_dict[prot].append(psm)

    return prot_psm_dict

def txt_all_psms(txt_file_list):
    psm_list = []
    for txt_file in txt_file_list:
        print (txt_file)
        with open(txt_file,'r',newline='\r\n') as f_o:

            for line in f_o:
                line_split = line.split('\t')
                psm = line_split[0][2:-2]
                psm_list.append(psm)
    return psm_list


def all_matrix_psms(txt_file_list,matrix_protein_dict):
    psm_count = 0
    for txt_file in txt_file_list:
        print(txt_file)
        with open(txt_file, 'r', newline='\r\n') as f_o:

            for line in f_o:
                line_split = line.split('\t')

                if ';' in line_split[1]:

                    prot_list = [each.split('|')[1] for each in line_split[1].split(';')]
                else:

                    prot_list = [line_split[1].split('|')[1]]
                for prot in prot_list:
                    if prot in matrix_protein_dict:
                        psm_count+=1
                        break

    return psm_count

if __name__=="__main__":
    import pickle
    protein_seq_dict = fasta_reader('F:/matrisomedb2.0/mat.fasta')
    annotation_dict = json.load(open('F:/matrisomedb2.0/annotation/matdb_dict.json'))
    sample_type_files_dict = defaultdict(set)
    base_path = 'F:/matrisomedb2.0/MDB2/result/'

    files = glob(base_path + '/**/*.txt', recursive=True)
    psm_list = txt_all_psms(files)
    pickle.dump(psm_list,open('F:/matrisomedb2.0/all_psm.p','wb'),protocol=5)
    all_prot_psm_dict = txt_reader_all(files,protein_seq_dict)


    # for f_path in annotation_dict:
    #     if '.RAW' in f_path:
    #         txt_file = f_path.replace('.RAW','.txt')
    #     else:
    #         txt_file = f_path.replace('.raw','.txt')
    #     sample_type = annotation_dict[f_path]["Description"]
    #     sample_type_files_dict[sample_type].add(txt_file.split('/')[-1])
    # print (sample_type_files_dict['Normal lung ECM (Quantitative detergent solubility profiling)'])
    # count = 0
    # sample_type_psm_dict = json.load(open('F:/matrisomedb2.0/sample_protein_psm_dict_3.json','r'))
    # all_prot_psm_dict = defaultdict(list)
    # for samp in sample_type_psm_dict:
    #     for prot in sample_type_psm_dict[samp]:
    #         count+=1
    # print (count)
    #         for psm in sample_type_psm_dict[samp][prot]:
    #             all_prot_psm_dict[prot].append(psm)
    # print (len(all_prot_psm_dict['P02751']))
    # json.dump(all_prot_psm_dict,open('F:/matrisomedb2.0/global_protein_psm.dict_fromsample.json','w'))

    # all_prot_psm_dict = json.load(open('F:/matrisomedb2.0/global_protein_psm.dict.json','r'))
    # print (len(all_prot_psm_dict['P02751']))

