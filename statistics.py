import pickle
import json
from raw_data_process import protein_info_from_fasta,fasta_reader
from sample_psm_extract import all_matrix_psms
import numpy as np
import pandas as pd
from glob import glob
import matplotlib.pyplot as plt
import seaborn as sns

## sequence coverage calculation

matrix_info_dict = json.load(open('F:/matrisomedb2.0/annotation/mat_dict.json'))
protein_info_dict = protein_info_from_fasta('F:/matrisomedb2.0/mat.fasta')
matrix_protein_dict = fasta_reader('F:/matrisomedb2.0/mat.fasta')
gene_prot_dict = {protein_info_dict[prot][0]:prot for prot in protein_info_dict}

"""
global_prot_freq_dict = pickle.load(open('F:/matrisomedb2.0/glob_prot_freq_dict.p','rb'))

df = pd.DataFrame(index=[k for k in global_prot_freq_dict.keys()],columns=['Gene','Species','Category','Sub','Sequence coverage'])
for prot in global_prot_freq_dict:
    print (prot)
    gene = protein_info_dict[prot][0]
    species, cate, sub_cate = matrix_info_dict[gene]["Species"], matrix_info_dict[gene]["Category"],matrix_info_dict[gene]["Sub"]
    freq_array = global_prot_freq_dict[prot]
    seq_cov = np.count_nonzero(freq_array)/len(matrix_protein_dict[prot])*100
    df.loc[prot,'Gene'] = gene
    df.loc[prot,'Species'] = species
    df.loc[prot, 'Category'] = cate
    df.loc[prot, 'Sub'] = sub_cate
    df.loc[prot, 'Sequence coverage'] = seq_cov
df.to_csv('F:/matrisomedb2.0/statistics/glob_seq_coverage_1.tsv', sep='\t')
"""

# gene numbers each category
"""
df = pd.read_csv('F:/matrisomedb2.0/statistics/glob_seq_coverage.tsv',delimiter='\t',index_col=0)
sub_list = df['Sub'].unique()

# for sub in sub_list:
#     df_slice = df[(df['Species']=='Human')&(df['Sub']==sub)]
#     print (sub,len(df_slice['Gene'].unique()))

# sequence coverage between MD1 and MD2
matrisome_db_1 = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx',index_col=0)
matrisome_db_1['cov_percentage'] = [each*100 for each in matrisome_db_1['cov']]

md1_group = matrisome_db_1.groupby(['category']).mean()
print (md1_group)

md2_group = df.groupby(['Sub']).mean()
print (md2_group)
"""

# PSM number check
"""
# all_psm_list = pickle.load(open('F:/matrisomedb2.0/all_psm.p','rb'))
# print (len(all_psm_list))

base_path = 'F:/matrisomedb2.0/MDB2/result/'
files = glob(base_path+'/**/*.txt',recursive=True)


print (all_matrix_psms(files,matrix_protein_dict))
"""
## PTM check
# glob_id_ptm_map = pickle.load(open('F:/matrisomedb2.0/glob_prot_ptm_ind_dict.p','rb'))

# df = pd.DataFrame(index=[k for k in glob_id_ptm_map.keys()],
#                   columns=['Gene','Species','Category','Sub',
#                            'n[42.0106]','M[15.9949]','R[0.9840]',
#                            'K[15.9949]','P[15.9949]','S[79.9663]',
#                            'T[79.9663]','Y[79.9663]'])
# for prot in glob_id_ptm_map:
#     print (prot)
#     gene = protein_info_dict[prot][0]
#     species, cate, sub_cate = matrix_info_dict[gene]["Species"], matrix_info_dict[gene]["Category"], \
#                               matrix_info_dict[gene]["Sub"]
#     df.loc[prot, 'Gene'] = gene
#     df.loc[prot, 'Species'] = species
#     df.loc[prot, 'Category'] = cate
#     df.loc[prot, 'Sub'] = sub_cate
#     for ptm in glob_id_ptm_map[prot]:
#         df.loc[prot,ptm.replace('\\','')] = len(glob_id_ptm_map[prot][ptm])
# df.to_csv('F:/matrisomedb2.0/statistics/glob_ptms_count.tsv',sep='\t')


## protein ptms average to gene-centric
# df_ptm = pd.read_csv('F:/matrisomedb2.0/statistics/glob_ptms_count.tsv',delimiter='\t',index_col=0)
# df_ptm.groupby(['Sub']).mean().to_csv('F:/matrisomedb2.0/statistics/glob_category_ptms_count.tsv',sep='\t')
# df_ptm_gene = df_ptm.groupby(['Gene']).mean()
# df_ptm_gene.to_csv('F:/matrisomedb2.0/statistics/glob_gene_ptms_count.tsv',sep='\t')

# df_ptm_gene = pd.read_csv('F:/matrisomedb2.0/statistics/glob_gene_ptms_count.tsv',delimiter='\t',index_col=0) # normalize ptm counts by length?
# df_ptm_gene_normalize = pd.DataFrame(index=df_ptm_gene.index,columns=df_ptm_gene.columns)
# for gene in df_ptm_gene.index:
#     ptms = df_ptm_gene.loc[gene,:]/len(matrix_protein_dict[gene_prot_dict[gene]])
#     df_ptm_gene_normalize.loc[gene,:] = ptms
# df_ptm_gene_normalize.to_csv('F:/matrisomedb2.0/statistics/glob_gene_ptms_count_normalize.tsv',sep='\t')

## plot heatmaps to show PTMs

# df_ptm_gene_normalize = pd.read_csv('F:/matrisomedb2.0/statistics/glob_gene_ptms_count_normalize.tsv',delimiter='\t',index_col=0)
# fig,ax = plt.subplots(1,1, figsize=(8,15))

# df_plot = pd.DataFrame(index=df_ptm_gene_normalize.index, columns=df_ptm_gene_normalize.columns)
# for col in df_ptm_gene_normalize:
#     col_max = df_ptm_gene_normalize[col].max()
#     df_plot[col] = df_ptm_gene_normalize[col]/col_max


# g = sns.heatmap(data=df_plot,cbar_kws={'label': 'normalized average PTM counts',
#                                                      'ticks': [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,0.9,1.0],
#                                                      'shrink': 0.3},vmin=0,vmax=1, cmap="viridis",yticklabels=False)
# ax.set_xticklabels(labels=df_plot.columns, rotation=30, fontsize=10, ha='right')
# plt.setp(g.ax_heatmap.get_xticklabels(), rotation=30, ha='right')
# plt.show()

## category ptms
# fig,ax = plt.subplots(1,1, figsize=(10,5))
# df_cat_ptms = pd.read_csv('F:/matrisomedb2.0/statistics/glob_category_ptms_count.tsv',delimiter='\t',index_col=0)
# g = sns.heatmap(data=df_cat_ptms,cbar_kws={'label': 'average PTM counts','shrink': 0.5},vmin=0,vmax=200,
#                 cmap="viridis",yticklabels=True)
# ax.set_xticklabels(labels=df_cat_ptms.columns, rotation=30, fontsize=10, ha='right')
# plt.show()

## seq coverage each sample
"""
data = []
sample_data = pickle.load(open('F:/matrisomedb2.0/sample.data','rb'))
for sample in sample_data:
    id_freq_array = sample_data[sample]['freq']
    sned1_freq = id_freq_array['Q8TER0']
    seq_cov = np.count_nonzero(sned1_freq)/len(matrix_protein_dict['Q8TER0'])*100
    print (sample, seq_cov)
    data.append([sample,seq_cov])

df_sned1 = pd.DataFrame(data,columns=['sample types', 'Q8TER0_sequence_cov'])
df_sned1.to_csv('F:/matrisomedb2.0/Q8TER0_seq_cov.tsv', sep='\t')
"""

## seq covearage histogram
sort_category = ["ECM Glycoproteins","Collagens","Proteoglycans","ECM-affiliated Proteins","ECM Regulators",
                  "Secreted Factors"]
ecm_class_color_dict = {"Collagens": '#0584B7', 'ECM-affiliated Proteins':'#F4651E',
                        'ECM Regulators':"#F9A287","Secreted Factors":"#FFE188",
                        "ECM Glycoproteins":"#13349D", "Proteoglycans":"#59D8E6"}

md1_df = pd.read_excel('D:/data/Naba_deep_matrisome/matrisome coverage_norepeat.xlsx')
md2_df = pd.read_csv('F:/matrisomedb2.0/statistics/glob_seq_coverage.tsv',sep='\t')

# fig,axs = plt.subplots(2,3,figsize=(10,5))
props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
## each category
# for each, ax in zip(sort_category,[[0,0],[0,1],[0,2],[1,0],[1,1],[1,2]]):
#     color = ecm_class_color_dict[each]
#     md1_sub_df = md1_df[md1_df['category']==each]
#     md2_sub_df = md2_df[md2_df['Sub']==each]
#
#     ave_cov_md1 = md1_sub_df['cov'].mean()*100
#     ave_cov_md2 = md2_sub_df['Sequence coverage'].mean()
#     text = 'Ave Seq Cov in MD1: %.2f%%\nAve Seq Cov in MD2: %.2f%%' % (ave_cov_md1, ave_cov_md2)
#
#     axs[ax[0],ax[1]].hist(md1_sub_df['cov']*100,50,color=color, alpha=0.3)
#     axs[ax[0],ax[1]].hist(md2_sub_df['Sequence coverage'], 50, color=color)
#     axs[ax[0],ax[1]].text(0.10, 0.95, text, transform=axs[ax[0],ax[1]].transAxes, fontsize=8,
#         verticalalignment='top', bbox=props)
#     axs[ax[0],ax[1]].set_xlabel(each+' Seq coverage')
#     axs[ax[0],ax[1]].set_ylabel('Frequency')
#     axs[ax[0],ax[1]].set_xlim([-10, 110])
# plt.tight_layout()
# plt.savefig('F:/matrisomedb2.0/statistics/MD1_MD2_ave_cov.png', dpi=300)
# plt.show()

## overall ECM proteins
fig,axs = plt.subplots(figsize=(6,6))
md1_ave_cov = md1_df['cov'].mean()*100
md2_ave_cov = md2_df['Sequence coverage'].mean()
plt.hist(md1_df['cov']*100,50,color='#6d706e',alpha=0.5)
plt.hist(md2_df['Sequence coverage'], 50, color='#111211')
text = 'Ave Seq Cov in MD1: %.2f%%\nAve Seq Cov in MD2: %.2f%%' % (md1_ave_cov, md2_ave_cov)
plt.text(0.10, 0.95, text, transform=axs.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)
plt.ylabel('Frequency')
plt.xlabel('Sequence Coverage in MD1 and MD2')
plt.tight_layout()
plt.savefig('F:/matrisomedb2.0/statistics/all_seq_cov.png', dpi=300)
# plt.show()