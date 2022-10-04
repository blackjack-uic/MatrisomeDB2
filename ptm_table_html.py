"""
output PTM table in html
"""

ptm_map_dict = {'Q\\[129\\]':'Gln deamidation','N\\[115\\]':'ASN deamidation',
                'Q\\[111\\]':'Gln to pyroGln','C\\[143\\]':'selenocysteine',
                'M\\[15\\.9949\\]':'Met oxidation','P\\[15\\.9949\\]':'Pro hydroxylation',
                'K\\[15\\.9949\\]':'Lys hydroxylation','n\\[42\\.0106\\]':'N-term acetylation',
                'C\\[57\\.0215\\]':'Cys redu-alky','R\\[0\\.9840\\]':'Arg deamidation','Y\\[79\\.9663\\]':'Phospho-Tyr',
                'T\\[79\\.9663\\]':'Phospho-Thr', 'S\\[79\\.9663\\]':'Phospho-Ser'}

html_header = '<!DOCTYPE html><html lang="en"><head><meta name="charset" content="UTF-8"><meta name="author" content="Gao lab"><meta name="viewport" content="width=device-width, initial-scale=1"><meta name="description" content="Just testing some table sorting."><title>sortable - the example page</title><link rel="stylesheet" type="text/css" href="sortable.css" /></head><body><script type="text/javascript" src="./js_downloadbutton.js"></script><button onclick="export2csv();">Download CSV</button><table class="sortable"><thead><tr><th>Position</th><th>AA</th><th>PTM</th></tr></thead>'
html_tail = '</tbody><tfoot></tfoot></table><script type="text/javascript" src="sortable.js"></script></body></html>'


def ptm_table_global(id_ptm_idx_dict,output_base):
    """
    output ptm table for global data
    :param id_ptm_idx_dict:
    :param output_base:
    :return:
    """

    for prot in id_ptm_idx_dict:
        html_content = '<tbody>'

        print (prot)
        ptm_ind_dict = id_ptm_idx_dict[prot]
        for ptm in ptm_ind_dict:
            aa = ptm[0]
            ptm_name = ptm_map_dict[ptm]
            for ind in ptm_ind_dict[ptm]:
                html_content+='<tr><td>'+str(ind+1)+'</td><td>'+aa+'</td><td>'+ptm_name+'</td></tr>'
        f = open(output_base+prot+'_ptmtable.html','w')
        f.write(html_header+html_content+html_tail)
        f.close()


if __name__=='__main__':
    import pickle

    glob_ptm_map = pickle.load(open('F:/matrisomedb2.0/glob_prot_ptm_ind_dict.p', 'rb'))
    outbase = 'F:/matrisomedb2.0/ptm_table_revise/'

    ptm_table_global(glob_ptm_map,output_base=outbase)
