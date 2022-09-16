from bokeh.models import HoverTool, ColumnDataSource, FactorRange, LinearColorMapper,ColorBar,BasicTicker,PrintfTickFormatter, Plot, Rect, Legend, LegendItem,SingleIntervalTicker, Label, LabelSet, TableColumn, DataTable, HTMLTemplateFormatter
from bokeh.palettes import Spectral7, Viridis, Plasma, Blues9, Turbo256
from bokeh.transform import factor_cmap
from bokeh.plotting import figure,output_file,save
from bokeh.io import save, output_file, show
from bokeh.embed import components


ptm_map_dict = {'Q\\[129\\]':'Gln deamidation','N\\[115\\]':'ASN deamidation',
                'Q\\[111\\]':'Gln to pyroGln','C\\[143\\]':'selenocysteine',
                'M\\[15\\.9949\\]':'Met oxidation','P\\[15\\.9949\\]':'Pro hydroxylation',
                'K\\[15\\.9949\\]':'Lys hydroxylation','n\\[42\\.0106\\]':'N-term acetylation',
                'C\\[57\\.0215\\]':'Cys redu-alky','R\\[0\\.9840\\]':'Arg deamidation','Y\\[79\\.9663\\]':'Phospho-Tyr',
                'T\\[79\\.9663\\]':'Phospho-Thr', 'S\\[79\\.9663\\]':'Phospho-Ser'}
regex_list = ['M\\[15\\.9949\\]','P\\[15\\.9949\\]','K\\[15\\.9949\\]','n\\[42\\.0106\\]','C\\[57\\.0215\\]',
              'R\\[0\\.9840\\]','T\\[79\\.9663\\]','S\\[79\\.9663\\]','Y\\[79\\.9663\\]']


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

def my_replace(match_obj):
    match_obj = match_obj.group()
    matched_aa = match_obj[0]
    if matched_aa != 'n':
        return matched_aa  # gives back the first element of matched object as string
    else:
        # if first match is n, then n acetylation, get rid of n
        return ''


def peptide_map(psm_dict,protein_dict):
    """
    map peptides to proteome and return freq array dict
    :param psm_dict: {'peptide': [psm1,psm2,...]}
    :return:
    """
    prot_freq_dict = {}
    coverage_dict = {}
    prot_psm_dict = defaultdict(int)
    prot_psm_list_dict = defaultdict(list)

    start = time.time()
    id_list,seq_list = extract_UNID_and_seq(protein_dict)
    seqline = creat_total_seq_line(seq_list,sep='|')
    zeroline = zero_line_for_seq(seqline)
    sep_pos_array = separator_pos(seqline)
    pos_id_dict = read_position_ID_into_dict(id_list,seq_list,zeroline)

    map_start = time.time()
    aho_result = automaton_matching(automaton_trie([pep for pep in psm_dict.keys()]),seqline)
    print (f'aho tree building and mapping took {time.time()-map_start}')
    for pos in aho_result:
        map_pep = pos[2]
        # zeroline[pos[0]:pos[1] + 1] += psm_dict[pos[2]]  # map PSMs instead of peptides
        zeroline[pos[0]:pos[1] + 1] += len(psm_dict[map_pep])
        prot_psm_dict[pos_id_dict[pos[0]]] += len(psm_dict[map_pep])
        prot_psm_list_dict[pos_id_dict[pos[0]]] += psm_dict[map_pep]
    for i in range(len(sep_pos_array) - 1):  # iterate from the left of zeroline all the way to the right
        freq_array = zeroline[sep_pos_array[i] + 1:sep_pos_array[i + 1]]

        prot_freq_dict[id_list[i]] = freq_array
        coverage_dict[id_list[i]] = np.count_nonzero(freq_array)/len(freq_array)*100

    print (f'script took {time.time()-start}s')

    return prot_freq_dict, prot_psm_dict, prot_psm_list_dict, coverage_dict


def ptm_map(psm_list,protein_dict):
    """
    map psms with ptm to protein sequence
    :param psm_list: a list of PSMs with PTM, P[100]EPTIDE
    :param protein_dict:
    :return:
    """
    time_start = time.time()
    regex_pat = '\w{1}\[\d+\.?\d+\]'  # universal ptm pattern

    regex_set = set()
    peptide_psm_dict = defaultdict(list)  # append all psm into a dictionary, {peptide:[psm1,psm2,...]}
    id_ptm_idx_dict = {}  # {protein_id:{ptm1:nonzero_index_array,ptm2:nonzero_index_array,...}}
    id_ptm_freq_dict = {}  # {protein_id:{ptm1:freq_array,ptm2:freq_array,...}}

    for each in psm_list:
        each_reg_sub = re.sub(regex_pat, my_replace, each)
        peptide_psm_dict[each_reg_sub].append(each)
        match = re.findall(regex_pat, each)
        if match:
            for ptm in match:
                regex_set.add(ptm.replace('[','\[').replace(']','\]').replace('.','\.'))
    if 'C\\[57\\.0215\\]' in regex_set:
        regex_set.remove('C\\[57\\.0215\\]')

    # print (regex_set)
    # aho mapping
    id_list, seq_list = extract_UNID_and_seq(protein_dict)
    seq_line = creat_total_seq_line(seq_list, sep="|")
    ptm_index_line_dict = {each:zero_line_for_seq(seq_line) for each in regex_set}
    separtor_pos_array = separator_pos(seq_line)
    aho_result = automaton_matching(automaton_trie([pep for pep in peptide_psm_dict]), seq_line)

    # ptm mapping, n-term mod would give 2 index, need to fix. ---> might fix by add n? at the begining of regex str
    ptm_start = time.time()
    for tp in aho_result:
        matched_pep = tp[2]  # without ptm site
        for psm in peptide_psm_dict[matched_pep]:
            for ptm in regex_set:  # check each ptm, mask other ptms
                new_psm = re.sub('n?\[\d+\.?\d+\]', '', psm.replace(ptm.replace('\\', ''), '*')).\
                    replace('*', ptm.replace('\\', ''))
                ptm_mod = set(re.findall(ptm, new_psm))
                if ptm_mod:
                    for ele in ptm_mod:
                        ### count multiple ptms in a peptide seq
                        num_of_mod = len(
                            re.findall(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm))

                        PTM_index = [m.start() for m in
                                     re.finditer(ele.replace('[', '\[').replace(']', '\]').replace('.', '\.'), new_psm)]
                        PTM_index_clean = [ind - num * (len(ele) - 1) for ind, num in zip(PTM_index, range(num_of_mod))]
                        for indx in PTM_index_clean:
                            ptm_index_line_dict[ptm][tp[0] + indx] += 1
    print (f'ptm mapping took {time.time()-ptm_start}')

    # get ptm freq array and index
    for i in range(len(separtor_pos_array) - 1):
        id_ptm_idx_dict[id_list[i]] = {ptm:
            np.nonzero(ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]])[0]
            for ptm in regex_set}
        id_ptm_freq_dict[id_list[i]] = {ptm:
            ptm_index_line_dict[ptm][separtor_pos_array[i] + 1:separtor_pos_array[i + 1]]
            for ptm in regex_set}
    print (f'ptm script took {time.time()-time_start}s')
    return id_ptm_idx_dict, id_ptm_freq_dict


def seq_cov_gen(freq_array,ptm_dict,protein_seq):
    """
    generate seq coverage map in html format (inside body)
    :param freq_array: mapped 1d np zero array for one protein
    :param protein_seq:
    :param protein_id:
    :param ptm_dict: {ptm1:[index list],ptm2:[index list]}
    :param gene_name:
    :return:
    """
    # marks for html
    """mark {
	background-color: #2C3E50;
	color: white;
}
mark1 {
	background-color: #2471A3;
	color: white;
}
mark2 {
	background-color: #AED6F1;
	color: white;
}
mark3 {
	background-color: #95A5A6;
	color: white;
}
mark4 {
	background-color: #D7DBDD;
	color: white;
}
mark5 {
	background-color: white;
	color: red;
}
    """
    print (len(freq_array),len(protein_seq))
    seq_cov = np.count_nonzero(freq_array)/len(freq_array)*100
    split_seq = np.arange(0, len(protein_seq), 165)
    split_seq = np.append(split_seq,len(protein_seq))
    max_freq = np.max(freq_array)
    ptm_freq_array = set([idx for ptm in ptm_dict for idx in ptm_dict[ptm]])
    # print (ptm_freq_array)

    seq_cov_str = ''
    for i in range(len(split_seq) - 1):
        for j in range(split_seq[i], split_seq[i + 1]):
            if j not in ptm_freq_array:
                if freq_array[j] == 0:
                    seq_cov_str += protein_seq[j]
                elif 1 <= freq_array[j] < 0.2 * max_freq:
                    seq_cov_str += '<mark4>' + protein_seq[j] + '</mark4>'
                elif 0.2 * max_freq <= freq_array[j] < 0.4 * max_freq:
                    seq_cov_str += '<mark3>' + protein_seq[j] + '</mark3>'
                elif 0.4 * max_freq <= freq_array[j] < 0.6 * max_freq:  # color legend changeable
                    seq_cov_str += '<mark2>' + protein_seq[j] + '</mark2>'
                elif 0.6 * max_freq <= freq_array[j] < 0.8 * max_freq:  # color legend changeable
                    seq_cov_str += '<mark1>' + protein_seq[j] + '</mark1>'
                else:
                     seq_cov_str += '<mark>' + protein_seq[j] + '</mark>'
            else:
                seq_cov_str += '<mark5>' + protein_seq[j] + '</mark5>'
        seq_cov_str += '\n'
    return seq_cov_str,seq_cov


def hashcolor(s):

    return Turbo256[hash(s) % 256]


def domain_cov_ptm(prot_freq_dict, ptm_map_result, domain_pos_dict,protein_entry:str, data_source='sample'):
    """
    -----
    draw rectangular box as protein domains and alpha as coverage on bokeh,
    and label PTMs.
    -----
    :param protein_freq_dict:
    :param domain_pos_dict:
    :param protein_entry:
    :return:
    """
    if protein_entry not in domain_pos_dict:
        return '','No SMART domain available'
    else:
        time_start = time.time()
        freq_array = prot_freq_dict[protein_entry]
        domain_dict = domain_pos_dict[protein_entry]
        protein_len = len(freq_array)

        ## coverage every 10 aa
        pos_cov_dict = {}
        bin_width = 10
        bar_shrink_raio = 5 # make bar shorter
        bar_bottom = 0.8
        interval=np.arange(0,protein_len,bin_width)
        for i in interval[:-1]:
            coverage = np.count_nonzero(freq_array[i:i+bin_width])/bin_width
            pos_cov_dict[i+bin_width/2] = coverage/bar_shrink_raio + bar_bottom  # move bar up
        pos_cov_dict[interval[-1]+bin_width/2] = np.count_nonzero(freq_array[interval[-1]:protein_len])/(protein_len-interval[-1])/bar_shrink_raio+bar_bottom
        source_cov = ColumnDataSource(dict(x=[k for k in pos_cov_dict.keys()],y=[v for v in pos_cov_dict.values()],
                                           label=['{:.1f}%'.format((v-0.8)*100) for v in pos_cov_dict.values()]))

        ## extract domain position and calculate domain coverage
        info_list = []

        for each_domain in domain_dict:
            for each_tp in domain_dict[each_domain]:
                start, end = each_tp[0], each_tp[1]
                coverage = np.count_nonzero(freq_array[start - 1:end]) / (end - start + 1)
                info_list.append((start,end,coverage,each_domain))
        if info_list ==[]:
            return '', 'No SMART domain available'
        else:
            start, end, coverage, domain_list = zip(*info_list)


            # x coordinates and widths of rectangular
            x,width = zip(*[((end_-start_)/2+start_,end_-start_) for end_,start_ in zip(end,start)])
            # hash color to show each domain
            color_map_dict = {domain:hashcolor(domain) for domain in set(domain_list)}
            color_list = [color_map_dict[each] for each in domain_list]

            source = ColumnDataSource(dict(x=x,w=width,color=color_list,domain=domain_list,
                                           start=start,end=end,
                                           coverage=coverage))

            ## prepare data for PTM labeling

            ptm_index_dict = ptm_map_result[protein_entry]
            ptm_index_sort = sorted([(idx,each) for each in ptm_index_dict for idx in ptm_index_dict[each]])
            ptm_count = sum([len(ptm_index_dict[ptm]) for ptm in ptm_index_dict])
            # bokeh plot, hovertool
            hover = HoverTool(names=['rec'],tooltips=[('domain', '@domain'), ('start position', '@start'),('end position','@end'),('domain coverage','@coverage{:.1%}'),])
            # initiate bokeh figure
            p = figure(x_range=(-10,protein_len),
                       y_range=(0,2),
                       tools=['pan', 'box_zoom', 'wheel_zoom', 'save',
                              'reset', hover],
                       plot_height=500, plot_width=700,
                       toolbar_location='right',
                       title='',
                       x_axis_label='amino acid position')

            # plot domains as rectangular and alpha shows coverage
            p.rect(x="x", y=0.6, width='w', height=50,
                   source=source,
                   fill_color='color',
            #       fill_alpha='coverage',
                   line_width=2,
                   line_color='black',
                   height_units="screen",
                   name='rec'
                   )
            # reference domains, alpha=1
            # p.rect(x="x", y=1, width='w', height=10,
            #        source=source,
            #        fill_color='color',
            #        line_width=2,
            #        line_color='black',
            #        height_units="screen",
            #        name='rec'
            #        )

            # sequence coverage bar charts
            p.vbar(x='x',width=bin_width,top='y',bottom=bar_bottom,source=source_cov,color='#D3D3D3', name='seq cov')
            # cov_label = LabelSet(x='x',y='y',text='label',text_font_size='8px',
            #                      x_offset=-13.5, y_offset=0, source=source_cov)
            # p.add_layout(cov_label)
            cov_bar_legend_top, cov_bar_legend_bottom = bar_bottom+1/bar_shrink_raio,bar_bottom
            cov_bar_x_coor = -8
            # p.vbar(x=[cov_bar_x_coor],width=bin_width,top=[cov_bar_legend_top],bottom=[cov_bar_legend_bottom],
            #        color='#D3D3D3',name='seq_cov_legend')
            # label annotations
            for y_coor, text in zip([cov_bar_legend_bottom,cov_bar_legend_top],['0%','100%']):
                label_cov = Label(x=cov_bar_x_coor,y=y_coor,x_offset=0, y_offset=-5,text=text,text_font_size='10px',text_align='left')
                p.add_layout(label_cov)
            seq_cov_title = Label(x=-5, y=1.1, text='Sequence coverage binned by every 10 aa',text_font_size='12px',text_align='left',text_color='#A9A9A9')
            p.add_layout(seq_cov_title)

            # line shows whole protein length
            p.line(x=[0,protein_len],y=[0.6,0.6],line_width=10,color='#000000',alpha=0.8,name='line')

            if data_source=='sample': # if plotting from sample-specific data, label PTMs below domains
                # adjusted PTM text coordinates calculation
                numpy_zero_array = np.zeros((1500, protein_len)) # mask numpy array for text plotting
                ptm_x,ptm_y = [], []
                new_ptm_x, new_ptm_y = [], [] # adjusted text coordinates to prevent overlap
                ptms = []
                for tp in ptm_index_sort:
                    each_idx, ptm = tp
                    # ptms.append(ptm.replace('\\', ''))
                    ptms.append(ptm)
                    ptm_x.append(each_idx)
                    ptm_y.append(0.5)
                    x_offset,y_offset = 0,0
                    x_move = int(protein_len/8)  # 12 to 8 when plot width is 1200 to 700
                    while True: # keep moving down if text are too close
                        nonzero_count = np.count_nonzero(numpy_zero_array[1490+y_offset:1500+y_offset,each_idx+x_offset:each_idx+x_move+x_offset])
                        if nonzero_count == 0:
                            # print (ptm,each_idx,x_offset,y_offset)
                            new_ptm_x.append(each_idx+x_offset)
                            new_ptm_y.append((25+y_offset)/200*2)
                            numpy_zero_array[1490+y_offset:1500+y_offset,each_idx+x_offset:each_idx+x_move+x_offset] += 1
                            break
                        else:
                            # print ('moving down')
                            # x_offset += 50 # value to move right
                            y_offset -= 12 # value to move down

                # label ptm and connect to protein domains
                for x,y,x_,y_,ptm in zip(new_ptm_x,new_ptm_y,ptm_x,ptm_y,ptms):
                    p.line(x=[x_+1,x+1],y=[y_,y+0.1],line_width=1,color='black',alpha=0.3) # connect domain with text
                    label = Label(x=x,y=y,text=ptm_map_dict[ptm]+'\n'+str(x_+1),text_font_size='10px', text_align='center', text_font='Tahoma')
                    p.add_layout(label)
            else: # if from global data, do not label PTMs
                label = Label(x=1,y=0.10,text=str(ptm_count)+' PTMs occurrence in total\nSee details from PTMs table', text_font_size = '15px', text_font = 'Tahoma')
                p.add_layout(label)

            # dummy glyphs to help draw legend
            legend_gly = [p.line(x=[1, 1], y=[1, 1], line_width=15, color=c, name='dummy_for_legend')
                            for c in [v for v in color_map_dict.values()]]

            legend = Legend(title='Domains', background_fill_color='white',
                            border_line_color='black',border_line_width=3,
                            border_line_alpha=0.7,
                            items=[LegendItem(label=lab, renderers=[gly])
                                   for lab, gly in zip([d for d in color_map_dict.keys()],legend_gly)])
            # alpha color bar, domain coverage
            # color_mapper = LinearColorMapper(palette=Blues9[::-1], low=0, high=1)
            # color_bar = ColorBar(color_mapper=color_mapper,location=(0, 0),ticker=SingleIntervalTicker(interval=0.1))
            # p.add_layout(color_bar,'right')

            p.add_layout(legend)
            p.xgrid.visible = False
            p.ygrid.visible = False
            p.yaxis.visible = False
            print (f'bokeh graph took {time.time()-time_start}s')
            # show(p)
            return components(p)


def ptm_table_bokeh3(sample_data, protein_dict, output_base):
    """
    show ptm positions for each protein
    :param id_ptm_idx_dict:
    :param protein_dict:
    :param uniprot_gene_dict:
    :return:
    """
    import pandas as pd

    for sample in sample_data:


        ptm_id_index_dict = sample_data[sample]['ptm']

        for prot in ptm_id_index_dict:
            print (prot)
            # seq = protein_dict[prot]
            info_dict = {}
            # info_dict['Amino acid sequence'] = [aa for aa in seq]
            # info_dict['Amino acid position'] = range(1,len(seq)+1)

            for ptm in ptm_id_index_dict[prot]:

                info_dict[ptm_map_dict[ptm]] = [each+1 for each in ptm_id_index_dict[prot][ptm]]


            df = pd.DataFrame(dict(PTMs=[k for k in info_dict.keys()],Positions=[v for v in info_dict.values()]))
            # print (df)
            source = ColumnDataSource(df)

            columns = [TableColumn(field=each,title=each)
                        for each in df.columns]
            table = DataTable(source=source,columns=columns, width=1000, height=400, editable=True)
            output_file(output_base+sample.replace('/','-')+'_'+prot+'_ptmtable.html')
            save(table)


def psmlist_todict(psm_list):
    regex_pat = '\w{1}\[\d+\.?\d+\]'
    psm_dict = defaultdict(list)

    for psm in psm_list:
        reg_sub = re.sub(regex_pat, my_replace, psm)
        psm_dict[reg_sub].append(psm)
    return psm_dict


if __name__ == '__main__':
    import json
    import pickle
    import os
    with open('F:/matrisomedb2.0/smart_domain_0908.json') as f_o:
        info_dict = json.load(f_o)

    protein_dict = fasta_reader('F:/matrisomedb2.0/mat.fasta')

    protein_info_dict = protein_info_from_fasta('F:/matrisomedb2.0/mat.fasta')

    glob_prot_freq_dict = pickle.load(open('F:/matrisomedb2.0/glob_prot_freq_dict.p', 'rb'))
    glob_ptm_map = pickle.load(open('F:/matrisomedb2.0/glob_prot_ptm_ind_dict.p', 'rb'))
    html_tempalte = open(r'F:\matrisomedb2.0\test/domain_seq_cov_html_template_0909.html', 'r')
    html_tempalte_read = html_tempalte.read()

    # one protein test
    start = time.time()
    out_put = r'F:\matrisomedb2.0/test/test2.html'
    sample = 'Pancreatic Ductal Adenocarcinoma Xenograft (BxPC3)'
    prot = 'P21980'
    smart_url = 'https://smart.embl.de/smart/show_motifs.pl?ID=' + prot
    sample_psm_list = pickle.load(open(r'F:\matrisomedb2.0/data_for_test/sample_psm_list.p', 'rb'))
    sample_psm_dict = psmlist_todict(sample_psm_list)
    # protein_dict = {prot:protein_dict[prot],'XXX':'XXX','AAA':'AAA'}
    # global_psm_list = global_protein_psm_dict[prot]
    # glob_psm_dict = psmlist_todict(global_psm_list)

    sample_prot_freq_dict = peptide_map(sample_psm_dict, protein_dict)[0]
    # print (prot,len(sample_prot_freq_dict[prot]),len(protein_dict[prot]))
    # glob_prot_freq_dict = peptide_map(glob_psm_dict,protein_dict)[0]

    sample_ptm_map = ptm_map(sample_psm_list, protein_dict)[0]
    # glob_ptm_map = ptm_map(global_psm_list,protein_dict)[0]

    sample_seq_cov = seq_cov_gen(sample_prot_freq_dict[prot], sample_ptm_map[prot], protein_dict[prot])
    glob_seq_cov = seq_cov_gen(glob_prot_freq_dict[prot], glob_ptm_map[prot], protein_dict[prot])

    sample_domain_cov = domain_cov_ptm(sample_prot_freq_dict, sample_ptm_map, info_dict, prot, data_source='sample')
    glob_domain_cov = domain_cov_ptm(glob_prot_freq_dict, glob_ptm_map, info_dict, prot, data_source='global')

    new_html = html_tempalte_read.replace('<!-- COPY/PASTE domain coverage SCRIPT HERE -->', sample_domain_cov[0]). \
        replace('<!-- COPY/PASTE domain coverage global SCRIPT HERE -->', glob_domain_cov[0]). \
        replace('<!-- INSERT domain DIVS HERE -->', sample_domain_cov[1]). \
        replace('<!-- INSERT domain DIVS global HERE -->', glob_domain_cov[1]). \
        replace('<!-- COPY/PASTE seq coverage value HERE-->', str(sample_seq_cov[1])). \
        replace('<!-- COPY/PASTE seq coverage value global HERE-->', str(glob_seq_cov[1])). \
        replace('<!-- COPY/PASTE seq coverage str HERE-->', sample_seq_cov[0]). \
        replace('<!-- COPY/PASTE seq coverage str global HERE-->', glob_seq_cov[0]). \
        replace('<!-- UniprotID -->', protein_info_dict['P21980'][0] + '  (' + protein_info_dict['P21980'][
        1].rstrip(' ') + ')').replace('<!-- sample_type -->', sample). \
        replace('<!-- 3d cov URL -->', sample + '_' + prot + '_3dcov.html').replace('<!-- 3d cov global URL -->',
                                                                                    prot + '_3dcov.html'). \
        replace('<!-- PTM table URL -->', sample + '_' + prot + '_ptmtable.html').replace(
        '<!-- PTM table global URL -->', prot + '_ptmtable.html').replace('<!-- SMART URL -->', smart_url)
    with open(out_put, 'w') as f_o:
        f_o.write(new_html)
    print(f'time for {prot}: {time.time() - start}')