
def getdictfromtext(text_file):
    import re
    from collections import defaultdict
    word_freq_dict = defaultdict(int)
    stop_pattern = '^a$|^A$|^the$|this|these|those|one|an|The|to|^in$|for|of|or|by|with|is|on|that|be|from|here|^can$|^cannot$|\(|\)|yet|and|Objective|Results|Conclusions|Methods|Background|Conclusions|Approach'

    with open(text_file,'r',encoding='utf8') as f_o:
        for line in f_o:
            line = line.replace('\n','').replace('.','').replace(',','').replace(':','').lstrip('').rstrip('')
            line_split = line.split(" ")
            for word in line_split:
                if re.match(stop_pattern,word):
                    continue
                else:
                    word_freq_dict[word] += 1

    return word_freq_dict


def get_freq_from_str(list_str:list):
    from collections import defaultdict
    import re
    # multiple strings
    word_freq_dict = defaultdict(int)
    stop_pattern = 'as|if|long|much|soon|though|because|before|by the time|even if|in case|in order that|in the event that|least|only|only if|provided|' \
                   'that|once|after|After|we|are|We|Were|as|was|which|in|changes|change|identified|identify|using|Both|both|^a$|^A$|^the$|^\d{1}|this|This|' \
                   'these|These|those|Those|one|an|The|to|To|^in$|In|for|For|of|or|by|with|is|on|On|that|be|from|here|^can$|^cannot$|\(|\)' \
                   '|yet|and|Objective|Results|results|result|Conclusions|Methods|methods|method|Background|Conclusions|Approach|have|Have|Had|Has|has|their|' \
                   'despite|Here|including|its|also|revealed|than|other|role|well|but|may|show|at|showed|compared|until|' \
                   'performed|biology|study|nor|By|Though|though|Furthermore|challenging|although|since|Since|supposing|till|' \
                   'when|whenever|where|whereas|wherever|whether or not|while|unless|group|Group|further|Further|rather|finally|Finally|studies'

    for each_str in list_str:
        each_str_clean = each_str.replace('\n','').replace('.','').replace(',','').replace(':','').lstrip('').rstrip('')
        str_split = each_str_clean.split(" ")
        for word in str_split:
            if re.match(stop_pattern, word):
                continue
            else:
                if len(word)>1:
                    word = re.sub(r'\>|\<|\(|\)|\"','',word)
                    word_freq_dict[word] += 1

    return word_freq_dict


def word_cloud_enrich(total_text_file, single_text_file):
    """
    -----
    plot enrichment word cloud
    -----
    :param total_text_file: all files into one text file
    :param single_text_file: individual text file
    :return: enrichment word frequency dictionary
    """
    enrich_freq_dict = {}

    total_wordfreq_dict = getdictfromtext(total_text_file)
    single_wordfreq_dict = getdictfromtext(single_text_file)
    print (total_wordfreq_dict,single_wordfreq_dict)

    sum_total = sum([v for v in total_wordfreq_dict.values()])
    normalized_total = {w:total_wordfreq_dict[w]/sum_total for w in total_wordfreq_dict}

    sum_single = sum([v for v in single_wordfreq_dict.values()])
    normalized_single = {w:single_wordfreq_dict[w]/sum_single for w in single_wordfreq_dict}

    for w in normalized_single:
        # if w in normalized_total:
        enrich_freq_dict[w] = (normalized_single[w]-normalized_total[w])/normalized_total[w]*100
        # else:
        # enrich_freq_dict[w] = normalized_single[w]*100
    return enrich_freq_dict


def textcloud_from_freq(freq_dict,output_png=None):
    wc = WordCloud(background_color="white",width=1600, height=800).generate_from_frequencies(freq_dict)
    plt.figure(figsize=(20, 10))
    plt.imshow(wc, interpolation='bilinear')
    plt.axis('off')
    plt.tight_layout(pad=0)
    # plt.show()
    if output_png:
        # wc.to_file(output_png)
        plt.savefig(output_png)

def get_abstract(rep_id_list:list,df):

    df_slice = df[df.PXD.isin(rep_id_list)]
    projectid_abstract_dict = {projectid:abstract for projectid, abstract in zip(df_slice['PXD'],df_slice['Abstract'])}
    str_list_subset = [v for v in projectid_abstract_dict.values()]
    word_freq_subset = get_freq_from_str(str_list_subset)

    str_list_total = df['Abstract'].tolist()
    word_freq_total = get_freq_from_str(str_list_total)

    return word_freq_subset, word_freq_total


def gen_wordcloud_md2(rep_id_list:list, df, output_png='test.png'):
    enrich_freq_dict = {}

    word_freq_subset, word_freq_total = get_abstract(rep_id_list,df)

    sum_total = sum([v for v in word_freq_total.values()])
    normalized_total = {w: word_freq_total[w] / sum_total for w in word_freq_total}

    sum_single = sum([v for v in word_freq_subset.values()])
    normalized_single = {w: word_freq_subset[w] / sum_single for w in word_freq_subset} # delete word that only shows once

    for w in normalized_single:
        # if w in normalized_total:
        enrich_freq_dict[w] = (normalized_single[w] - normalized_total[w]) / normalized_total[w] * 100
    print ('enrichment score:',sorted(enrich_freq_dict.items(), key=lambda x: x[1],reverse=True))
    textcloud_from_freq(enrich_freq_dict,output_png=output_png)
    return output_png


if __name__=='__main__':
    import matplotlib.pyplot as plt
    from wordcloud import WordCloud
    import pandas as pd

    # matrisomeDB2.0

    df = pd.read_excel(r'F:\matrisomedb2.0/Abstracts for Word Cloud.xlsx')
    project_ids = ['PXD005130','MSV000082639','PXD020187']

    gen_wordcloud_md2(project_ids,df,output_png='F:/matrisomedb2.0/statistics/3projects_0914_1.png')