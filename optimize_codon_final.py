from Bio import SeqIO
# from Bio.Data import CodonTable
from collections import Counter
"""
This is a python program which is to calculate the codes usage of a species including termination codon 
I used the fasta file of genome and protein whose id was used to regard as cds_dict and pron_dict
I need remember how to cycle the id of genome and protein's file, which can ensure more information and accuracy.
"""
codon_table = [["TTT", "F"], ["TTC", "F"], ["TTA", "L"], ["TTG", "L"],
          ["TCT", "S"], ["TCC", "S"], ["TCA", "S"], ["TCG", "S"],
          ["TAT", "Y"], ["TAC", "Y"], ["TGT", "C"], ["TGC", "C"],
          ["TGG", "W"], ["CTT", "L"], ["CTC", "L"], ["CTA", "L"],
          ["CTG", "L"], ["CCT", "P"], ["CCC", "P"], ["CCA", "P"],
          ["CCG", "P"], ["CAT", "H"], ["CAC", "H"], ["CAA", "Q"],
          ["CAG", "Q"], ["CGT", "R"], ["CGC", "R"], ["CGA", "R"],
          ["CGG", "R"], ["ATT", "I"], ["ATC", "I"], ["ATA", "I"],
          ["ATG", "M"], ["ACT", "T"], ["ACC", "T"], ["ACA", "T"],
          ["ACG", "T"], ["AAT", "N"], ["AAC", "N"], ["AAA", "K"],
          ["AAG", "K"], ["AGT", "S"], ["AGC", "S"], ["AGA", "R"],
          ["AGG", "R"], ["GTT", "V"], ["GTC", "V"], ["GTA", "V"],
          ["GTG", "V"], ["GCT", "A"], ["GCC", "A"], ["GCA", "A"],
          ["GCG", "A"], ["GAT", "D"], ["GAC", "D"], ["GAA", "E"],
          ["GAG", "E"], ["GGT", "G"], ["GGC", "G"], ["GGA", "G"],
          ["GGG", "G"], ["TGA", "&"], ["TAG", "&"], ["TAA", "&"]]

def getDict(inputfile):
    dict_name = {}
    for record in SeqIO.parse(inputfile, "fasta"):
        dict_name[record.id] = str(record.seq)
    return dict_name

def getCdsPronCount(cds_dict,pron_dict):
    cds_stop = []
    table_list = []
    getsameid_list = [gene_id for gene_id in cds_dict if gene_id in pron_dict]
    for gene_id in getsameid_list:
        cds_tmp = cds_dict[gene_id]
        pro_tmp = pron_dict[gene_id]
        step = 3
        # 将字符串三个一组切分
        split_cds = [cds_tmp[i:i + step] for i in range(0, len(cds_tmp), step)]
        pron_string = pro_tmp.replace('*', '')
        split_pron = []
        for arr in pron_string:
            split_pron.append(arr)
        cds_pron = [[a, b] for a, b in zip(split_cds, split_pron)]
        # 将上面匹配得到的列表cds_pron与标准密码子表table1比对得到相同的结果并放到列表same中，x代表每个匹配的内容比如： ["ATG", "M"]
        same = [x for x in cds_pron if x in codon_table]
        same_getcds = []  # 存放密码子
        for arr in same:
            same_getcds.append(arr[0])  # 每一项只取第一位比如：["ATG", "M"]=>ATG
        # 得到每个密码子出现次数
        count =Counter(same_getcds)
        table_list.append(count)
    # 获取终止密码子
    get_id_list = [gene_id for gene_id in cds_dict]
    for gene_id in get_id_list:
        cds_seq = cds_dict[gene_id]
        cds_seq = cds_seq[-3:]
        if cds_seq == 'TAA':
            cds_stop.append(cds_seq)
        elif cds_seq == 'TGA':
            cds_stop.append(cds_seq)
        elif cds_seq == 'TAG':
            cds_stop.append(cds_seq)
    table_list.append(Counter(cds_stop))
    result = sum(map(Counter, table_list), Counter())
    #获取标准密码子表转化字典
    # standard_table = CodonTable.unambiguous_dna_by_name['Standard'].forward_table
    tmpdict = {} #存放输出的密码子和氨基酸
    for standard_table_key in codon_table:
        for result_cds in dict(result).keys():
            if standard_table_key[0] == result_cds:
                tmpdict[result_cds] = standard_table_key[1]
    csv_dict = {}
    for k, v in tmpdict.items():
        csv_dict[(v,result[k])] = k
    return  csv_dict

def read_cds_split(txtfile):
    with open(txtfile,'r') as cds_data:
        tmp_string = ''
        for single_cds in cds_data:
            single_cds = single_cds.strip()
            tmp_string += single_cds
        split_cds = [tmp_string[i:i + 3] for i in range(0, len(tmp_string), 3)]
        return split_cds


def get_cdspronmax(csv_dict):
    sort_pron_count = sorted(csv_dict.keys(), key=lambda x: (x[0], x[1]))
    tmpdict_pron_count = {k: v for k, v in sort_pron_count}
    max_pron_count = []
    for keys in tmpdict_pron_count:
        max_pron_count.append((keys, tmpdict_pron_count[keys]))
    # print(max_pron_count)

    result_list = []
    for pron_count in csv_dict.keys():
        for max in max_pron_count:
            if pron_count == max:
                result = [csv_dict[pron_count], max[0]]
                result_list.append(result)
    # print(result_list)
    return result_list

def replace_cds(cds_list,optimize_bestusage_codon):

    cds_pron = [] #将cds_list的密码子对应相应的氨基酸然后存放cds_pron中
    for tmp_cds in cds_list:
        for standard_table in codon_table:
            if tmp_cds == standard_table[0]:
                tmp_cds = standard_table
                cds_pron.append(tmp_cds)

    result_list = [] #存放密码子和氨基酸，密码子替换成使用率最高的
    for list_cds in cds_pron:
        for table_key in optimize_bestusage_codon:
            if table_key[1] == list_cds[1]:
                list_cds[0] = table_key[0]
                result_list.append(list_cds)

    result_string = ''
    for get_cds in result_list:
        result_string += get_cds[0]
    return result_string


if __name__ == '__main__':
    cdsfile = "cds.fasta"
    pronfile = "protein.fasta"
    txtfile = "HPH Gene.txt"
    output = "result.fasta"

    cds_dict = getDict(cdsfile)
    pron_dict = getDict(pronfile)
    csv_dict = getCdsPronCount(cds_dict, pron_dict)
    cds_list = read_cds_split(txtfile)
    optimize_bestusage_codon = get_cdspronmax(csv_dict)
    result_fasta = replace_cds(cds_list, optimize_bestusage_codon)
    with open(output, 'w') as outputfile:
        outputfile.write(str(result_fasta))


