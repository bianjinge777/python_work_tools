from Bio import SeqIO
from Bio.Data import CodonTable
from collections import Counter
import pandas as pd
cds_dict = {}
pron_dict = {}
table_list = []
"""
This is a python program which is to calculate the codes usage of a species
I used the fasta file of genome and protein whose id was used to regard as cds_dict and pron_dict
I need remember how to cycle the id of genome and protein's file, which can ensure more information and accuracy.
"""
def getDict(inputfile,dict_name):
    for record in SeqIO.parse(inputfile, "fasta"):
        dict_name[record.id] = str(record.seq)
getDict("D:\pylib\study\cds.fa", cds_dict)
getDict("D:\pylib\study\protein.fa", pron_dict)
getsameid_list = [gene_id for gene_id in cds_dict if gene_id in  pron_dict]
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
    table1 = [["TTT", "F"], ["TTC", "F"], ["TTA", "L"], ["TTG", "L"],
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
              ["GGG", "G"]]

    # 将上面匹配得到的列表cds_pron与标准密码子表table1比对得到相同的结果并放到列表same中，x代表每个匹配的内容比如： ["ATG", "M"]
    same = [x for x in cds_pron if x in table1]
    same_getcds = []  # 存放密码子
    for arr in same:
        same_getcds.append(arr[0])  # 每一项只取第一位比如：["ATG", "M"]=>ATG
    # 得到每个密码子出现次数
    count =Counter(same_getcds)
    table_list.append(count)
result = sum(map(Counter, table_list), Counter())
#获取标准密码子表转化字典
standard_table = CodonTable.unambiguous_dna_by_name['Standard'].forward_table
tmpdict = {} #存放输出的密码子和氨基酸
for standard_table_key in standard_table.keys():
    for result_cds in dict(result).keys():
        if standard_table_key == result_cds:
            tmpdict[result_cds] = standard_table[standard_table_key]
result_list = [] #存放循环得到的每个密码子对应的氨基酸和次数
for k, v in tmpdict.items():
    single_dict = {'密码子':k,'氨基酸':v,'次数':result[k]}
    result_list.append(single_dict)
printfile = pd.DataFrame(result_list)
printfile.to_csv('out.csv')
