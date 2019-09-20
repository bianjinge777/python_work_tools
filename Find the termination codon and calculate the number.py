from Bio import SeqIO
from collections import Counter
import pandas as pd

cds_dict = {}
# all_cds = []
cds_stop = []
"""
a program to look for the number of termination codons.
"""

def getDict(inputfile, dict_name):  #获取id和序列存到字典
    for record in SeqIO.parse(inputfile, "fasta"):
        dict_name[record.id]=str(record.seq)
getDict("D:\pylib\study\cds.fa", cds_dict)

get_id_list = [gene_id for gene_id in cds_dict]
for gene_id in get_id_list:
    cds_seq = cds_dict[gene_id]
        # print(cds_seq)
    cds_seq = cds_seq[-3:]
    # cds_seq = cds_seq[-4:-1]
    # print(cds_seq)
    if cds_seq =='TAA':
        cds_stop.append(cds_seq)
    elif cds_seq == 'TGA':
        cds_stop.append(cds_seq)
    elif cds_seq == 'TAG':
        cds_stop.append(cds_seq)
    # elif cds_seq != 'TAA, TGA,TAG':
    #     cds_stop.append(cds_seq)
printfile = pd.DataFrame(columns=['密码子', '次数'], data=Counter(cds_stop).most_common())
print(printfile)
