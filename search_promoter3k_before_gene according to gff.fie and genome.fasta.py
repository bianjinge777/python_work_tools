from BCBio import GFF
from Bio import SeqIO
out = open('promoter_out.fasta', 'w')

def getSequenceFragment(gff_file, gene_id):
    candi_gene = ''
    with open("D:\pylib\study\A188.gff", 'r') as in_handle:
        OK_flag = 0
        for rec in GFF.parse(in_handle):
            # print(rec)
            if OK_flag == 1:
                break
            for gene in rec.features:
                # print(gene)
                if gene.type == 'gene':
                    key_value = ''
                    if 'Name' in gene.qualifiers:
                        key_value = 'Name'
                    elif 'ID' in gene.qualifiers:
                        key_value = 'ID'
                    if gene_id == gene.qualifiers[key_value][0]:
                        candi_gene = gene
                        OK_flag = 1
                        break
    if candi_gene == '':
        ValueError("candi gene not found")

    candi_mRNA = ['mRNA']

    abs_coord = [] #存放起始位置、终止位置、基因id
    for feature_tmp in candi_gene.sub_features:
        if feature_tmp.type in candi_mRNA:   # 判断类型是否是mRNA，然后得到起始位置、终止位置
            start_tmp = feature_tmp.location.start + 1
            end_tmp = feature_tmp.location.end
            abs_coord.append((start_tmp, end_tmp, feature_tmp.id))

    strand = gene.strand
    site_tmp = [] #创建临时字典存放起始位置、终止位置
    for each_point in abs_coord:
        site_tmp.extend([each_point[0], each_point[1]])

    rel_coord = []      # rel_coord存放片段的位置信息 # 如果是正链(目标位置，起始位置)，如果是负链(终止位置,目标位置)
    for each_point in abs_coord:
        length = 3100
        if strand == 1:
            base_point = min(site_tmp) #正链取左端位置即小的值
            get_site = [abs(base_point - length), abs(base_point)]
        elif strand == -1:  #负链取右端位置即大的值
            base_point = max(site_tmp)
            get_site = [abs(base_point), abs(base_point + length)]
        # rel_coord.append((get_site, i[2]))
        rel_coord.append(get_site)
    rel_coord = sorted(rel_coord, key=lambda x: x[1])
    gene_name = gene_id + '-t1'  #得到mRNA的id
    all_dict = {} #创建字典存放fastaid和序列
    def getDict(inputfile, dict_name):
        for record in SeqIO.parse(inputfile, "fasta"):
            dict_name[record.id] = str(record.seq)
    getDict("D:\pylib\study\A188.fasta", all_dict)
    for id in all_dict.keys():
        if id == rec.id:
            result_seq = all_dict[id][rel_coord[0][0]:rel_coord[0][1]]  # 如果对应id相同然后得到序列信息，最后截取片段
            print(result_seq)
            out.write(str(gene_name)+':'+str(result_seq)+'\n')
# getSequenceFragment('D:\pylib\study\A188.gff', 'A188G15514')
with open('D:\pylib\study\geneid.txt','r') as geneid:  # geneid.txt是用到的id，可以循环读取得到每个id对应的序列片段
    for gene_id in geneid:
        gene_id = gene_id.strip()
        getSequenceFragment('D:\pylib\study\A188.gff', gene_id)

