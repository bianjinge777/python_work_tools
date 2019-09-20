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
        raise ValueError("candi gene not found")

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
        
        
####################################################################################################
from BCBio import GFF
from Bio import SeqIO

def read_list_file(input_file):
    with open(input_file,'r') as geneid:
    output_list = []
    for gene_id in geneid:
        gene_id = gene_id.strip()
        output_list.append(gene_id)
    return output_list
    
def gff_reader(gff_file):
    gff_dict = {}
    with open(gff_file, 'r') as in_handle:
        for rec in GFF.parse(in_handle):
            for gene in rec.features:
                if gene.type == 'gene':
                    gff_dict[gene.id] = (gene,rec.id)
    return gff_dict

def read_fasta(fasta_file):
    dict_name = {}
    for record in SeqIO.parse(fasta_file, "fasta"):
        dict_name[record.id] = str(record.seq)
    return dict_name

def reverse_complement(seq):
    """
    #>>> seq = "TCGGinsGCCC"
    #>>> print "Reverse Complement:"
    #>>> print(reverse_complement(seq))
    #GGGCinsCCGA
    """
    alt_map = {'ins': '0'}
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    # for k, v in alt_map.iteritems():
    for k, v in alt_map.items():
        seq = seq.replace(k, v)
    bases = list(seq)
    bases = reversed([complement.get(base, base) for base in bases])
    bases = ''.join(bases)
    # for k, v in alt_map.iteritems():
    for k, v in alt_map.items():
        bases = bases.replace(v, k)
    return bases

def get_promotor(gene_id, gff_dict, fasta_dict):
    gene,rec_id = gff_dict[gene_id]
    
    length = 3100
    if gene.strand == 1:
        base_point = gene.location.start + 1
        get_site = [base_point - length, base_point]
        promotor_seq = fasta_dict[rec_id][get_site[0]:get_site[1]]
    else:
        base_point = gene.location.end
        get_site = [base_point, base_point + length]
        promotor_seq = fasta_dict[rec_id][get_site[0]:get_site[1]]
        promotor_seq = reverse_complement(promotor_seq)
    
    return promotor_seq

if __name__ == '__main__':

    id_file = "D:\pylib\study\geneid.txt"
    out_file = "D:\pylib\study\out_file"
    fasta_file = "D:\pylib\study\A188.fasta"
    gff_file = "D:\pylib\study\A188.gff"
    
    gene_id_list = read_list_file(id_file)
    gff_dict = gff_reader(gff_file)
    fasta_dict = read_fasta(fasta_file)
    
    pro_output_dict = {}
    for gene_id in gene_id_list:
        pro_output_dict[gene_id] = get_promotor(gene_id,gff_dict,fasta_dict)
    
    with open(out_file,'w') as f:
        for gene_id in pro_output_dict.keys():
            f.write(">%s\n%s\n" % (gene_id, pro_output_dict[gene_id]))
    
    
    






