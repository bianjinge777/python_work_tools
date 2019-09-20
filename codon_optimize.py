def read_cds_split(txtfile):
    with open(txtfile, 'r') as cds_data:
        tmp_string = ''
        for single_cds in cds_data:
            single_cds = single_cds.strip()
            tmp_string += single_cds
        split_cds = [tmp_string[i:i + 3] for i in range(0, len(tmp_string), 3)]
        return split_cds


def replace_cds(cds_list):
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
    table2 = [["TTC", "F"], ["TTG", "L"], ["TCT", "S"], ["TAC", "Y"],
              ["TGC", "C"], ["TGG", "W"], ["CTC", "L"], ["CCT", "P"],
              ["CAT", "H"], ["CAG", "Q"], ["CGC", "R"], ["ATC", "I"],
              ["ATG", "M"], ["ACC", "T"], ["AAC", "N"], ["AAG", "K"],
              ["AGC", "S"], ["AGG", "R"], ["GTC", "V"], ["GCT", "A"],
              ["GAT", "D"], ["GAG", "E"], ["GGA", "G"]]

    cds_pron = []  # 将cds_list的密码子对应相应的氨基酸然后存放cds_pron中
    for tmp_cds in cds_list:
        for standard_table in table1:
            if tmp_cds == standard_table[0]:
                tmp_cds = standard_table
                cds_pron.append(tmp_cds)

    result_list = []  # 存放密码子和氨基酸，密码子替换成使用率最高的
    for list_cds in cds_pron:
        for table_key in table2:
            if table_key[1] == list_cds[1]:
                list_cds[0] = table_key[0]
                result_list.append(list_cds)

    result_string = ''
    for get_cds in result_list:
        result_string += get_cds[0]
    return result_string


if __name__ == '__main__':
    txtfile = "D:\pylib\study\HPH Gene.txt"
    output = "result.txt"

    cds_list = read_cds_split(txtfile)

with open(output, 'w') as outputfile:
    outputfile.write(str(replace_cds(cds_list)))
