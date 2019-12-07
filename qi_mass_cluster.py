import os
import pandas as pd
from decimal import Decimal
import csv

def get_filelist(path):
    filelist = os.listdir(path)  # 列出文件夹下所有的目录与文件
    all_RT = []
    all_dist_list = []
    for filename in filelist:
        open_file = path + '/'+ filename
        with open(open_file, 'r', encoding='utf-8') as fr:
            imp_dist = {}
            splice_str = ''
            for line in fr.readlines()[2:]:
                line = line.strip()
                splice_str += line + ','
                out_list = splice_str.split('RT ; Compound Name ; SI ; RSI ; Probability ; Molecular Weight ; Molecular Formula ; Cas # ; Library',1)[0].split(',')
            out_list.pop()
            for data in out_list:
                data = data.split(';')
                imp_dist[data[0]] = data[1]
            for each_key in imp_dist.keys():
                all_RT.append(each_key)
            all_dist_list.append(imp_dist)
    # print(all_RT)
    # print(all_dist_list)
    return all_RT, all_dist_list

def get_RT_intetval(all_RT):
    interval_list = []
    for each_Rt in all_RT:
        min_each_key = Decimal(each_Rt) - Decimal('0.05')
        max_each_key = Decimal(each_Rt) + Decimal('0.05')
        interval_list.append((float(min_each_key), float(max_each_key)))
    return interval_list

def merge_intervals(input_list, int=False):
    intervals = []
    for i in input_list:
        intervals.append(tuple(i))
    sorted_by_lower_bound = sorted(intervals, key=lambda tup: min(tup))
    merged = []

    for higher in sorted_by_lower_bound:
        if not merged:
            merged.append(higher)
        else:
            lower = merged[-1]
            # test for intersection between lower and higher:
            # we know via sorting that lower[0] <= higher[0]
            if int is False:
                if min(higher) <= max(lower):
                    upper_bound = max(lower + higher)
                    merged[-1] = (min(lower), upper_bound)  # replace by merged interval
                else:
                    merged.append(higher)
            elif int is True:
                if min(higher) <= max(lower):
                    upper_bound = max(lower + higher)
                    merged[-1] = (min(lower), upper_bound)  # replace by merged interval
                elif max(lower) + 1 == min(higher):
                    upper_bound = max(lower + higher)
                    merged[-1] = (min(lower), upper_bound)  # replace by merged interval
                else:
                    merged.append(higher)
    # print(merged)
    return merged

def RT_list_by_intervals(all_RT,merged,all_dist_list):
    merge_dist = {}
    for (min, max) in merged:
        for each_RT in all_RT:
            if min <= float(each_RT) <= max:
                merge_dist.setdefault((min, max), []).append(each_RT)
    print(merge_dist)

    result_dist = {}
    for each_dist in all_dist_list:
        for each_dist_RT, each_dist_Area in each_dist.items():
            for merge_dist_value in merge_dist.values():
                for each_merge_dist_value in merge_dist_value:
                    if each_dist_RT == each_merge_dist_value:
                        result_dist.setdefault(tuple(merge_dist_value), []).append(each_dist_Area)

    print(result_dist)
    return result_dist
def dict_to_csv(input_dict):
    result_list = []  
    for k, v in input_dict.items():
        single_dict = {'RT': k, 'Area': v}
        result_list.append(single_dict)
    printfile = pd.DataFrame(result_list)
    printfile.to_csv('out.csv')
    print(result_list)

if __name__ == '__main__':
    path = "D:\pylib\study\qi_mass"

    all_RT = get_filelist(path)[0]
    all_dist_list = get_filelist(path)[1]
    input_list = get_RT_intetval(all_RT)
    merged = merge_intervals(input_list, int=False)
    input_dict = RT_list_by_intervals(all_RT, merged,all_dist_list)
    dict_to_csv(input_dict)
