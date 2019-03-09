#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

# size match: total length of chromosome
import os
import glob
import datetime
import pandas as pd
import numpy as np
from multiprocessing import Pool
from sklearn.utils import shuffle

"""
Function: This script was used to generate 100 sets of randomly sampled size-matched input control region according to
            m6a bed of the same tissue.
Input/Output data: as described in choose_control_samples.txt.
"""

#################################################################################
CYCLE_NUMBER = 100
m6a_bed_dir = ""
control_bed_dir = ""
result_dir = ""
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
####################################################################################


def get_input_list():
    os.chdir(control_bed_dir)
    bed_list = glob.glob("*_discrete.bed")
    bed_list = [os.path.abspath(bed) for bed in bed_list]
    return bed_list


def get_sample_name_from_input(infile):
    sample_name = os.path.basename(infile).split("_discrete.bed")[0].lower()
    return sample_name


def get_m6a_list():
    os.chdir(m6a_bed_dir)
    bed_list = glob.glob("*.bed")
    bed_list = [os.path.abspath(bed) for bed in bed_list]
    return bed_list


def get_sample_name_from_m6a(infile):
    sample_name = os.path.basename(infile).split(".bed")[0].lower()
    return sample_name


def pick_control_peak(ip_bed, input_bed):
    chromosome_len_dict = count_ip_length_byChrom(ip_bed)
    df = pd.read_table(input_bed, sep="\t", header=None, comment="#", names=["chr", "start", "end"])
    df.loc[:, "length"] = df["end"] - df["start"]
    sample_name = get_sample_name_from_input(input_bed)
    pool = Pool(processes=40)
    for i in range(CYCLE_NUMBER):
        df = shuffle(df)
        # df_chr = df_chr.sample(frac=1).reset_index(drop=True)
        pool.apply_async(pick_control_length_byChrom, (i+1, sample_name, df, chromosome_len_dict))
        # pick_control_length_byChrom(i + 1, sample_name, df, chromosome_len_dict)
    pool.close()
    pool.join()


def count_ip_length_byChrom(ip_bed):
    chromosome_len_dict, df_list = {}, []
    df = pd.read_table(ip_bed, sep="\t", header=None, comment="#", names=["chr", "start", "end"])
    df.loc[:, "length"] = df["end"] - df["start"]
    chr_list = list(set(df.loc[:, "chr"].tolist()))
    for chr_name in chr_list:
        df_i = df[df["chr"] == chr_name]
        chr_sum_len = df_i["length"].sum()
        chromosome_len_dict[chr_name] = chr_sum_len
    return chromosome_len_dict


def pick_control_length_byChrom(index, sample_name, df, chr_len_dict):
    total_list = [[np.nan, 0, 0]]
    chr_list = list(chr_len_dict.keys())
    for chr_name in chr_list:
        df_chr = df[df["chr"] == chr_name]
        chr_length = chr_len_dict[chr_name]
        peak_list = pick_equal_length_byChrom(chr_length, df_chr)
        total_list += peak_list
    df_result = pd.DataFrame(total_list, columns=["chr", "start", "end"]).dropna(how="any")
    result_file = os.path.join(result_dir, "control-%s_%d.bed" % (sample_name, index))
    df_sorted = df_result.sort_values(["chr", "start"])
    df_sorted.to_csv(result_file, sep="\t", header=None, index=False)
    

def pick_equal_length_byChrom(chr_length, df_chr):
    df_chr.columns = ["chr", "start", "end", "length"]
    count, peak_list = 0, []
    if df_chr["length"].sum() < chr_length:
        print("Input chr length less than IP!")
    for name, values in df_chr.iterrows():
        peak_list.append([values["chr"], values["start"], values["end"]])
        count += values["length"]
        if count >= chr_length:
            break
    chr_name, start, end = peak_list[-1]
    chr_name, start, new_end = chr_name, int(start), (int(end) - (count - chr_length))
    peak_list[-1] = [chr_name, start, new_end]
    return peak_list


def get_input_according_m6a(m6a_bed):
    i_m6a_name = get_sample_name_from_m6a(m6a_bed)
    i_control_file = ""
    input_list = get_input_list()
    for i_input in input_list:
        input_name = get_sample_name_from_input(i_input)
        if input_name == i_m6a_name:
            print(input_name, i_m6a_name)
            i_control_file = i_input
            break
    return i_control_file


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    all_m6a_list = get_m6a_list()
    all_input_list = get_input_list()
    print("m6a list length %d" % len(all_m6a_list))
    for m6a_file in all_m6a_list:
        control_file = get_input_according_m6a(m6a_file)
        pick_control_peak(m6a_file, control_file)
    print(start_time)
    end_time = datetime.datetime.now()
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))

