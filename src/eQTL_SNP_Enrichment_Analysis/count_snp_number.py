#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import glob
import datetime
from multiprocessing import Pool

"""
Runtime: ~8min
Function: This script was used to count how many eQTL SNPs existing in m6A peaks and control peaks(100 random samples).
Input/Output data: as described in count_snp_number.txt.
"""

######################################################################################################################
eQTL_snp_dir = ""
m6a_bed_dir = ""
con_bed_dir = ""
result_dir = ""
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
#######################################################################################################


def get_input_list():
    os.chdir(con_bed_dir)
    bed_list = glob.glob("control-*_*.bed")
    bed_list = [os.path.abspath(bed) for bed in bed_list]
    return bed_list


def get_tissue_name_from_input(input_file):
    tissue_name = os.path.basename(input_file).split("-")[1].split(".bed")[0].lower()
    return tissue_name


def get_m6a_list():
    os.chdir(m6a_bed_dir)
    bed_list = glob.glob("*.bed")
    bed_list = [os.path.abspath(bed) for bed in bed_list]
    return bed_list


def get_tissue_from_m6a(m6a_file):
    tissue_name = os.path.basename(m6a_file).split(".bed")[0].lower()
    return tissue_name


def get_snp_list():
    eqtl_snp_list = glob.glob("%s/*.bed" % eQTL_snp_dir)
    return eqtl_snp_list


def get_tissue_name_from_eQTL(eqtl_snp_bed):
    tissue_name = os.path.basename(eqtl_snp_bed).split(".bed")[0].split("_")[0].lower()
    return tissue_name


def intersect_bed_snp(peak_type, peak_bed, snp_bed):
    snp_name = os.path.basename(snp_bed).split(".bed")[0]
    intersect_result_file = ""
    if peak_type == "m6a":
        m6a_name = get_tissue_from_m6a(peak_bed)
        intersect_result_file = "%s_%s.bed" % (snp_name, m6a_name)
    elif peak_type == "control":
        con_name = get_tissue_name_from_input(peak_bed)
        intersect_result_file = "%s_%s.bed" % (snp_name, con_name)
    os.system("bedtools intersect -a {snp_bed} -b {peak_bed} -wa > {intersect_result_file}"
              .format(snp_bed=snp_bed, peak_bed=peak_bed, intersect_result_file=intersect_result_file))
    print("%s\t%s done!" % (peak_type, snp_name))


def run_all():
    m6a_list = get_m6a_list()
    input_list = get_input_list()
    snp_list = get_snp_list()
    print(len(m6a_list), len(input_list), len(snp_list))
    #
    os.chdir(m6a_bed_dir)
    condition_result_dir = "%s%scondition" % (result_dir, os.sep)
    if not os.path.exists(condition_result_dir):
        os.makedirs(condition_result_dir)
    os.chdir(condition_result_dir)
    for snp in snp_list:
        t_name = get_tissue_name_from_eQTL(snp)
        for m6a in m6a_list:
            name = get_tissue_from_m6a(m6a)
            if name == t_name:
                intersect_bed_snp("m6a", m6a, snp)
    #
    os.chdir(con_bed_dir)
    control_result_dir = "%s%scontrol" % (result_dir, os.sep)
    if not os.path.exists(control_result_dir):
        os.makedirs(control_result_dir)
    os.chdir(control_result_dir)
    pool = Pool(processes=40)
    for snp in snp_list:
        t_name = get_tissue_name_from_eQTL(snp)
        i_input_list = [x for x in input_list if t_name in x]
        for con in i_input_list:
            pool.apply_async(intersect_bed_snp, ("control", con, snp))
    pool.close()
    pool.join()


if __name__ == "__main__":
    start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    run_all()
    print(start_time)
    end_time = datetime.datetime.now()
    print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
