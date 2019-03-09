#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import re
import math
import glob
import subprocess
from scipy import stats

"""
Function: This script was used to calculate the observed-to-expected CpG ratio in the promoter.
Input/Output data: as described in calculate_CpG_OE_ratio.txt.
"""

#################################################################################################################
type_one, type_two = "m6a", "free"
m6a_promoter_dir = ""
free_promoter_dir = ""
reference_genome = ""
result_dir = ""
if not os.path.exists(result_dir):
	os.makedirs(result_dir)
result_data_file = "%s/%s_vs_%s.txt" % (result_dir, type_one, type_two)
result_stat_file = "%s/stat_%s_vs_%s.txt" % (result_dir, type_one, type_two)
for i_file in [result_data_file, result_stat_file]:
	if os.path.exists(i_file):
		os.remove(i_file)
##################################################################################################################


def compare_CpG_ratio():
	m6a_dict, free_dict = generate_dict(m6a_promoter_dir), generate_dict(free_promoter_dir)
	for tissue, m6a_promoter_bed in m6a_dict.items():
		free_promoter_bed = get_same_tissue_bed(tissue, free_dict)
		print(m6a_promoter_bed)
		print(free_promoter_bed)
		m6a_seq, free_seq = get_fasta(m6a_promoter_bed), get_fasta(free_promoter_bed)
		m6a_score_list, free_score_list = calculate_CpG_density(m6a_seq), calculate_CpG_density(free_seq)
		write_score_into_file(tissue, m6a_score_list, free_score_list)
		statistic, pvalue = statistic_significance(m6a_score_list, free_score_list)
		write_stat_result(tissue, statistic, pvalue)


def generate_dict(data_dir):
	result_dict = {}
	bed_list = glob.glob("%s/*.bed" % data_dir)
	for bed in bed_list:
		result_dict[os.path.basename(bed).split(".bed")[0].lower()] = bed
	return result_dict


def get_same_tissue_bed(tissue, free_dict):
	tissue_name, promoter_bed = "", ""
	try:
		tissue_name = tissue.lower()
		promoter_bed = free_dict[tissue_name]
	except KeyError:
		print("%s not in free_dict" % tissue_name)
	return promoter_bed


def get_fasta(bed):
	command = "bedtools getfasta -fi %s -bed %s" % (reference_genome, bed)
	sub_p = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
	fa_sequence = str(sub_p.communicate()[0].decode("utf-8"))
	return fa_sequence


def calculate_CpG_density(fa_sequence):
	score_list = []
	for line in fa_sequence.strip().split("\n"):
		if not line.startswith(">"):
			seq = line.strip().lower()
			cg_num, c_num, g_num = len(re.findall("cg", seq)), len(re.findall("c", seq)), len(re.findall("g", seq))
			score = cg_num / (math.pow((c_num + g_num)/2.0, 2) / (len(seq)*1.0))
			print(score)
			score_list.append(score)
	return score_list


def statistic_significance(m6a_scores, free_scores):
	statistic, pvalue = stats.mannwhitneyu(m6a_scores, free_scores)
	return statistic, pvalue


def write_score_into_file(tissue_name, pos_list, neg_list):
	with open(result_data_file, 'a') as fw:
		for score in pos_list:
			fw.write("%s\t%s\t%f\n" % (tissue_name, type_one, score))
		for score in neg_list:
			fw.write("%s\t%s\t%f\n" % (tissue_name, type_two, score))


def write_stat_result(tissue_name, stat, pvalue):
	with open(result_stat_file, 'a') as fw:
		fw.write("%s\t%f\t%f\n" % (tissue_name, stat, pvalue))


if __name__ == '__main__':
	compare_CpG_ratio()