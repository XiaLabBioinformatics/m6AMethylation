#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

"""
Runtime: ~25min
Function: This script was used to do enrichment analysis for eQTL SNPs,
          which may enriched in m6A peaks comparing with control peaks(100 random samples).
Input/Output data: as described in enrichment_analysis.txt.
"""

import os
import glob
import datetime
from scipy import stats
from multiprocessing import Pool
from multiprocessing import Manager

######################################################################################################################
m6a_dir = ""
m6a_snp_dir = ""
con_dir = ""
con_snp_dir = ""
result_dir = ""
if not os.path.exists(result_dir):
	os.makedirs(result_dir)
outfile = ""
######################################################################################################################


def get_input_list():
	os.chdir(con_dir)
	bed_list = glob.glob("control-*_*.bed")  # "control-kid_78.bed"
	bed_list = [os.path.abspath(bed) for bed in bed_list]
	return bed_list


def get_tissue_name_from_input(input_file):
	tissue_name = "_".join(os.path.basename(input_file).split("-")[1].split(".bed")[0].split("_")[:-1]).lower()
	return tissue_name


def get_cycle_name_from_input(input_file):
	cycle_name = os.path.basename(input_file).split("-")[1].split(".bed")[0].split("_")[-1].lower()
	return cycle_name


def get_snp_input_list():
	os.chdir(con_snp_dir)
	# bed_list = glob.glob("*_*_*.bed")  # "response-to-drug_kid_92.bed"
	bed_list = glob.glob("*.bed")
	bed_list = [os.path.abspath(bed) for bed in bed_list]
	return bed_list


def get_tissue_and_class_from_snp_input(snp_input_file):
	class_name = "_".join(os.path.basename(snp_input_file).split("_")[:-2]).lower()
	tissue_name = os.path.basename(snp_input_file).split("_")[-2].lower()
	cycle_name = os.path.basename(snp_input_file).split("_")[-1].split(".bed")[0]
	return class_name, tissue_name, cycle_name


def get_m6a_list():
	os.chdir(m6a_dir)
	bed_list = glob.glob("*.bed")  # pla.bed
	bed_list = [os.path.abspath(bed) for bed in bed_list]
	return bed_list


def get_tissue_name_from_m6a(m6a_file):
	tissue_name = os.path.basename(m6a_file).split(".bed")[0].lower()
	return tissue_name


def get_snp_m6a_list():
	os.chdir(m6a_snp_dir)
	bed_list = glob.glob("*_*.bed")  # "inflammatory-measurement_testis.bed"
	bed_list = [os.path.abspath(bed) for bed in bed_list]
	return bed_list


def get_tissue_and_class_from_snp_m6a(snp_m6a_file):
	tissue_name = os.path.basename(snp_m6a_file).split("_")[-1].split(".bed")[0].lower()
	class_name = "_".join(os.path.basename(snp_m6a_file).split("_")[:-1]).lower().lower()
	return class_name, tissue_name


def statics_bed_sum(input_bed):
	with open(input_bed, 'r') as f:
		sum_interval = 0
		for line in f.readlines():
			if not line.startswith("#"):
				info = line.strip().split("\t")
				start, end = int(info[1]), int(info[2])
				interval_len = end - start
				sum_interval += interval_len
	return sum_interval


def stats_fisher_exact(m6a_bed, con_bed, m6a_snp, con_snp):
	print(os.path.basename(m6a_bed), os.path.basename(con_bed), os.path.basename(m6a_snp), os.path.basename(con_snp))
	m6a_all, con_all = statics_bed_sum(m6a_bed), statics_bed_sum(con_bed)
	m6a_snp = statistic_line_num(m6a_snp)
	con_snp = statistic_line_num(con_snp)
	m6a_unsnp, con_unsnp = (m6a_all - m6a_snp), (con_all - con_snp)
	print(m6a_snp, con_snp, m6a_unsnp, con_unsnp)
	#        m6a   control
	# snp    a       b
	# unsnp  c       d
	if m6a_snp == con_snp == 0 or m6a_unsnp == con_unsnp == 0 or m6a_snp == m6a_unsnp == 0 or con_snp == con_unsnp == 0:
		oddsratio, pvalue = 0, 999
	elif con_snp == 0 or m6a_unsnp == 0:  # Haldane-Anscombe correction
		oddsratio, pvalue = stats.fisher_exact([[m6a_snp, con_snp], [m6a_unsnp, con_unsnp]])
		m6a_snp, con_snp, m6a_unsnp, con_unsnp = m6a_snp + 0.5, con_snp + 0.5, m6a_unsnp + 0.5, con_unsnp + 0.5
		oddsratio = (m6a_snp * con_unsnp) / (con_snp * m6a_unsnp)
	else:
		oddsratio, pvalue = stats.fisher_exact([[m6a_snp, con_snp], [m6a_unsnp, con_unsnp]])
	# sample_name = get_tissue_name_from_m6a(m6a_bed)
	# print("########################")
	# print("sample_name\t%s" % sample_name)
	# print("oddsratio\t%f" % oddsratio)
	# print("pvalue\t%f" % pvalue)
	# print("#\tm6a\tcontrol")
	# print("#snp\t%d\t%d" % (m6a_snp, con_snp))
	# print("#unsnp\t%d\t%d" % (m6a_unsnp, con_unsnp))
	print(oddsratio, pvalue)
	return oddsratio, pvalue


def statistic_line_num(infile):
	with open(infile, 'r') as f:
		line_num = len(f.readlines())
	return line_num


def all_fisher_exact_test(snp_m6a_list, snp_input_list, m6a_list, input_list, queue, lock):
	pool = Pool()
	for snp_input_file in snp_input_list:
		pool.apply_async(fisher_exact_test, args=(snp_m6a_list, snp_input_file, m6a_list, input_list, queue, lock))
	pool.close()
	pool.join()


def fisher_exact_test(snp_m6a_list, snp_input_file, m6a_list, input_list, queue, lock):
	i_result_line = one_hundred_circle(m6a_list, input_list, snp_m6a_list, snp_input_file)
	lock.acquire()
	queue.put(i_result_line)
	lock.release()


def one_hundred_circle(m6a_list, input_list, snp_m6a_list, snp_input_file):
	m6a_file, control_file, con_snp_file, m6a_snp_file = "", "", snp_input_file, ""
	class_name, tissue_name, cycle_name = get_tissue_and_class_from_snp_input(snp_input_file)
	i_key = "%s\t%s" % (tissue_name, class_name)

	for i_m6a in m6a_list:
		m6a_tissue_name = get_tissue_name_from_m6a(i_m6a)
		if m6a_tissue_name == tissue_name:
			m6a_file = i_m6a
			break

	for i_input in input_list:
		input_tissue_name, input_cycle_name = get_tissue_name_from_input(i_input), get_cycle_name_from_input(i_input)
		if input_tissue_name == tissue_name and input_cycle_name == cycle_name:
			control_file = i_input
			break

	for i_m6a_snp in snp_m6a_list:
		i_class_name, i_tissue_name = get_tissue_and_class_from_snp_m6a(i_m6a_snp)
		if i_class_name == class_name and i_tissue_name == tissue_name:
			m6a_snp_file = i_m6a_snp
			break
	i_oddsratio, i_pvalue = stats_fisher_exact(m6a_file, control_file, m6a_snp_file, con_snp_file)
	result_line = "%s\t%f\t%f" % (i_key, i_oddsratio, i_pvalue)
	if len(result_line.split("\t")) != 4:
		print(m6a_file, control_file, m6a_snp_file, con_snp_file)
	return result_line


def write_to_file(result_dict):
	global outfile
	with open(outfile, 'w') as fw:
		fw.write("tissue\ttrait\tmeanOddsratio\n")
		for i_key, oddsratio_list in result_dict.items():
			if len(oddsratio_list) == 0:
				print(i_key)
			new_oddsratio_list = oddsratio_list
			odd_average = sum(new_oddsratio_list) / len(new_oddsratio_list)
			fw.write("%s\t%f\n" % (i_key, odd_average))


# if pvalue == 999 or pvalue > 0.05, ignore oddsratio;
def process_raw_result(queue):
	queue_list, result_dict = [], {}
	while not queue.empty():
		value = queue.get()
		queue_list.append(value)
	queue_list.sort()
	print("queue_list length: %d" % len(queue_list))
	for line in queue_list:
		tissue, trait, oddsratio, pvalue = line.strip().split("\t")
		if int(float(pvalue)) == 999:
			oddsratio = 0
		elif float(pvalue) > 0.05:
			oddsratio = 0
		i_key = "%s\t%s" % (tissue, trait)
		result_dict[i_key] = result_dict.get(i_key, []) + [float(oddsratio)]
	return result_dict


if __name__ == "__main__":
	start_time = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
	all_input_list = get_input_list()
	print("all_input_list: %d" % len(all_input_list))
	all_m6a_list = get_m6a_list()
	print("all_m6a_list: %d" % len(all_m6a_list))
	all_snp_input_list = get_snp_input_list()
	print("all_snp_input_list: %d" % len(all_snp_input_list))
	all_class_list = list(set([get_tissue_and_class_from_snp_input(bed) for bed in all_snp_input_list]))
	print("all_class_list: %d" % len(all_class_list))
	all_snp_m6a_list = get_snp_m6a_list()
	print("all_snp_m6a_list: %d" % len(all_snp_m6a_list))
	manager = Manager()
	m_queue = manager.Queue()
	m_lock = manager.Lock()
	all_fisher_exact_test(all_snp_m6a_list, all_snp_input_list, all_m6a_list, all_input_list, m_queue, m_lock)
	print("queue size: %d" % m_queue.qsize())
	final_dict = process_raw_result(m_queue)
	write_to_file(final_dict)
	print("Fisher exact test has completed!")
	print(start_time)
	end_time = datetime.datetime.now()
	print(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
