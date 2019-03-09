#!/usr/bin/env python
# -*- coding: utf-8 -*-
# @author xlj

import os
from multiprocessing import Pool

"""
Function: This script was used to calculate the m6A methylation of gene.
Input/Output data: as described in get_gene_methylation.txt.
"""

######################################################################
#directory of peak methylation level in each tissue
peak_m6A_dir = ''
#peak and related gene in each tissue
peak_gene_dir = ''
#Length of longest insform each gene in gencode v27 annotation file
len_file = ''
Result_file = ''
######################################################################


#get tissue name
def get_peafix():
    Prefix = []
    for i in os.listdir(peak_gene_dir):
        if i.endswith('.bed'):
            Prefix.append(i.split('.')[0])
    return Prefix

def run_step(file):
        a = open("%s" %(os.path.join(peak_m6A_dir,file)),'r+').readlines()
        b = open("%s.bed" %(os.path.join(peak_gene_dir,file)),'r+').readlines()
        c = open("%s" %(len_file),'r+').readlines()
        d = open("%s" %(os.path.join(Result_file,file)),'w')
        # peak methylation dic,gene length dic,gene methylation dic
        dica,dicb,dicc= {},{},{}
        for i in a:
            i_ = i.split()
            dica[i_[0]] = i_[1]  
        for j in c:
            j_ = j.split()
            dicb[j_[0]] = int(j_[1])
        for k in b:
            k_ = k.split()
            if k_[3] not in dicc:
                dicc[k_[3]] = float(dica[k_[4]]) * (int(k_[2]) - int(k_[1])) / dicb[k_[3]]
            else:
                dicc[k_[3]] += float(dica[k_[4]]) * (int(k_[2]) - int(k_[1])) / dicb[k_[3]]
        for m in dicc:
            s = str(m + '\t' + str(dicc[m])).replace('\'', '')
            s += '\n'
            d.writelines(s)
        d.close()


if __name__ == "__main__":
    pool = Pool(8)
    Prefix = get_peafix()
    for file in Prefix:
        pool.apply_async(run_step, (file,))
    pool.close()
    pool.join()

