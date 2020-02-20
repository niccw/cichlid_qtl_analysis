#!/usr/bin/env python3

import sys
from collections import OrderedDict

def read_qtl_regions(qtl_regions_path:str)->dict:
    d = OrderedDict()
    last_id_root =""
    cnt = 0
    with open(qtl_regions_path,"r") as f:
        next(f)
        for line in f:
            col = line.strip().split("\t")
            if col[4] == "-":
                continue
            qtl_id_root =  f"{col[0]}-{col[1]}-{col[2]}"
            if qtl_id_root == last_id_root:
                cnt += 1
                qtl_id = f"{col[0]}-{col[1]}-{col[2]}-{cnt}"
                last_id_root = qtl_id_root
            else:
                cnt = 0
                qtl_id = f"{col[0]}-{col[1]}-{col[2]}-{cnt}"
                last_id_root = qtl_id_root
            d[qtl_id] = [col[3],int(col[4]),int(col[5])]
    return d

def filter_collinearity_block(qtl_d:dict,collinearity_block_path:str):
    col_d = {}
    target_bed = open("target_block.bed","w+")
    query_bed = open("query_block.bed", "w+")
    with open(collinearity_block_path,"r") as f:
        for line in f:
            col = line.strip().split("\t")
            if col[0] not in col_d: # initiate k,v pair
                col_d[col[0]] = [[col[4]],[col[7]]] # chr,start,end in species A and B 
            if col[1] == "start":
                col_d[col[0]][0].append(int(col[5]))
                col_d[col[0]][1].append(int(col[8]))
            if col[1] == "end":
                col_d[col[0]][0].append(int(col[6]))
                col_d[col[0]][1].append(int(col[9]))
                col_d[col[0]][0].append(col[10])
                col_d[col[0]][1].append(col[10])
            # also mark down the alignment sign

    # decide wheather to keep the collineary block or not
    for k,v in col_d.items():
        c_chr = v[0][0]
        c_start = v[0][1]
        c_end = v[0][2]
        for x,y in qtl_d.items():
            q_chr = y[0]
            # TODO: add flag control to show gene duplication in one species
            if q_chr != c_chr or c_chr.split("_")[0] == v[1][0].split("_")[0]: # the fasta header has to be <prefix>_<#> !!!
                continue
            q_start = y[1]
            q_end = y[2]
            if check_overlap(c_start,c_end,q_start,q_end) == 0:
                continue
            # Overlap!
            alignment_sign = v[0][3]
            overlap_size = check_overlap(c_start,c_end,q_start,q_end)
            target_bed.write(f"{c_chr}\t{c_start}\t{c_end}\t{k}\t{x}\t{overlap_size}\t{alignment_sign}\n")
            query_bed.write(f"{v[1][0]}\t{v[1][1]}\t{v[1][2]}\t{k}\t{x}\t{overlap_size}\t{alignment_sign}\n")

    target_bed.close()
    query_bed.close()

def check_overlap(min1,max1,min2,max2):
    return max(0, min(max1, max2) - max(min1, min2))

if __name__ == "__main__":
    qtl_d = read_qtl_regions(sys.argv[1])
    filter_collinearity_block(qtl_d=qtl_d, collinearity_block_path=sys.argv[2])

