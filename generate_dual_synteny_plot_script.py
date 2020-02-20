#!/usr/bin/env python3

import sys
import os
from collections import OrderedDict

# hard-coded
target_species_chr_prefix = "ScC7KKH"
query_species_chr_prefix = "Seq_"

# plot setting
plot_w = 600
plot_h = 800

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


def read_collinearity_range(collinearity_range_path:str)->dict:
    d = OrderedDict()
    with open(collinearity_range_path,"r") as f:
        for line in f:
            col = line.strip().split("\t")
            if col[4] not in d:
                d[col[4]] = [[],[]] # target-target list, target-query list
            if col[7].startswith(target_species_chr_prefix):
                d[col[4]][0].append(col[7])
            else:
                d[col[4]][1].append(col[7])
    return d


def generate_dual_synteny_plot_script(qtl_d:dict,col_d:dict)->None:
    if not os.path.exists("dual_synteny_ctl"):
        os.makedirs("dual_synteny_ctl")
    # one script for each qtl
    for k,v in qtl_d.items():
        if v[0] not in col_d:
            print(f"{k}: {v[0]} not in collineraity file.",file=sys.stderr)
            continue
        else:
            # create ctl file
            file_name = f"{k}.ctl"
            print(file_name)
            with open(f"./dual_synteny_ctl/{file_name}","w+") as f:
                target_chr = [v[0]]
                target_chr.extend(col_d[v[0]][0])
                query_chr = col_d[v[0]][1]
                # remove duplicate
                target_chr = list(set(target_chr))
                query_chr = list(set(query_chr))
                target_chr_str = ",".join(target_chr)
                query_chr_str = ",".join(query_chr)

                f.write(f"{plot_w}\n")
                f.write(f"{plot_h}\n")
                f.write(f"{target_chr_str}\n")
                f.write(f"{query_chr_str}\n")


if __name__ == "__main__":
    qtl_regions_path = sys.argv[1]
    collinearity_range_path = sys.argv[2]

    qtl_d = read_qtl_regions(qtl_regions_path=qtl_regions_path)
    collinearity_chr_d = read_collinearity_range(collinearity_range_path=collinearity_range_path)
    generate_dual_synteny_plot_script(qtl_d=qtl_d, col_d=collinearity_chr_d)
