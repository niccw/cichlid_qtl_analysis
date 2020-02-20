#!/usr/bin/env python3

import sys
import re
import pdb

def parse_gene_gff(gene_gff_path:str)->dict:
    gene_d = {}
    with open(gene_gff_path, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.strip().split("\t")
            gene_d[col[1]] = [col[0], int(col[2]), int(col[3])]
    return gene_d

def minmaxloc(num_list:list)->list:
    return [num_list.index(min(num_list)), num_list.index(max(num_list))]

def parse_mscanx(collinearity_path:str, gene_d:dict):
    """
    Annotate the collinearity file with coordinate
    Extract the 'range' of each alignment block 
    """
    alignment_idx = None
    alignment_size = 0
    alignment_sign = None
    # tmp lists to hold the value for comparison
    alignment_gene_target_ls_tmp = []
    alignment_gene_query_ls_tmp = []
    alignment_start_target_ls_tmp = []
    alignment_start_query_ls_tmp = []
    alignment_end_target_ls_tmp = []
    alignment_end_query_ls_tmp = []
    
    # file handler
    annotated_collinearity = open("annotated_collinearity.tsv","w")
    block_range = open("collinearity_block_range.tsv", "w")
    with open(collinearity_path, "r") as f:
        for line in f: # comment line, check whether it is header or alignment line
            if line.startswith("# "):
                annotated_collinearity.write(f"{line.strip()}\n")
                continue
            if line.startswith("## "):
                alignment_idx = re.search(r"(?<=Alignment\s)\d+(?=:)",line)[0]
                alignment_size = re.search(r"(?<=N=)\d+",line)[0]
                alignment_sign = re.search("minus",line)
                if alignment_sign is None:
                    alignment_sign = "+"
                else:
                    alignment_sign = "-"
                annotated_collinearity.write(f"{line.strip()}\n") 
                continue
            if line.startswith("###"):
                annotated_collinearity.write(f"{line.strip()}\n")
                continue
                        
            col = line.strip().split("\t")
            # annotate the collinearity file
            try:
                col.extend(gene_d[col[1]])
                col.extend(gene_d[col[2]])
            except Exception as e:
                print(f"{col}", file = sys.stdout)
                print(f"{e}", file = sys.stdout)
            col = [str(x) for x in col]
            col_str = "\t".join(col)
            annotated_collinearity.write(f"{col_str}\n")

            # extract the range of alignment
            alignment_id_g = re.split(r"-\s+",col[0])
            if len(alignment_id_g) == 1:
                alignment_id_g = re.split(r"-",col[0])
            alignment_id_g[1] = alignment_id_g[1].replace(":","")
            ## check whether the row is the start or end of the alignment block
            if alignment_idx != None and int(alignment_id_g[1]) <= (int(alignment_size) - 1):
                # record the gene id
                alignment_gene_target_ls_tmp.append(col[1])
                alignment_gene_query_ls_tmp.append(col[2])
                # record the start pos
                alignment_start_target_ls_tmp.append(gene_d[col[1]][1])
                alignment_start_query_ls_tmp.append(gene_d[col[2]][1])
                # record the end pos
                alignment_end_target_ls_tmp.append(gene_d[col[1]][2])
                alignment_end_query_ls_tmp.append(gene_d[col[2]][2])
                if int(alignment_id_g[1]) < (int(alignment_size) - 1):
                    continue
            ## reach the last record of this alignment block, check the min and max
            if int(alignment_id_g[1]) == int(alignment_size)-1 :
                # find the index of 1) min(target_start) 2) max(target_end) 3) min(query_start) 4) max(query_end)
                min_target_start_i,_ = minmaxloc(alignment_start_target_ls_tmp)
                _,max_target_end_i = minmaxloc(alignment_end_target_ls_tmp)
                min_quert_start_i,_= minmaxloc(alignment_start_query_ls_tmp)
                _,max_query_end_i = minmaxloc(alignment_end_query_ls_tmp)
                min_target_gene = alignment_gene_target_ls_tmp[min_target_start_i]
                max_target_gene = alignment_gene_target_ls_tmp[max_target_end_i]
                min_query_gene = alignment_gene_query_ls_tmp[min_quert_start_i]
                max_query_gene = alignment_gene_query_ls_tmp[max_query_end_i]
                
                block_range.write(f"{alignment_idx}\tstart\t{min_target_gene}\t{min_query_gene}\t{gene_d[min_target_gene][0]}\t{gene_d[min_target_gene][1]}\t{gene_d[min_target_gene][2]}\t{gene_d[min_query_gene][0]}\t{gene_d[min_query_gene][1]}\t{gene_d[min_query_gene][2]}\t{alignment_sign}\n")
                block_range.write(f"{alignment_idx}\tend\t{max_target_gene}\t{max_query_gene}\t{gene_d[max_target_gene][0]}\t{gene_d[max_target_gene][1]}\t{gene_d[max_target_gene][2]}\t{gene_d[max_query_gene][0]}\t{gene_d[max_query_gene][1]}\t{gene_d[max_query_gene][2]}\t{alignment_sign}\n")
                # reset variables
                alignment_idx = None
                alignment_size = 0
                alignment_gene_query_ls_tmp = []
                alignment_gene_target_ls_tmp = []
                alignment_start_target_ls_tmp = []
                alignment_start_query_ls_tmp =[]
                alignment_end_target_ls_tmp = []
                alignment_end_query_ls_tmp = []
    annotated_collinearity.close()
    block_range.close()

if __name__ == "__main__":
    collinearity_path = sys.argv[1]
    gene_gff_path = sys.argv[2]

    d = parse_gene_gff(gene_gff_path)
    parse_mscanx(collinearity_path=collinearity_path,gene_d=d)
