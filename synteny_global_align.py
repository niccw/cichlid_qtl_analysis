#!/usr/bin/env python3
import sys
import os
import argparse
from collections import OrderedDict
from Bio import Seq
from Bio import SeqIO
from Bio import Align

""" read the genome file (index it, don't read in the memory)
read the bed files (open the target and query bed simutanuosly and put the pair record into LIST)
for each record, extreact the DNA sequences and perform pairwise alignment
write the alignment to the output
 > rename the ID! include the collinearity ID, qtl ID and coordinate info
"""

def read_beds(bed_a_path:str, bed_b_path:str)->list:
    l = []
    with open(bed_a_path,"r") as f, open (bed_b_path,"r") as g:
        for x,y in zip (f,g):
            cola = x.strip().split("\t")
            colb = y.strip().split("\t")
            l.append([[cola[0],int(cola[1]),int(cola[2])],[colb[0],int(colb[1]),int(colb[2])],[cola[3],cola[4],cola[5],cola[6]]])
    print(f"Finish reading {len(l)} pairs of collineary region.",file=sys.stderr)
    return l

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

def find_seqB_point_coord(seqA,seqB,seqA_coord:list,seqB_coord:list,alignment_sign:str,seqA_point_coord:int,aligner)->int:
    # alignment
    if alignment_sign == "-":
        seqB = seqB.reverse_complement()
    alignments = aligner.align(seqA.seq, seqB.seq)
    a = alignments[0]
    [seqA_aligned,_,seqB_aligned,_] = str(a).split("\n")
    base_cnt = seqA_point_coord - seqA_coord[1]
    pos_cnt = 0
    # 1. find the aligned position in seqA (base -> pos in seqA)
    for i in range(len(seqA_aligned)):
        if base_cnt == 0:
            break
        pos_cnt += 1
        if seqA[i] != "-":
            base_cnt -=1
    # 2. find the base in seqB (pos -> base in seqB)
    base_cnt = 0
    for i in range(len(seqB_aligned)):
        if pos_cnt == 0:
            break
        pos_cnt -= 1
        if seqB[i] != "-":
            base_cnt += 1
    # now the base cnt is the coord in seqB :) (0-base)
    base_cnt += 1  # convert to 1-base
    if alignment_sign == "+":
        return seqB_coord[1] + base_cnt
    else:
        return seqB_coord[2] - base_cnt

def global_alignment(ga_path:str, gb_path:str, collinearity_list:list, qtl_d:dict):
    # out put file
    qtl_map = open("qtl_collinearity_map.tsv","w+")
    # index genomes
    # build index name
    print("Indexing genome...", file=sys.stderr)
    if not os.path.exists(f"{ga_path}.idx"): 
        print(f"Writing index to {ga_path}.idx...", file=sys.stderr)
        ga_dict = SeqIO.index_db(f"{ga_path}.idx",ga_path, format="fasta")
    else:
        print(f"Find {ga_path}.idx", file=sys.stderr)
        ga_dict = SeqIO.index_db(f"{ga_path}.idx")
    if not os.path.exists(f"{gb_path}.idx"):
        print(f"Writing index to {gb_path}.idx...", file=sys.stderr)
        gb_dict = SeqIO.index_db(f"{gb_path}.idx",gb_path, format="fasta")
    else:
        print(f"Find {gb_path}.idx", file=sys.stderr)
        gb_dict = SeqIO.index_db(f"{gb_path}.idx")
    print("Finish indexing genome.", file=sys.stderr)
    #  s = genome_dict[i.seqid][i.start:i.end]
    # initiate aligner
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"
    # tune the aligner parameter, high hanging open gap penalty, high mismatch score, small gap penalty
    aligner.mismatch = -1.5
    aligner.open_gap_score = -2.5
    aligner.extend_gap_score = -5

    for l in collinearity_list:
        target_coord = l[0]
        query_coord = l[1]
        pair_info = l[2]
        # don't do alignment if it is "within"! save time!!!
        qtl_l = qtl_d[pair_info[1]]
        # print the qtl corresponding coord in species B
        if target_coord[1] >= qtl_l[1] and target_coord[2] <= qtl_l[2]:
            # don't have to do the alignment :)
            qtl_map.write(f"{pair_info[1]}\t{pair_info[0]}\twithin\t{target_coord[0]}\t{target_coord[1]}\t{target_coord[2]}\t{query_coord[0]}\t{query_coord[1]}\t{query_coord[2]}\n")
        else:
            """
            Have to do alignment, this is a bit complicated. There are three scenario we have to do the alignment: L hanging, R hanging, over. It is hard to align the whole sequence since the synteny is too long. The two speices are phylogentically closed so the synteny size is quite conserve among two species. Therefore, here we extract 10000bp around the point coord of interest for the alignment.
            """
            seqA_start = seqB_start = seqA_end = seqB_end = overlap_type = None

            # classify the type of overlapping (L hanging, R hanging or over)
            # subset sequences and compute the coord(s) on seqB
            if target_coord[1] >= qtl_l[1]: # R hanging
                overlap_type = "r_hanging"
                seqA_start = target_coord[1]
                seqB_start = query_coord[1]
                seqA_end = qtl_l[2]
                seqA_end_subseq = ga_dict[target_coord[0]][seqA_end-5000:seqA_end+5000]
                seqA_end_subseq_coord = [target_coord[0],seqA_end-5000,seqA_end+5000]
                seqB_end = query_coord[1] + abs(seqA_end-target_coord[1]) # tentative seqB_end
                seqB_end_subseq = gb_dict[query_coord[0]][seqB_end-5000:seqB_end+5000]
                seqB_end_subseq_coord = [query_coord[0],seqB_end-5000, seqB_end+5000]

                seqB_end = find_seqB_point_coord(seqA=seqA_end_subseq,seqB=seqB_end_subseq,seqA_coord = seqA_end_subseq_coord,seqB_coord=seqB_end_subseq_coord,alignment_sign=pair_info[3],seqA_point_coord=seqA_end, aligner=aligner)

            if target_coord[2] <= qtl_l[2]: # L hanging
                overlap_type = "l_hangning"
                seqA_end = target_coord[2]
                seqB_end = query_coord[2]
                seqA_start = qtl_l[1]
                seqA_start_subseq = ga_dict[target_coord[0]][seqA_start-5000:seqA_start+5000]
                seqA_start_subseq_coord = [target_coord[0],seqA_start-5000, seqA_start+5000]
                seqB_start = query_coord[1] + abs(seqA_start-target_coord[1])
                seqB_start_subseq = gb_dict[query_coord[0]][seqB_start-5000:seqB_start+5000]
                seqB_start_subseq_coord = [query_coord[0],seqB_start-5000, seqB_start+5000]  
                seqB_start = find_seqB_point_coord(seqA=seqA_start_subseq,seqB=seqB_start_subseq, seqA_coord=seqA_start_subseq_coord, seqB_coord=seqB_start_subseq_coord,alignment_sign=pair_info[3],seqA_point_coord=seqA_start,aligner=aligner)


            if target_coord[1] <= qtl_l[1] and target_coord[2] >= qtl_l[2]: # over (synteny region > qtl region)
                overlap_type = "over"
                seqA_start = qtl_l[1]
                seqA_end = qtl_l[2]
                seqA_start_subseq = ga_dict[target_coord[0]][seqA_start-5000:seqA_start+5000]
                seqA_start_subseq_coord = [target_coord[0],seqA_start-5000, seqA_start+5000]
                seqA_end_subseq = ga_dict[target_coord[0]][seqA_end-5000:seqA_end+5000]
                seqA_end_subseq_coord = [target_coord[0],seqA_end-5000,seqA_end+5000]
                seqB_start = query_coord[1] + abs(seqA_start-target_coord[1])
                seqB_end = query_coord[1] + abs(seqA_end-target_coord[1]) # tentative seqB_end
                seqB_start_subseq = gb_dict[query_coord[0]][seqB_start-5000:seqB_start+5000]
                seqB_start_subseq_coord = [query_coord[0],seqB_start-5000, seqB_start+5000]
                seqB_end_subseq = gb_dict[query_coord[0]][seqB_end-5000:seqB_end+5000]
                seqB_end_subseq_coord = [query_coord[0],seqB_end-5000, seqB_end+5000]
                
                seqB_start = find_seqB_point_coord(seqA=seqA_start_subseq,seqB=seqB_start_subseq, seqA_coord=seqA_start_subseq_coord, seqB_coord=seqB_start_subseq_coord,alignment_sign=pair_info[3],seqA_point_coord=seqA_start,aligner=aligner)
                seqB_end = find_seqB_point_coord(seqA=seqA_end_subseq,seqB=seqB_end_subseq,seqA_coord = seqA_end_subseq_coord,seqB_coord=seqB_end_subseq_coord,alignment_sign=pair_info[3],seqA_point_coord=seqA_end, aligner=aligner)
        
            # reverse seqB coord if "-"
            qtl_map.write(f"{pair_info[1]}\t{pair_info[0]}\t{overlap_type}\t{target_coord[1]}\t{seqA_start}\t{seqA_end}\t{query_coord[0]}\t{seqB_start}\t{seqB_end}\n")
    qtl_map.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-ba', type=str, help='Species A BED file', required=True)
    parser.add_argument('-bb', type=str, help='Species B BED file', required=True)
    parser.add_argument('-ga', type=str, help='Species A genome fasta', required=True)
    parser.add_argument('-gb', type=str, help='Species B genome fasta', required=True)
    parser.add_argument('-q', type=str, help='qtl_regions_reformat.tsv', required=True)
    args = parser.parse_args()

    pair_l = read_beds(bed_a_path=args.ba, bed_b_path=args.bb)
    qtl_d = read_qtl_regions(qtl_regions_path=args.q)
    global_alignment(ga_path=args.ga, gb_path=args.gb,collinearity_list=pair_l,qtl_d=qtl_d)

