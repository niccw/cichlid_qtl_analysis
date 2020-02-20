#!/usr/bin/env python3

import gffutils
import argparse
from collections import OrderedDict
import os
import sys
import subprocess

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import colorlover as cl

def read_collinearity_map(m_path:str)->list:
    """
    [qtl_id,col_id,overlap_type,target_chr,target_start,target_end,query_chr,query_start,query_end]
    """
    l = []
    with open(m_path, "r") as f:
        for line in f:
            col = line.strip().split("\t")
            l.append(col)
    return l

def load_gff_db(ga_path:str, gb_path:str):
    global ga_db
    global gb_db
    # check if the sqlite db already exists
    if not os.path.exists(f"{ga_path}.db"):
        print(f"Creating sqlite db for {ga_path}",file=sys.stderr)
        ga_db = gffutils.create_db(ga_path,f"{ga_path}.db",keep_order=True)
    else:
        ga_db = gffutils.FeatureDB(f"{ga_path}.db",keep_order=True)
    if not os.path.exists(f"{gb_path}.db"):
        print(f"Creating sqlite db for {gb_path}",file=sys.stderr)
        gb_db = gffutils.create_db(gb_path,f"{gb_path}.db",keep_order=True)
    else:
        gb_db = gffutils.FeatureDB(f"{gb_path}.db",keep_order=True)
    print("Finish importing gff db")
        

def all_te_family(ga_path:str, gb_path:str):
    ga_te_byte =subprocess.check_output(f"cut -f3 {ga_path} |sort|uniq|grep -v '^#'", shell=True)
    ga_te_ls = ga_te_byte.decode().split("\n") 
    gb_te_byte =subprocess.check_output(f"cut -f3 {gb_path} |sort|uniq|grep -v '^#'", shell=True)
    gb_te_ls = gb_te_byte.decode().split("\n") 
    all_te_ls = set(ga_te_ls).union(gb_te_ls)
    all_te_superfamily_set = set([x.split("/")[0] for x in all_te_ls if len(x) != 0])
    return all_te_superfamily_set


def check_overlap(min1,max1,min2,max2):
    return max(0, min(max1, max2) - max(min1, min2))

def cnt_te(scaffold:str,start_pos:int,end_pos:int,gff_db):
    te_d = OrderedDict()
    g = gff_db.region(region=(scaffold,start_pos,end_pos))
    # g is a generator, containing gffutils.feature.Feature
    for f in g:
        # f.featuretype = family name
        superfamily = f.featuretype.split("/")[0]
        if superfamily not in te_d:
            te_d[superfamily] = check_overlap(min1=start_pos,min2=f.start,max1=end_pos,max2=f.end)
        else:
            te_d[superfamily] += check_overlap(min1=start_pos,min2=f.start,max1=end_pos,max2=f.end)
    return te_d

def plot_pie(col_map:list,superfamily_list:set):
    # prepare a universal color palatte
    p = cl.scales['12']['qual']['Paired']
    con_p = cl.interp(p, len(superfamily_list))
    color_d = dict(zip(superfamily_list,con_p))

    # create output folder
    if not os.path.exists("pie_plots"):
        os.makedirs("pie_plots")
    for p in col_map:
        print(f"Plotting {p[0]} {p[1]}...", file=sys.stderr)
        # p is already a list 
        # cnt the TE composition in species
        target_te_d = cnt_te(scaffold=p[3], start_pos=int(p[4]), end_pos=int(p[5]), gff_db=ga_db)
        query_te_d = cnt_te(scaffold=p[6], start_pos=int(p[7]), end_pos=int(p[8]), gff_db=gb_db)
        # plot
        # build titles
        target_title = f"Ampcit {p[3]}:{p[4]}-{p[5]}"
        query_title = f"ArcCen {p[6]}:{p[7]}-{p[8]}"
        main_title = f"qtl_id:{p[0]} collinearity_block:{p[1]}"
        # build color list for each pie
        target_color = [color_d[x] for x in target_te_d.keys()]
        query_color = [color_d[x] for x in query_te_d.keys()]
        
        fig = make_subplots(rows=1, cols=2, specs=[[{'type':'domain'}, {'type':'domain'}]], subplot_titles=[target_title,query_title])
        fig.add_trace(go.Pie(labels=list(target_te_d.keys()), values=list(target_te_d.values()), marker_colors=target_color),1,1)
        fig.add_trace(go.Pie(labels=list(query_te_d.keys()), values=list(query_te_d.values()), marker_colors=query_color),1,2)
        
        fig.update_layout(title_text=main_title)
        # change the subplot_titles font size, they are annotation, just change it manually
        for i in fig['layout']['annotations']:
            i['font'] = dict(size=10)
        fig.write_image(f"pie_plots/{p[0]}_{p[1]}.png")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--map", type = str, dest= "map",help = "qtl region map")
    parser.add_argument("-ga", "--gffa", type = str, dest= "ga",help = "reformated TE annotation (gff3) of target species")
    parser.add_argument("-gb", "--gffb", type = str, dest= "gb",help = "reformated TE annotation (gff3) of query species")
    parser.add_argument("-p","--pie",dest = "pie",action="store_true", help = "Pie Plot")
    parser.add_argument("-b","--bar",dest = "bar",action="store_true", help = "Bar Plot")
    args = parser.parse_args()

    col_map = read_collinearity_map(m_path=args.map)
    ga_db = gb_db = None
    load_gff_db(ga_path=args.ga, gb_path=args.gb)    
    all_te_superfamily_set = all_te_family(ga_path=args.ga, gb_path=args.gb)

    if args.pie:
        plot_pie(col_map=col_map,superfamily_list=all_te_superfamily_set)

    if args.bar:
        pass

