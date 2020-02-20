#!/usr/bin/env python3
import sys

def parse_qtl_region(file_path:str, genome_file:str):
    # parse the genome file as dict
    genome_d = {}
    with open(genome_file,"r") as f:
        for line in f:
            scaffold, s_len = line.strip().split("\t")
            genome_d[scaffold] = s_len

            
    with open(file_path,"r") as f:
        # save the header
        header =  f.readline().strip().split("\t")
        print(*header, sep="\t", file=sys.stdout)

        traits = ""
        lg = ""
        qtl_region = ""

        for line in f:
            if not line.strip():
                continue
            col = line.strip().split("\t")

            if len(col) > 3: #  A QTL region
                traits = col[0]
                lg = col[1]
                qtl_region = col[2]
                s_chr = col[3]
                s_start = col[4]
                s_end = col[5]
                if s_end == "end":
                    # find the corresponding end of the scaffold
                    s_end = genome_d[s_chr]
            else: # sub record of a QTL region, we have to fill the Traits, LG and QTL_region
                try:
                    s_chr = col[0]
                    s_start = col[1]
                    s_end = col[2]
                except Exception as e:
                    print(line)
                    print(e)
                if s_end == "end":
                    s_end = genome_d[s_chr]
            
            # print the parsed record
            print(f"{traits}\t{lg}\t{qtl_region}\t{s_chr}\t{s_start}\t{s_end}", file = sys.stdout)

if __name__ == "__main__":
    # hard code the absoulte path...
    parse_qtl_region(file_path = "/scratch/nicola/other/cichlids/qtl_region/qtl_regions.txt",
    genome_file= "/proj/nicola/raw/genome/cichlid/AmpCit.dna.fa.genome")