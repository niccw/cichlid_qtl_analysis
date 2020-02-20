#!/usr/bin/env python3
import sys
import re

def parse_gff_for_mcscanx(gff):
    print("## gff for mcscanx")
    with open(gff,"r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            col = line.strip().split("\t")
            if col[2] != "mRNA":
                continue
            # extract gene name from ID
            gene_id = re.search(r"(?<=ID=).*?(?=;)",col[8])[0]
            print(f"{col[0]}\t{gene_id}\t{col[3]}\t{col[4]}" ,file = sys.stdout)

if __name__ == "__main__":
    parse_gff_for_mcscanx(sys.argv[1])

