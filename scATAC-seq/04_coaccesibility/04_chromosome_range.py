from collections import defaultdict

#f = open("data/peak_annotation_expand_newformat.tsv")
import getopt
import os
import sys


TimePoint = "Healthy"

try:
    options,args = getopt.getopt(sys.argv[1:],"t:")
except getopt.GetoptError:
    print("Erorr Parametes")
    sys.exit()
for name,value in options:
    if name in "-t":
        TimePoint = value

f = open(f"save/{TimePoint}/corr.tsv")
f.readline()
last_chromosome = "xxx"

d = defaultdict(list)
genes_dict = defaultdict(set)

for idx, line in enumerate(f):
    items = line.strip("\n").split("\t")
    peak = items[1]
    gene = items[2]
    chromosome = peak.split("_")[0]
    genes_dict[chromosome].add(gene)
    if last_chromosome == "xxx":
        d[chromosome].append(idx+1)
        last_chromosome = chromosome
        continue
        
    if last_chromosome != chromosome and last_chromosome != "xxx":
        d[last_chromosome].append(idx)
        d[chromosome].append(idx+1)
        last_chromosome = chromosome


d[chromosome].append(idx)


fw = open(f"save/{TimePoint}/chromosome_range.txt", "w")
fw.write("chromosome start end genes\n")
for k, v in d.items():
    genes = genes_dict[k]
    fw.write("%s %d %d %s\n" %(k, v[0], v[1], ",".join(genes)))
fw.close()
