from collections import defaultdict 
from pathlib import Path
import getopt
import sys

"""
python detail_corr.py
python detail_corr.py -f save/corr_pval_final_0.1.tsv
python detail_corr.py -f save/corr_pval_final_0.2.tsv
python detail_corr.py -f save/corr_pval_final_topn_10000.tsv
python detail_corr.py -f save/corr_pval_final.tsv
python detail_corr.py -f save/corr_pval_filtered.tsv
python detail_corr.py -f save/corr_pval.tsv
"""

corr = "0.6"
pval = "0.05"
TimePoint = "Healthy"


try:
    options,args = getopt.getopt(sys.argv[1:],"c:p:t:")
except getopt.GetoptError:
    print("Erorr Parametes")
    sys.exit()
for name,value in options:
    if name in "-c":
        corr = value
        print("", corr) 
    if name in "-p":
        pval = value
        print("", pval) 
    if name in "-t":
        TimePoint = value
        print("On:", TimePoint) 


fname = f"save/{TimePoint}/corr_pval_final_{corr}_{pval}.tsv"


p = Path(fname)
if not p.exists():
    print("No this file")
    sys.exit()



out_file = p.parent / f"detail_{p.name}" 

fname_TimePoint = f"save/{TimePoint}/peak2gene_putative.tsv" 
f_TimePoint = open(fname_TimePoint)
TimePoint_first_line = f_TimePoint.readline()


f = open(fname)
first_line = f.readline()
other_header = "\t".join(first_line.strip("\n").split("\t")[3:]) 

fw = open(out_file, "w")
fw.write("%s\t%s\n" % (TimePoint_first_line.strip("\n"), other_header))


d = defaultdict(str)
for line in f_TimePoint:
    items_TimePoint = line.strip("\n").split("\t")
    peak_gene_TimePoint = items_TimePoint[0:2]
    d[tuple(peak_gene_TimePoint)] = "\t".join(items_TimePoint)



for line in f:
    items = line.strip("\n").split("\t")
    other_info = "\t".join(items[3:])
    peak_gene = items[0:2]
    
    out = d.get(tuple(peak_gene), "")
    if not out:
        continue            
    fw.write("%s\t%s\n" % (out, other_info))
fw.close()        

