import os
import sys
import getopt
from glob import glob
from datetime import datetime
from tempfile import TemporaryFile

env = os.environ.copy()
pjoin = os.path.join


gender="Female"
condition="Day3"
celltype="Cardiomyocytes"
step=1


celltype_dict = {
     "Macrophages" : ["Macrophages", "inflammatory macrophages", "Anti-inflammatory macrophages"],
     "Fibroblasts" : ["Fibroblasts", "Myofibroblasts"]
 }

## Parameters
try:
    options,args = getopt.getopt(sys.argv[1:], "t:c:s:x:")
except getopt.GetoptError:
    print("Erorr Parametes")
    sys.exit()
for k, v in options:
    if k in "-t":
        gender = v 
        print("gender: ", gender)

    if k in "-c":
        celltype = v 
        print("celltype: ", celltype)

    if k in "-x":
        condition = v
        print("condition: ", condition)

    if k in "-s": # step1 prepare, step2 differ
        step = int(v)
        print("step: ", step)
#end parameter


## run bash shell script
def get_out(*args):
    print("Executing:", " ".join(args))
    os.system(" ".join(args)) 
    print("Finished:", " ".join(args))


def samtools_merge_bam(indir,Dir, condition, gender, celltype, threads):
    CMD="samtools merge"
    ndir =  pjoin(Dir, celltype, "BAM")
    if not os.path.exists(ndir):
        os.makedirs(ndir)

    Genderbams=glob(f"{indir}/{condition}_{gender}_{celltype}.bam")

    if celltype in celltype_dict.keys():
        stackGenderbams = [glob(f"{indir}/{condition}_{gender}_{x}.bam") for x in celltype_dict[celltype]]
        Genderbams = [f"{item}" for sublist in stackGenderbams for item in sublist] 

    bams = Genderbams

    #if  len(bams) == 1:
    #    return

    outputbam = pjoin(ndir, f"{condition}_{gender}_{celltype}.bam")
    inputbams = " ".join([f"'{x}'" for x in bams])

    PARAM = f"\
            -f --thread {threads} \
            '{outputbam}' \
            {inputbams}"

    get_out(CMD, PARAM)

#endf 


def exist_cells(indir, gender, celltype):
    BAMdir = pjoin(indir, f"*_{gender}_{celltype}.bam")
    return os.path.exists(BAMdir)  
#endf



def samtools_index_bam(Dir,condition, gender, celltype, isOthers=False):
    #mkdir
    CMD="samtools index"
    input_file = pjoin(Dir,celltype, "BAM", f"{condition}_{gender}_{celltype}.bam")
    PARAM = f"'{input_file}'"
    get_out(CMD, PARAM) 

#endf 


def bigWig(Dir,condition, gender, celltype, isOthers=False, bs=24, p=24, norm="CPM"):
    #mkidir
    CMD="bamCoverage"
    ndir =  pjoin(Dir, celltype, "BigWig")
    if not os.path.exists(ndir):
        os.makedirs(ndir)

    inbam = pjoin(Dir, celltype, "BAM", f"{condition}_{gender}_{celltype}.bam")
    outbw = pjoin(Dir, celltype, "BigWig", f"{condition}_{gender}_{celltype}.bw") 

    PARAM=f"\
          -bs {bs} \
          -b '{inbam}' \
          -o '{outbw}' \
          -p {p} \
          --normalizeUsing {norm}"

    get_out(CMD, PARAM)
#endf 


def peak_calling(Dir, condition, gender, celltype, isOthers=False, spec="mm", q=0.05):
    CMD="macs2 callpeak" 
    ndir =  pjoin(Dir, celltype, "Peaks")
    if not os.path.exists(ndir):
        os.makedirs(ndir)
    peakdir = pjoin(Dir,celltype, "Peaks") 
    inbam = pjoin(Dir, celltype, "BAM", f"{condition}_{gender}_{celltype}.bam")
    n_ = f"{condition}_{gender}_{celltype}"
    PARAM=f"\
          -t '{inbam}' \
          -n '{n_}' \
          --outdir '{peakdir}' \
          -g {spec} \
          --nomodel \
          --nolambda \
          -f BAMPE \
          -q {q} \
          --keep-dup all"

    get_out(CMD, PARAM)
#endf

def del_chr(Dir, condition,gender, celltype, isOthers=False):
    CMD="sed -i"

    narrowpeak_file = pjoin(Dir,celltype, "Peaks", 
                           f"{condition}_{gender}_{celltype}_peaks.narrowPeak")
    PARAM=f"\
          '/random\|chrUn\|chrM/d' \
          '{narrowpeak_file}'"

    get_out(CMD, PARAM)
#endf 


def bed_intersect(Dir, condition,gender, celltype, blacklist_file, isOthers=False):
    CMD = "bedtools intersect"
    narrowpeak_file = pjoin(Dir, celltype,"Peaks", f"{condition}_{gender}_{celltype}_peaks.narrowPeak")
    outbed = pjoin(Dir,celltype, "Peaks", f"{condition}_{gender}_{celltype}_peaks_wo_blacklist.narrowPeak")
 
    PARAM = f"\
            -a '{narrowpeak_file}' \
            -b '{blacklist_file}' \
            -wa \
            -v \
            > '{outbed}'"


    get_out(CMD, PARAM)
#endf

def footprint(Dir, condition,gender, celltype, isOthers=False, organ="mm10"):
    #mkdir
    CMD="rgt-hint footprinting"
    ndir =  pjoin(Dir, celltype, "Footprints")
    if not os.path.exists(ndir):
        os.makedirs(ndir)

    inbam = pjoin(Dir,celltype, "BAM", f"{condition}_{gender}_{celltype}.bam")
    inpeak=pjoin(Dir,celltype, "Peaks", f"{condition}_{gender}_{celltype}_peaks_wo_blacklist.narrowPeak")
    outpre=f"{condition}_{gender}_{celltype}"

    PARAM=f"\
            --organism={organ} \
            --atac-seq \
            --paired-end \
            --output-location='{Dir}/{celltype}/Footprints' \
            --output-prefix='{outpre}' \
            '{inbam}' \
            '{inpeak}'"

    get_out(CMD, PARAM)
#endf

def match_motif(Dir, condition,gender, celltype, isOthers=False, organ="mm10"):
    #mkdir
    CMD = "rgt-motifanalysis matching"
    ndir =  pjoin(Dir, celltype, "MotifMatching")
    if not os.path.exists(ndir):
        os.makedirs(ndir)

    out_loc = pjoin(Dir,celltype, "MotifMatching")
    input_file =pjoin(Dir, celltype, "Footprints", f"{condition}_{gender}_{celltype}.bed")
   
    PARAM= f"\
            --organism={organ} \
            --output-location='{out_loc}' \
            --input-files='{input_file}'"
    get_out(CMD, PARAM)
#endf



def hint_diff_celltype_merge_gendervs(Dir, condition, celltype, gendervs_list, organ="mm10", nc=10):
    CMD="rgt-hint differential"

    gendervs = f"{gendervs_list[0]}-vs-{gendervs_list[1]}"
    ndir = pjoin(Dir,  f"{gendervs}_merge", condition, celltype)
    if not os.path.exists(ndir):
        os.makedirs(ndir)


    motif_list = []
    bam_list = []
    ct = celltype
    a_motif_list = glob(f"{Dir}/{ct}/MotifMatching/{condition}_{gendervs_list[0]}_*_mpbs.bed")
    b_motif_list = glob(f"{Dir}/{ct}/MotifMatching/{condition}_{gendervs_list[1]}_*_mpbs.bed")
    motif_list.extend(a_motif_list)
    motif_list.extend(b_motif_list)
    a_bam_list = glob(f"{Dir}/{ct}/BAM/{condition}_{gendervs_list[0]}*bam")
    b_bam_list = glob(f"{Dir}/{ct}/BAM/{condition}_{gendervs_list[1]}*bam")
    bam_list.extend(a_bam_list)
    bam_list.extend(b_bam_list)

    motif_list = ["'%s'" % st for st in motif_list]
    motif_list.sort()
    bam_list = ["'%s'" % st for st in bam_list]
    bam_list.sort()
    bams = ["%s" % st.split("/")[-1]  for st in bam_list ] ###BUG, fibroblasts
    cases = ["_".join(i.split("_")[0:2]) for i in bams]
    length = len(cases)
    if length <= 1:
        print("Too few cases", length)
        sys.exit(1)

    colors = ["#8dd3c7", "#bebada", "#fb8072",
              "#80b1d3", "#b3de69", "#fccde5",
              "#d9d9d9", "#bc80bd", "#ccebc5",
              "#ffed6f"]


    for i in range(len(motif_list)):
        bn_bam = os.path.basename(bam_list[i])
        bn_mot = os.path.basename(motif_list[i])
        ## same courterpart predix
        assert(os.path.splitext(bn_bam)[0] == os.path.splitext(bn_mot)[0][:-5])


    motif_param = ",".join(motif_list)
    bam_param = ",".join(bam_list)
    cases_param = ",".join(cases)
    color_param = ",".join(colors[0:length])


    PARAM=f"\
          --organism={organ} \
          --bc \
          --nc={nc} \
          --mpbs-files={motif_param} \
          --reads-files={bam_param} \
          --conditions={cases_param} \
          --output-location='{ndir}' \
          --output-prefix='{celltype}' \
          --color={color_param}"
    get_out(CMD, PARAM)
#endf



blacklist_file="/home/mingbo/mingbo/project/Mouse_scATAC//Blacklists/mm10-blacklist.v2.bed"
indir = "/home/mingbo/mingbo/project/Mouse_scATAC/HINT_myo/data/SplitBamFilesByClustering/BAM"
Dir=f"/home/mingbo/mingbo/project/Mouse_scATAC/HINT_myo/data/DiffMergeCelltypes_condition_gender"


print(f"checking gender: [{gender}], celltype: [{celltype}]", datetime.now())
#go = exist_cells(indir, gender, celltype)
#print("existed? ", go)
#if not go:
#    print("No cell exists")
#    sys.exit(1)

msg = f"{gender} {celltype}"
#"""

if step == 1:
    print("merge", msg, datetime.now())
    samtools_merge_bam(indir,Dir, condition, gender, celltype, 100)


    print("index", msg, datetime.now())
    samtools_index_bam(Dir,condition, gender, celltype)


    print("bw", msg, datetime.now())
    bigWig(Dir,condition, gender, celltype, bs=24, p=24, norm="CPM")


    print("peak calling", msg, datetime.now())
    peak_calling(Dir,condition, gender, celltype, spec="mm", q=0.05)



    print("del_chr", msg, datetime.now())
    del_chr(Dir,condition, gender, celltype)


    print("intersect", msg, datetime.now())
    bed_intersect(Dir,condition,gender, celltype, blacklist_file)


    print("footprint", msg, datetime.now())
    footprint(Dir,condition, gender, celltype, organ="mm10")



    print("motif", msg, datetime.now())
    match_motif(Dir,condition, gender, celltype, organ="mm10")

#"""
elif step == 2:
    print("hint diff merge celltypes gender vs", msg, datetime.now())
    gendervs_list = ["Female", "Male"]
    hint_diff_celltype_merge_gendervs(Dir,condition, celltype, gendervs_list, organ="mm10", nc=100)

