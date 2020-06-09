
targetfile = 'list.txt'
"""
ENSMUSG00000075602      A       TGCGTCTATTTAGTGGAGCC    ACTG    Ly6a
ENSMUSG00000005672      C       TGCGTCTATTTAGTGGAGCC    GATC    Kit
ENSMUSG00000015316      T       TGCGTCTATTTAGTGGAGCC    CGAT    Slamf1
ENSMUSG00000015355      G       TGCGTCTATTTAGTGGAGCC    TCGA    Cd48
ENSMUSG00000042817      A       TGCGTCTATTTAGTGGAGCC    ATGC    Flt3
ENSMUSG00000042817      C       TGCGTCTATTTAGTGGAGCC    TCAG    Flk2
ENSMUSG00000034664      T       TGCGTCTATTTAGTGGAGCC    GACT    Itga2b
ENSMUSG00000001930      G       TGCGTCTATTTAGTGGAGCC    GTCA    Vwf
ENSMUSG00000022847      A       TGCGTCTATTTAGTGGAGCC    CATG    Thpo
"""

outputprefix = 'BM'
# bowtie2 index files, generate by modified cDNA fasta
index = '/mnt/data8/zhaoyuancun/GRCm38/GRCm38_idx'
# OverlapModeVal
OverlapModeVal = False
# verbose
verbocity = False

# The minimum allowed probe length
l = 36
# The maximum allowed probe length
L = 41
# The minimum allowed percent G + C
gcPercent = 20
# The maximum allowed percent G + C
GCPercent = 80
# The minimum allowed Tm
tm = 54
# The maximum allowed Tm
TM = 60
# Prohibited sequence list
X = 'AAAAA,TTTTT,CCCCC,GGGGG'
# The mM Na+ concentration
sal = 390
# The percent formamide being used
form = 50
# The minimum spacing between adjacent probes
sp = 0
# probelength
probelength = 70
# Shannon entropy to evalute probe
entropy = 1.0
# The temperature at which you want to hybridize your probes
hytemp = 37.0
# mfold mode 
mfold = False
# thread
thread = 4
# Accept detG value filtering for secondary structure check using mfold
detG = 0.0
# cDNA mode, if true, padlock probes will be hybridized to cDNA instead of RNA
cDNA = False
# vargibbs parameters
vargibbs = '/mnt/data8/zhaoyuancun/DNA-RNA_hybrids/VarGibbs-2.2/vargibbs'
par = '/mnt/data8/zhaoyuancun/DNA-RNA_hybrids/VarGibbs-2.2/data/AOP-DRFT.par'
saltscheme = 'chen13eq21'
ct = 0.3
# Concentration of higher concentration strand [nM] -typically the probe- to use for thermodynamic calculations
conc1 = 25
# Concentration of lower concentration strand [nM] -typically the probe- to use for thermodynamic calculations
conc2 = 25
nn_table = 'DNA_NN3'
