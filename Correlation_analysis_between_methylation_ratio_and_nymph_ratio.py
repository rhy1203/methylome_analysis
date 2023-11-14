import csv
import pandas as pd
import numpy as np

file = open(./split_cg_context/methylBase_CG_SK.txt, "r")
f = csv.reader(file, delimiter="\t", skipinitialspace=True)

oPK1_nymph_ratio = float(0.217)
oPK2_nymph_ratio = float(0.192)
yPK1_nymph_ratio = float(0.459)
yPK2_nymph_ratio = float(0.355)
SK1_nymph_ratio = float(0.302)
SK2_nymph_ratio = float(0.219)

for line in f:
    if len(line) < 10:
            continue

    Scaffold = line[0]
    Locus = line[1]
    oPK1_meth = float(line[5])/float(line[4])
    oPK2_meth = float(line[8])/float(line[7])
    yPK1_meth = float(line[20])/float(line[19])
    yPK2_meth = float(line[23])/float(line[22])
    nymph1_meth = float(line[35])/float(line[34])
    nymph2_meth = float(line[38])/float(line[37])
    meth = [oPK1_meth, oPK2_meth, yPK1_meth, yPK2_meth, SK1_meth, SK2_meth]
    nymph_ratio = [oPK1_nymph_ratio, oPK2_nymph_ratio, yPK1_nymph_ratio, yPK2_nymph_ratio, SK1_nymph_ratio, SK2_nymph_ratio]
    meth_s = pd.Series(meth)
    nymph_ratio_s = pd.Series(nymph_ratio)
    res = meth_s.corr(nymph_ratio_s)
    res_a = float(res)*float(res)

file.close()