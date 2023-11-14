import numpy as np
import pandas as pd

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams

rcParams['font.size'] = 16
rcParams['figure.figsize'] = (18, 9)

fig = plt.figure()
ax = fig.add_subplot(1, 1, 1)
colors = ["#9eeb30","#b3734d","#888888"]
labels = ["hypermethylated sites in YPK", "hypermethylated sites in OPK", "n.s."]

df = pd.read_table("./split_cg_context/methylDiff_CG.txt")
down = df[(df["qvalue"] < 10**-5) \
          & (df["methdiff"] < -25)]
up = df[(df["qvalue"] < 10**-5) \
        & (df["methdiff"] > 25)]
non_sig = df[(df["qvalue"] >= 10**-5) \
             | ((df["methdiff"] <= 25) & (df["methdiff"] >= -25))]

for i, met in enumerate([down, up, non_sig]):
    ax.scatter(met["methdiff"],
               -np.log10(met["qvalue"]),
               label=labels[i],
               color=colors[i],
               alpha=0.5)

ax.set_xlabel("Difference of methylation level (%)")

ax.set_ylabel("$-\log_10($q-value$)$")

ax.legend(loc="upper left", bbox_to_anchor=(1, 1))

filename = "volcano_compress_new"
out_dir = "."
for fmt in ["pdf","png"]:
    plt.savefig("{}/{}.{}".format(out_dir,
                                  filename,
                                  fmt),
                format=fmt,
                bbox_inches="tight")