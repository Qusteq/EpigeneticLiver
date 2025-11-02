import re
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

def worker(f):
    mark = re.findall("Results_(.*?)/", f)[0]
    sample = re.findall("\/(\w*?)_Peak_Anno_knownGene.xls", f)[0]
    df = pd.read_csv(f, sep="\t")
    df["dist"] = abs(df["distanceToTSS"] / 1000)
    df["MACS peaks"] = f"{mark}.{sample}"
    return df

marks = ["H3K4me3", "H3K27ac", "H3K27me3"]
files = open("files.txt").readlines()
files = [f.strip() for f in files]
lst = []
for f in files:
    lst.append(worker(f))
df = pd.concat(lst)
df = df[df["dist"] < 10]
df = df.reset_index(drop=True)
plt.clf()
sns.kdeplot(data=df, x="dist", hue="MACS peaks", bw_adjust=1, palette="tab20", gridsize=200, hue_order=["H3K4me3.NASH", "H3K4me3.Normal", "H3K27ac.NASH", "H3K27ac.Normal", "H3K27me3.NASH", "H3K27me3.Normal"])
plt.xlim(-0.1, 10)
plt.axvline(x=2, color='black', linestyle='--')
plt.xlabel("Kb from the nearest gene")
plt.ylabel("Density")
plt.savefig("5.1.svg", dpi=300)
