import pandas as pd
import os
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

def load_deg():
    log2FC_threshold = 1
    padj_threshold = 0.05
    deg_file = "NASH_vs_Normal_raw_DEG.geneSymbol.xls"
    deg = pd.read_csv(deg_file, sep="\t") 
    deg["rna_Event"] = "Unchange"
    deg.loc[(deg["log2FoldChange"] > log2FC_threshold) & (deg["padj"] < padj_threshold), "rna_Event"] = "Up"
    deg.loc[(deg["log2FoldChange"] < log2FC_threshold) & (deg["padj"] < padj_threshold), "rna_Event"] = "Down"
    deg = deg[deg["rna_Event"].isin(["Up", "Down"])] 
    return deg

df = load_deg()
total = df.shape[0]
d = "../8.32/results/"
lst = []
for f in os.listdir(d):
    if f.endswith("txt"):
        prefix = f.split(".")[0]
        if prefix.startswith('CRE'):
            continue
        if prefix.startswith("random"):
            random = True
        else:
            random = False
        f = f"{d}/{f}"
        df_tmp = pd.read_csv(f, sep="\t", header=None)
        val = df[df["GeneID"].isin(df_tmp[0])]["NASH"].shape[0] / total * 100 
        if random:
            lst.append({"tag": "random", "val":val})
        else:
            lst.append({"tag": prefix, "val":val})
f = "SE.gene.txt"
df_tmp = pd.read_csv(f, sep="\t", header=None)
df = df[df["GeneID"].isin(df_tmp[0])]
df = df.drop('geneSymbol', axis=1)
df.to_csv("DE.overlap_SE.csv", sep="\t", index=False)
df = pd.DataFrame(lst)    
plt.clf()
sns.boxplot(data=df, x="tag", y="val", order=["SE", "random"]) 
plt.xlabel("")
plt.ylabel("%Gene\noverlapping")
plt.savefig("8.33.pdf")
