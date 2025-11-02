import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import os
import sys
from scipy import stats
import itertools
from statannotations.Annotator import Annotator

df_gene = pd.read_csv("gene_TPM.xls", sep="\t")
df_gene["NASH"] = (df_gene["NASH1"] + df_gene["NASH2"] + df_gene["NASH3"]) / 3
df_gene["CTRL"] = (df_gene["Ctrl1"] + df_gene["Ctrl2"] + df_gene["Ctrl3"]) / 3
df_dist = pd.read_csv(f"enhancer.noloop.nearTSS.bed", sep="\t", header=None)
df_dist["dist"] = df_dist[16] / 1000
df_dist = df_dist[df_dist["dist"] <= 1000]
df_dist["dist_group"] = pd.cut(x=df_dist["dist"], bins=[0, 10, 50, 100, 200, 1000], labels=["<10Kb", "10Kb-50Kb", "50Kb-100Kb", "100Kb-200Kb", "200Kb-1M"], include_lowest=True)
df = pd.merge(df_gene, df_dist, left_on="GeneID", right_on=3)
pairs = [("<10Kb", "10Kb-50Kb"), ("<10Kb", "50Kb-100Kb"), ("<10Kb", "100Kb-200Kb"), ("<10Kb", "200Kb-1M"), ("<10Kb", ">1M")]
pairs = itertools.combinations(set(df_dist["dist_group"]), 2)
groups = list(set(df_dist["dist_group"]))
labels=["<10Kb", "10Kb-50Kb", "50Kb-100Kb", "100Kb-200Kb", "200Kb-1M", ">1M"]
ref = sorted(groups, key=lambda x: labels.index(x))[0]
others = set(groups) - set([ref])
pairs = []
for o in others:
    pairs.append((ref, o))
plt.clf()
ax = sns.boxplot(data=df, x="dist_group", y="NASH", showfliers=False)
annotator = Annotator(ax, pairs, data=df, x="dist_group", y="NASH") 
annotator.configure(test='Mann-Whitney', text_format='simple', loc='outside')
annotator.apply_and_annotate()
plt.xlabel("")
plt.ylabel("Gene expression\n(log2(TPM))")
plt.ylim(0, 125)
plt.savefig(f"8.19.enhancer.pdf")

df = df.rename({"dist": "dist(Kb)"}, axis=1)
df.to_csv("gene.noloop_enhancer.dist_group.csv", sep="\t", index=False)
