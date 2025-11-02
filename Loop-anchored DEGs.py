import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

loop_annot = "loops.all.with_annot.bedpe"
genes = pd.read_csv(loop_annot, header=None, sep="\t").iloc[:, -3]

df = pd.read_csv("NASH_vs_Normal_raw_DEG.geneSymbol.xls", sep="\t")

fc_thr = 1.5
p_thr = 0.05
df["Group"] = 'Unchaged'
df.loc[(df["log2FoldChange"] > fc_thr) & (df["pvalue"]<p_thr), "Group"] = "Up"
df.loc[(df["log2FoldChange"] < -fc_thr) & (df["pvalue"]<p_thr), "Group"] = "Down"

df['neg_log10_pvalue'] = -np.log10(df['pvalue'])

df = df[df["GeneID"].isin(genes)]

df.to_csv("loop_gene.DEG.csv", sep="\t", index=False)
plt.clf()
plt.figure(figsize=(6, 4))
sns.scatterplot(data=df, x="log2FoldChange", y="neg_log10_pvalue", hue="Group", palette={"Unchaged": "gray", "Up":"red", "Down": "blue"}, alpha=0.5, s=10) #, size="size")
plt.axhline(y=-np.log10(p_thr), color='black', linestyle='--')  
plt.axvline(x=fc_thr, color='black', linestyle='--')  
plt.axvline(x=-fc_thr, color='black', linestyle='--')  
plt.xlabel('-log2(FoldChange)')
plt.ylabel('-log10(P value)')
plt.ylim(0, 30)
n1 = df[df["Group"] == "Down"].shape[0]
n2 = df[df["Group"] == "Up"].shape[0]
plt.text(-5, 15, f"Downregulated\n{n1}", ha="center")
plt.text(7.5, 15, f"Upregulated\n{n2}", ha="center")
plt.savefig("8.13.pdf")
