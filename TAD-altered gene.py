import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns

df = pd.read_csv("NASH_vs_Normal_raw_DEG.geneSymbol.xls", sep="\t")

resolution = 10000
expand = [i.strip() for i in open(f"expand.{resolution}.genes").readlines()]  
stable = [i.strip() for i in open(f"stable.{resolution}.genes").readlines()]  
shift = [i.strip() for i in open(f"shift.{resolution}.genes").readlines()]
shrink = [i.strip() for i in open(f"shrink.{resolution}.genes").readlines()]

df["Change"] = "Stable"
df.loc[df["GeneID"].isin(expand), "Change"] = "Expand"
df.loc[df["GeneID"].isin(shift), "Change"] = "Shift"
df.loc[df["GeneID"].isin(shrink), "Change"] = "Shrink"
df["Expand"] = 0
df["Stable"] = 0
df["Shift"] = 0
df["Shrink"] = 0
df.loc[df["GeneID"].isin(expand), "Expand"] = 1
df.loc[df["GeneID"].isin(stable), "Stable"] = 1
df.loc[df["GeneID"].isin(shift), "Shift"] = 1
df.loc[df["GeneID"].isin(shrink), "Shrink"] = 1
df["Change2"] = df["Expand"] + df["Stable"] + df["Shift"] + df["Shrink"]
df = df[df["Change2"] == 1]

df["size"] = 3
df.loc[df["GeneID"].isin(expand), "size"] = 15
df.loc[df["GeneID"].isin(shift), "size"] = 15
df.loc[df["GeneID"].isin(shrink), "size"] = 15
df['neg_log10_pvalue'] = -np.log10(df['pvalue'])
pval_thresh_lst = [0.01, 0.05]
fold_thresh_lst = [1]
df_ori = df.copy()
for pval_thresh in pval_thresh_lst:
	for fold_thresh in fold_thresh_lst:
		df = df_ori.copy()
		df.loc[df["padj"] > pval_thresh, "Change"] = "non-significant"
		df.loc[(df["log2FoldChange"] < 0) & (df["log2FoldChange"] > -fold_thresh), "Change"] = "non-significant"
		df.loc[(df["log2FoldChange"] > 0) & (df["log2FoldChange"] < fold_thresh), "Change"] = "non-significant"
		df.to_csv(f"fc{fold_thresh}_p{pval_thresh}.csv", sep="\t", index=False)
plt.clf()
plt.figure(figsize=(6, 4))
sns.scatterplot(data=df, x="log2FoldChange", y="neg_log10_pvalue", hue="Change", palette={"Stable": "#E7B800", "Expand":"#FC4E07", "Shift":"blue", "Shrink": "#00AFBB", "non-significant": "#8491B4FF"}, alpha=0.7, s=13) 
plt.axhline(y=-np.log10(0.01), color='black', linestyle='--')  
plt.axvline(x=1.5, color='black', linestyle='--')  
plt.axvline(x=-1.5, color='black', linestyle='--')  
plt.xlabel('log2 Fold Change')
plt.ylabel('-log10 p-value')
plt.ylim(0, 40)
plt.xlim(-10, 10)
plt.title('NASH vs Ctrl')
plt.savefig("8.7.pdf")
