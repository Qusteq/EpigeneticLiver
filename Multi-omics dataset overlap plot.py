import itertools
import pandas as pd
import os

marks = ["H3K4me3", "H3K27ac", "H3K27me3", "Methylation"]
files = [f"{mark}.diff.peaks" for mark in marks[:-1]]
files += ["NASH_vs_Ctrl_Filt_pval.0.01.xls"]
files2 = []
log2FC_threshold = 1
padj_threshold = 0.05
n = 0
for f, mark in zip(files, marks):
    n += 1
    out_f = f"{mark}.clean.peaks"
    cmd = f"awk '/^N/' {f} > {out_f}"
    df = pd.read_csv(f, sep="\t")
    if n <= 3:
        df = df[(abs(df["log2FC"]) > log2FC_threshold) & (df["padj"] < padj_threshold)]
        df = df[["Chrom", "Start", "End", "Id", "log2FC", "strand"]] 
    else:
        df["Id"] = df.index.to_series().apply(lambda x: f'Methylation_peak_{x+1}')
        df = df[["chr", "start", "end", "Id", "methyl_diff", "strand"]]
    df.to_csv(out_f, sep="\t", header=False, index=False)
    files2.append(out_f)
dic = dict(zip(marks, files2))
for mark1, mark2 in list(itertools.combinations(marks, 2)):
    f1 = dic[mark1]
    f2 = dic[mark2]
    cmd = f"bedtools intersect -a {f1} -b {f2} -wa -wb > {mark1}_{mark2}.intersect.bed"
    os.system(cmd)

dic = {}
for mark in marks:
    dic[mark] = []
for f in os.listdir("./"):
    if f.endswith(".intersect.bed"):
        mark1 = f.split(".")[0].split("_")[0] 
        mark2 = f.split(".")[0].split("_")[1] 
        df = pd.read_csv(f, sep="\t", header=None)
        dic[mark1] += list(df[3].values)
        dic[mark1] += list(df[9].values)
        dic[mark2] += list(df[3].values)
        dic[mark2] += list(df[9].values)

with open('upset_data.json', 'w') as json_file:
    json.dump(dic, json_file, indent=4)
