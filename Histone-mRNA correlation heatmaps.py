import pandas as pd
f = "mRNA/gene_count.xls"
'''
Group.1 Ctrl1   Ctrl2   Ctrl3   NASH1   NASH2   NASH3
A1BG    56984   38691   38841   35089   26093   52558
A1CF    15660   16323   13824   19640   9904    13174
A2ML1   0       0       0       0       0       0
A3GALT2 11      14      24      0       11      14
'''
df = pd.read_csv(f, sep="\t")
df["NASH"] = df["NASH1"] + df["NASH2"] + df["NASH3"]
df = df[["Group.1", "NASH"]]

f2 = "GCF_012559485.2_MFA1912RKSv2.ncbiRefSeq.gene.bed"
'''
NC_052273.1 54989800    54996351    A1BG    .   -
NC_052263.1 82401893    82478985    A1CF    .   +
NC_052265.1 10407076    10462280    A2ML1   .   -
NC_052255.1 189974072   189980348   A3GALT2 .   +
'''
df2 = pd.read_csv(f2, sep="\t", header=None)
df = pd.merge(df, df2, left_on="Group.1", right_on=3)
df = df.sort_values("NASH", ascending=False)

f3 = "chr.txt"
'''
chr1    NC_052255.1
chr2    NC_052256.1
chr3    NC_052257.1
'''
df3 = pd.read_csv(f3, header=None, sep="\t")

df = df[df[0].isin(df3[1])]
df4 = df.copy()
df[[0, 1, 2, 3, 4, 5]].to_csv("sorted_gene.bed", sep="\t", index=False, header=False)

mat = pd.DataFrame()
files = ["H3K4me3.NASH.matrix.gz", "H3K27ac.NASH.matrix.gz", "H3K27me3.NASH.matrix.gz"]
for f in files:
    df = pd.read_csv(f, compression="gzip", comment="@", sep="\t", header=None)
    mark = f.split(".")[0]
    mat[mark] = df.iloc[:, 6:].mean(axis=1).values
mat["RNA"] = df4["NASH"].apply(lambda x: np.log2(x/3 + 1)).values
mat.to_csv("mat.csv", sep="\t", index=False)
