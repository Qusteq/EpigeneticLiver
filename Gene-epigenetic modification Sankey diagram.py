import os
import numpy as np
import pandas as pd

# DEG
log2FC_threshold = 1
padj_threshold = 0.05
deg_file = "NASH_vs_Normal_raw_DEG.geneSymbol.xls"
deg = pd.read_csv(deg_file, sep="\t") 
deg["rna_Event"] = "Unchange" 
deg.loc[(deg["log2FoldChange"] > log2FC_threshold) & (deg["padj"] < padj_threshold), "rna_Event"] = "Up"
deg.loc[(deg["log2FoldChange"] < -log2FC_threshold) & (deg["padj"] < padj_threshold), "rna_Event"] = "Down"
df_deg = deg["rna_Event"].to_frame()
df_deg.index = deg["GeneID"]
df_deg = df_deg[df_deg["rna_Event"].isin(["Up", "Down"])]
df_deg.columns = ["RNA"]

# Meth
meth_file = "NASH_vs_Ctrl_DMR_Filt_Modified_Genes.xls"
meth = pd.read_csv(meth_file, sep="\t", index_col=0)
gain = []
loss = []
ambiguous = []
records = []
for gene, df_tmp in meth.groupby("GeneID"):
    tmp = df_tmp
    tmp2 = df_tmp.dropna()
    n1 = tmp.shape[0]
    n2 = tmp2.shape[0]
    events = set(tmp2["variation"])
    if n2 != 0:
        if ("Hypo" in events) and ("Hyper" in events):
            ambiguous.append(gene) 
            records.append(["Ambiguous", gene])
        elif "Hypo" in events:
            loss.append(gene)
            records.append(["Down", gene])
        elif "Hyper" in events:
            gain.append(gene)
            records.append(["Up", gene])

df2 = pd.DataFrame(records, columns=["Methylation", "Gene"])
df_meth_diff = df2[df2["Methylation"].isin(["Up", "Down"])]
df_meth_diff.set_index('Gene', inplace=True)
def chip_parser(mark, padj_threshold=0.05, log2FC_threshold=1):
    d = "2_3_4"
    annot_file = f"{d}/{mark}.diff.peaks.annot"
    peak_file = f"{d}/{mark}.diff.peaks"

    print(f"Processing {peak_file} {annot_file} ...")

    df = pd.read_csv(peak_file, sep="\t", comment="#")
    df['PeakID'] = df.index.to_series().apply(lambda x: f"{mark}_peak_{x + 1}")
    peak = df[["PeakID", "Event", "padj", "log2FC"]]
    peak["log2FC_abs"] = peak["log2FC"].apply(lambda x: np.abs(x))
    peak.loc[(peak["padj"] > padj_threshold) | (peak["log2FC_abs"] < log2FC_threshold), "Event"] = "Unchange"
    peak.loc[peak["Event"] == "Up", "Event"] = "Gain"
    peak.loc[peak["Event"] == "Down", "Event"] = "Loss"

    annot = pd.read_csv(annot_file, sep="\t")
    cols = ['PeakID'] + list(annot.columns)[1:]
    annot.columns = cols
    annot = annot[["PeakID", 'Gene Name', 'Distance to TSS', 'Nearest PromoterID']]
    annot["Distance to TSS"] = annot["Distance to TSS"].apply(lambda x: np.abs(x))
    annot = annot[annot["Distance to TSS"] < 2000]

    df = pd.merge(annot, peak, left_on="PeakID", right_on="PeakID") 

    unchange = []
    gain = []
    loss = []
    ambiguous = []
    records = []
    for gene, df_tmp in df.groupby("Gene Name"):
        tmp = df_tmp
        tmp2 = df_tmp.dropna()
        n1 = tmp.shape[0]
        n2 = tmp2.shape[0]
        events = set(tmp2["Event"])
        if n2 != 0:
            if ("Gain" in events) and ("Loss" in events):
                ambiguous.append(gene) 
                records.append(["Ambiguous", gene])
            elif "Loss" in events:
                loss.append(gene)
                records.append(["Down", gene])
            elif "Gain" in events:
                gain.append(gene)
                records.append(["Up", gene])
            else:
                unchange.append(gene)
                records.append(["Unchange", gene])

    df2 = pd.DataFrame(records, columns=[f"{mark}", "Gene"])
    df_mark_diff = df2[df2[f"{mark}"].isin(["Down", "Up"])]
    df_mark_diff.set_index('Gene', inplace=True)
    return df_mark_diff

df = pd.merge(df_deg, df_meth_diff, left_index=True, right_index=True, how="left")
marks = ["H3K4me3", "H3K27ac", "H3K27me3"]
for mark in marks:
    df_mark_diff = chip_parser(mark)
    df = pd.merge(df, df_mark_diff, left_index=True, right_index=True, how="left")
df["gene"] = df.index
df = df.melt(id_vars=["gene", "RNA"]).dropna()
df2 = pd.DataFrame()
df2["source"] = "RNA_" + df["RNA"]
df2["target"] = df["variable"] + "_" + df["value"] 
node_mapping = {node: idx for idx, node in enumerate(pd.concat([df2['source'], df2['target']], axis=0).unique())}
mapping_df = pd.DataFrame(list(node_mapping.items()), columns=['Node', 'Number'])
encoded_edges = df2.copy()
encoded_edges['source'] = df2['source'].map(node_mapping)
encoded_edges['target'] = df2['target'].map(node_mapping)
encoded_edges['value'] = 1
encoded_edges.to_csv("edges.txt", sep="\t", index=False)
mapping_df.to_csv("nodes.txt", sep="\t", index=False)

