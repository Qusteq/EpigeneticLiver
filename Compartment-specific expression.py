import os
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import sys

def worker(f1, f2, sample="NASH"):
    resolution = f1.split("/")[-1].split(".")[1]
    df1 = pd.read_csv(f1, sep="\t")
    df2 = pd.read_csv(f2, sep="\t")
    df2.columns = [i + ".1" for i in list(df2.columns)]
    df = pd.concat([df1, df2], axis=1)
    df = df.dropna(subset=["E1", "E1.1"])
    df = df[df["E1"]!=0]
    df = df[df["E1.1"]!=0]

    df["compartment_type_1"] = "A"
    df.loc[df["E1"] < 0, "compartment_type_1"] = "B"
    df["compartment_type_2"] = "A"
    df.loc[df["E1.1"] < 0, "compartment_type_2"] = "B"

    df.loc[(df["compartment_type_1"] == "A") & (df["compartment_type_2"] == "A"), "switch_type"] = "Stable A"
    df.loc[(df["compartment_type_1"] == "B") & (df["compartment_type_2"] == "B"), "switch_type"] = "Stable B"
    df.loc[(df["compartment_type_1"] == "A") & (df["compartment_type_2"] == "B"), "switch_type"] = "A-to-B"
    df.loc[(df["compartment_type_1"] == "B") & (df["compartment_type_2"] == "A"), "switch_type"] = "B-to-A"
  
    df["sample"] = sample
    df = df.astype({"start": str, "end": str})
    df["compartment_id"] = df["chrom"] + "-" + df["start"] + "-" + df["end"]
    return df

resolution = sys.argv[1] #10000   
d = "compartment"
f1 = f"Normal.{resolution}.cis.vecs.tsv"
samples = ["NASH-MboI-R1-filtered", "NASH-MboI-R2-filtered"] 
f2s = [f'{d}/{sample}.{resolution}.cis.vecs.tsv' for sample in samples] 
df_lst = []
for f2 in f2s:
    tmp = f2.split("/")[-1]
    try:
        sample = tmp.split("-")[0] + tmp.split("-")[2][1:]
    except:
        sample = tmp.split(".")[0]
    df_tmp = worker(f1, f2, sample)
    df_lst.append(df_tmp)
df = pd.merge(df_lst[0], df_lst[1], left_on="compartment_id", right_on="compartment_id")
df = df[df["switch_type_y"] == df["switch_type_x"]]
df.to_csv(f"compartment_switch.{sample}.{resolution}.csv", sep="\t", header=None, index=False)
