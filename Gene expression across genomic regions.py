import os
import pandas as pd
import random
import numpy as np

def sample(out_f = "./random.bed"):
    f1 = "MFA1912RKSv2.chrom.sizes"
    df = pd.read_csv("./results/SE.bed", sep="\t", header=None, on_bad_lines='skip')
    df = df.astype({1:int, 2:int})
    df["d"] = df[2] - df[1]
    l = np.median(df["d"])
    n = df.shape[0]
    cmd = f"bedtools random -n {n} -l {l} -g {f1} > {out_f}"
    os.system(cmd)

for i in range(100):
    out_f = f"./results/random_{i}.bed"
    sample(out_f)
