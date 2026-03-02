#!/usr/bin/env python
# coding: utf-8

# In[13]:


import pandas as pd
import sys
filename=sys.argv[1]
length_file=sys.argv[2]
out=filename + "_with_TPM"
df=pd.read_table(filename,sep="\t")
length=pd.read_table(length_file,sep="\t")
name_list=list(df[:0])[1:]
df= pd.concat([df,length["Length"]],axis=1)
def count_per_length(df):
    return float(df[name])/float(df["Length"])
def tpm_cal(df):
    return (float(df[counts_per_length])/float(cpt_sum)) * 1000000
df_out = pd.DataFrame()
df_out["Geneid"] = df["Geneid"]
for name in name_list:
    counts_per_length=name+"counts_per_length"
    count_sum = df[name].sum()
    df[counts_per_length] = df.apply(count_per_length, axis=1)
    cpt_sum = df[counts_per_length].sum()
    tpm = name + "TPM"
    df_out[tpm] = df.apply(tpm_cal, axis = 1)
df_out.to_csv(out, sep="\t", index=False)

