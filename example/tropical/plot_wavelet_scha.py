#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Created on Sat Jun  4 19:27:58 2022

@author: eli
"""
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt
from cycler import cycler
dfs=[]
for B in [2.5,4,8]:
    fname = f"cmor_{B}-0.96_result.csv"
    df0 =pd.read_csv(fname,index_col=0,header=0)
    if len(dfs)>0:
        cols=list(df0.columns)
        usecols = [col for col in cols if "_" in col ]
        df0 = df0[usecols]
    dfs.append(df0)

df = pd.concat(dfs,axis=1)
print(df)     

fname1 = "tropical_2_4/SCHA_2_4_param_amp.csv"

df1 = pd.read_csv(fname1,header=0,usecols=[1,2,3],skiprows=1).iloc[1:,]
df1.index=df.index
df1.columns=["SCHA_24_0","SCHA_24_1","SCHA_24_3"]
fname2 = "tropical_2_3/SCHA_2_3_param_amp.csv"
df2 = pd.read_csv(fname2,header=0,usecols=[1,2,3],skiprows=1).iloc[1:,]
df2.index=df.index
df2.columns=["SCHA_23_0","SCHA_23_1","SCHA_23_3"]

df = pd.concat((df,df1,df2),axis=1)
df.index.name="hours"

cases = [None,96,70]
mpl.rcParams['axes.prop_cycle'] = cycler(markevery=cases)

plt.style.use(["seaborn-colorblind",'seaborn-paper'])
fig,(ax0,ax1) = plt.subplots(2,sharex=True,figsize=(7,7))

colors = ["black","black","black","0.6","black","black"]
styles = ['-',':','--',"-.","o","^"]

df[["D1","B2.5_D1","B4_D1","B8_D1","SCHA_23_1"]].plot(ax=ax0,style=styles,
                         color=colors,markersize=6,markevery=96)
df[["SCHA_24_1"]].plot(ax=ax0,style=["^"],
                         color="black",markersize=7,markevery=48)

for line in ax0.get_lines():
    if line.get_label() == 'D1':
        line.set_linewidth(1.5)
        line.set_color("black")
    else:
        line.set_linewidth(1)


df[["D3","B2.5_D3","B4_D3","B8_D3","SCHA_23_3","SCHA_24_3"]].plot(ax=ax1,style=styles,
                                        color=colors,markersize=7,markevery=96)

df[["SCHA_24_3"]].plot(ax=ax1,style=["^"],
                         color="black",markersize=7,markevery=48)

for line in ax1.get_lines():
    if line.get_label() == 'D3':
        line.set_linewidth(1.5)
        line.set_color("black")
    else: 
        line.set_linewidth(1)
        
ax0.set_xlim(120,540)
ax0.set_xlabel("hours")
ax0.legend(["D1 Analytical","Morlet B=2.5","Morlet B=4",
            "Morlet B=8","SCHA(2,3)","SCHA (2,4)"],
           loc="upper right")
ax0.set_ylabel("|Amplitude|")
ax0.set_ylim(0.,1.)
ax1.set_ylabel("|Amplitude|")
ax1.set_ylim(0,0.4)
ax1.get_legend().remove()
plt.tight_layout()
plt.show()
fig.savefig("tropical_fig.png")

