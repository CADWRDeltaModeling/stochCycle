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

plt.style.use(["seaborn-colorblind",'seaborn-paper'])
fig,((ax00,ax01,ax10,ax11)) = plt.subplots(4,1,sharex=True,
                                             figsize=(6,6.5))
#fig,((ax00,ax01)) = plt.subplots(2,1,sharex=True,
#                                             figsize=(7,7))

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

fname3='tropical_vmd/tropical_vmd_amp.csv'
dfvmd = pd.read_csv(fname3,skipinitialspace=True,sep=",",header=0).iloc[0:3200]
dfvmd.columns = ["vmd_"+x for x in dfvmd.columns] 
print("shape")
print(dfvmd.shape)
dfvmd.index = df.index
df = pd.concat((df,df1,df2,dfvmd),axis=1)
print(df.columns)
df.index.name="hours"

cases = [None,96,70]
mpl.rcParams['axes.prop_cycle'] = cycler(markevery=cases)



colors = ["black","black","0.6","black","black","0.3"]
styles = ['-',':','--',"-.","o","^",":"]
styles0 = ["-","-o","-^"]
df['dummy'] = df.D1
df['dummy'] = np.nan

df[["D1","SCHA_23_1",'dummy']].plot(ax=ax00,style=styles0,
                            color=colors,
                            markersize=6,markevery=96)
      
df[["SCHA_24_1"]].plot(ax=ax00,style=["-^"],
                       color='0.6',
                       markersize=6,markevery=144)

lines00 = ax00.get_lines()
legend=ax00.legend(lines00,["Analytical","SCHA(2,3)","SCHA (2,4)"],
           loc="upper right",markerscale=1)
# legend.legendHandles[1]._legmarker.set_markersize(5)
# legend.legendHandles[1]._legmarker.set_alpha(1)           
# legend.legendHandles[2]._legmarker.set_markersize(5)
# legend.legendHandles[2]._legmarker.set_alpha(1) 
           
#df[["D1","vmd_D1_amp"]].plot(ax=ax01,style=["-X"],
#                         color="0.6",markersize=5,markevery=32)
colors2=["black","black","black","black","0.6"]
styles2=["-",":","--","-.","-X"]
df[["D1","B2.5_D1","B4_D1","B8_D1","vmd_D1_amp"]].plot(ax=ax01,style=styles2,
                         color=colors2,markersize=6,markevery=96)


for line in ax00.get_lines():
    if line.get_label() == 'D1':
        line.set_linewidth(1.5)
        line.set_color("black")
    else:
        line.set_linewidth(1)



df[["D3","SCHA_23_3","dummy"]].plot(ax=ax10,style=styles0,
                         color=colors,markersize=6,markevery=96)
df[["SCHA_24_3"]].plot(ax=ax10,style=["-^"],
                         color="0.6",markersize=6,markevery=144)

colors2=["black","black","black","black","0.6"]
styles2=["-",":","--","-.","-X"]
df[["D3","B2.5_D3","B4_D3","B8_D3","vmd_D3_amp"]].plot(ax=ax11,style=styles2,
                         color=colors2,markersize=6,markevery=96)


for line in ax01.get_lines():
    if line.get_label() == 'D3':
        line.set_linewidth(1.5)
        line.set_color("black")
    else: 
        line.set_linewidth(1)
        
ax00.set_xlim(120,540)
ax00.set_xlabel("hours")
#ax00.legend(["Analytical","Morlet B=2.5","Morlet B=4",
#            "Morlet B=8","SCHA(2,3)","SCHA (2,4)","VMD"],
#           loc="upper right")

ax00.set_ylabel("D1 Amplitude")
ax00.set_ylim(0.,0.95)
ax01.legend(["Analytical","Morlet B=2.5","Morlet B=4",
             "Morlet B=8","VMD"],
           loc="upper right")
ax01.set_ylabel("D1 Amplitude")
ax01.set_ylim(0,0.95)
ax10.get_legend().remove()
ax10.set_ylim(0,0.4)
ax10.set_ylabel("D3 Amplitude")
ax11.get_legend().remove()
ax11.set_ylim(0,0.4)
ax11.set_ylabel("D3 Amplitude")
plt.tight_layout()
plt.show()
fig.savefig("tropical_fig.png")

