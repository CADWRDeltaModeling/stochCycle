# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 19:27:58 2022

@author: eli
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt

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


plt.style.use(["seaborn-colorblind",'seaborn-paper'])
fig,(ax0,ax1) = plt.subplots(2,sharex=True)


df[["D1","B2.5_D1","B4_D1","B8_D1"]].plot(ax=ax0,style=['-','--','-.',':'],
                         color=["black","0.3","black","0.3"])

for line in ax0.get_lines():
    if line.get_label() == 'D1':
        line.set_linewidth(1)
        line.set_color("black")
    else:
        line.set_linewidth(1)
ax0.set_xlim(100,600)

df[["D3","B2.5_D3","B4_D3","B8_D3"]].plot(ax=ax1,style=['-','--','-.',':'],
                                        color=["black","0.3","black","0.6"])
#ax1.get_legend().remove()
for line in ax1.get_lines():
    if line.get_label() == 'D3':
        line.set_linewidth(1)
        line.set_color("black")
    else: 
        print("ddd")
        line.set_linewidth(1)
ax0.set_xlim(100,600)
ax0.set_ylabel("Amplitude")
ax1.set_ylabel("Amplitude")

plt.show()