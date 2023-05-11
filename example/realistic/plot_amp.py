#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  4 19:27:58 2022

@author: eli
"""
import os
import matplotlib.pyplot as plt
from matplotlib import gridspec
from cycler import cycler
import numpy as np
import pandas as pd
import pywt
from vtools import datetime_elapsed
import seaborn as sns

from dms_datastore.read_ts import *
CFS2CMS = 0.028316847
item = "amp"

cfreq = [0.2529340,0.5058680,0.7588021,1.0117361,1.2646701,1.5176041]

plt.style.use(["seaborn-deep","seaborn-paper"])

# File with time units of hours (but 15min time step)
data = pd.read_csv("./realistic_tide.csv",sep=",",index_col=0,header=0)

scha24_path = "./scha24_6spec"
scha23_path = "./scha23"
vmd_path = './vmd'
morlet_path ='./morelet'

species_list = [f"M1_{item}_parm","M2_{item}_parm","M3_{item}_parm","M4_{item}_parm"]
species_labels =  ["D1","D2","D3","D4"]
ylims = {"D1_amp":(0,0.75),"D2_amp":(0,1.5),"D3_amp":(0,0.5),"D4_amp":(0,0.4)}

species_integers = [1,2,3,4]



def process_phase(x,centradhr):
    return x
    t = x.index.to_numpy()
    z = x.copy()
    z[:] = np.unwrap(x.to_numpy()) - centradhr*t
    return z



def do():

    scha24 = pd.read_csv(os.path.join(scha24_path,"realistic_scha24.csv"),sep=",",header=0,index_col=0)
    scha23 = pd.read_csv(os.path.join(scha23_path,"realistic_scha23.csv"),sep=",",header=0,index_col=0)
    vmd = pd.read_csv(os.path.join(vmd_path,"realistic_vmd.csv"),sep=",",header=0,na_values="NA",skipinitialspace=True,index_col=0)
    analytical = pd.read_csv("realistic_tide.csv",sep=",",header=0,skipinitialspace=True,index_col=0,na_values="NA")
    B=4
    C=0.966
    wspec= f'cmor{B:.1f}-{C}'
    cwt_outfile = os.path.join(morlet_path,wspec+"_result.csv")
    mor = pd.read_csv(cwt_outfile,sep=",",index_col=0) 



    vmdcol='0.7'
    vmdline = '-.'
    anycol='maroon'

    

    nrows=5
    ncols=1
    fig = plt.figure(figsize=(6,7.5))
    gs = gridspec.GridSpec(nrows, ncols,
         left=0.1,right=.88,
         bottom=None,top=None,
         wspace=0.25,hspace=0.35)


    ax = fig.add_subplot(gs[0,0])
    ax.set_prop_cycle(cycler('color',sns.color_palette("deep",n_colors=7)))
    axes = [ax]       
       
    for i,(r,c) in enumerate([(r,c) for r in range(1,nrows) for c in range(ncols)]):
        axnew = fig.add_subplot(gs[r,c],sharex=axes[0])
        
        axnew.set_prop_cycle(cycler('color',sns.color_palette("deep",n_colors=7)))
        axes.append(axnew)
        
        


    ax0 = axes[0]
    anysub = analytical['sub']
    anysub.index = anysub.index/24.
    anysub.plot(ax=ax0,color=anycol)
    sc23sub = scha23["subtide_parm"]
    sc23sub.index = sc23sub.index/24.
    sc23sub.plot(ax=ax0,color="0.2",style=[':'])
    sc24sub = scha24["subtide_parm"]
    sc24sub.index = sc24sub.index/24.
    sc24sub.plot(ax=ax0,color="0.4")
    vmdsub = vmd['subtide']
    vmdsub.index = vmdsub.index/24.
    vmdsub.plot(ax=ax0,color=vmdcol,linestyle=vmdline)
    ax0.set_title("Sub",x=0.88,y=1,pad=-12)

    


    for ispecies in range(len(axes)):
        if (i < 5): plt.setp(axes[i].get_xticklabels(), visible=False)




    for ispecies in range(1,5):
        ax = axes[ispecies]
        lab=f"D{ispecies}_{item}"
        ax.set_title(f"D{ispecies}",x=0.88,y=1.,pad=-12)
        analy=analytical[lab]
        print(ispecies,cfreq[5-ispecies])

        col = f"M{ispecies}_{item}_parm"
        sc23 = scha23[col]
        #sc23 = process_phase(-sc23,cfreq[ispecies])
        sc24 = scha24[col]
        #sc24 = process_phase(-sc24,cfreq[ispecies])
        vm = vmd[lab]
        #vm = process_phase(vm,cfreq[ispecies])
        mo = mor[lab]
        #mo = process_phase(mo,cfreq[ispecies])

        for ndxitem in [sc24,sc23,vm,analy,mo]:
            ndxitem.index = analytical.index/24.
            ndxitem.index.name = "Time (days)"
        
        analy.plot(ax=ax,label="Analytical",color=anycol)        
        sc23.plot(ax=ax, label="SCHA(2,3)",color='0.2',style=[':'])        
        sc24.plot(ax=ax, label="SCHA(2,4)",color='0.4',style=['-'])
        if ispecies > 1: vm.plot(ax=ax,label="VMD",color=vmdcol,linestyle=vmdline,zorder=0)
        mo.plot(ax=ax,label="Morlet, B=4",color='0.4',style=['--'])
        if ispecies == 1:
            subax=ax.inset_axes([0.38,0.82,0.35,0.65])
            analy.plot(ax=subax,color=anycol)
            sc23.plot(ax=subax,color='0.2',style=[':'])
            sc24.plot(ax=subax,color="0.1")
            vm.plot(ax=subax,color=vmdcol,linestyle=vmdline,zorder=0)
            mo.plot(ax=subax,color="0.4",style=['--'])
    
            subax.set_xlim(71,85)
            subax.minorticks_off()
            subax.set_ylim(0.0,0.6)
            subax.get_xaxis().set_visible(False)
            subax.get_yaxis().set_visible(False)
            ax.indicate_inset_zoom(subax,edgecolor="black")

        
        if ispecies == 4:
            subax=ax.inset_axes([0.42,0.65,0.35,0.6])
            analy.plot(ax=subax,color=anycol)
            sc23.plot(ax=subax,color='0.2',style=[':'])
            sc24.plot(ax=subax,color="0.1")
            vm.plot(ax=subax,color=vmdcol,linestyle=vmdline,zorder=0)
            mo.plot(ax=subax,color="0.4",style=['--'])
    
            subax.set_xlim(84,95)
            subax.minorticks_off()
            subax.set_ylim(0.06,0.2)
            subax.get_xaxis().set_visible(False)
            subax.get_yaxis().set_visible(False)            
            ax.indicate_inset_zoom(subax,edgecolor="black")
        
     
        if lab in ylims: ax.set_ylim(*ylims[lab])
    axes[2].legend(bbox_to_anchor=[0.0,3.],loc="lower left")
    axes[0].set_xlim(0,3000/24.)
    axes[2].set_ylabel("Amplitude")
    axes[4].set_xlabel("Time (days)")
    


    
    plt.tight_layout()
    plt.savefig(f'realistic_amp.png')
    plt.savefig(f'realistic_amp.pdf')    
    
    plt.show()
        
if __name__ == "__main__":
    do()        
        