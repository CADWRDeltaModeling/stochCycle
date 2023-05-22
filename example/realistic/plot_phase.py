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
item = "phase"
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

ylims = {"D1_phase":(0,180),"D2_phase":(-55,50.),"D3_phase":(0,200),"D4_phase":(-120,80)}

species_integers = [1,2,3,4]


def process_phase(x,centradhr):
    print(centradhr)
    
    t = x.index.to_numpy()
    print(t[0:5])
    print(t[-1])
    z = x.copy()
    z[:] = np.unwrap(x.to_numpy()) - centradhr*t
    return z



def do():

    scha24 = pd.read_csv(os.path.join(scha24_path,"realistic_scha24.csv"),sep=",",header=0,index_col=0)
    scha23 = pd.read_csv(os.path.join(scha23_path,"realistic_scha23.csv"),sep=",",header=0,index_col=0)
    vmd = pd.read_csv(os.path.join(vmd_path,"realistic_vmd.csv"),sep=",",header=0,na_values="NA",skipinitialspace=True,index_col=0)
    analytical = pd.read_csv("realistic_tide.csv",sep=",",header=0,skipinitialspace=True,index_col=0,na_values="NA")
    
    vmdcol='0.7'
    vmdline = '-.'
    anycol='maroon'
    col23 = "0.2"
    col24 = "0.1"
    colmo = "0.4"
    t=analytical.index.to_numpy()

    

    nrows=5 if item == "amp" else 4
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
        
        

    if item == "amp":
        ax0 = axes[0]
        scha23["subtide_parm"].plot(ax=ax0,color="0.2")
        scha24["subtide_parm"].plot(ax=ax0,color="0.1")
        vmd["subtide"].plot(ax=ax0,color=vmdcol,linestyle=vmdline)
        ax0.legend(["SCHA(2,3)","SCHA(2,4)","VMD"],bbox_to_anchor=(0.7, 1.1),ncol=3)
        offset_spec = 0
    else: 
        offset_spec = 1
    
    #subax = ax0.inset_axes([0.6,0.6,0.3,0.3])
    #scha23["subtide_parm"].plot(ax=subax,color="0.5")
    #scha24["subtide_parm"].plot(ax=subax,color="0.1")
    #vmd["subtide"].plot(ax=subax,color=vmdcol)
    
    #subax.set_xlim(pd.Timestamp(2012,6,2),pd.Timestamp(2012,6,3))
    #subax.minorticks_off()
    #subax.set_ylim(2.6,2.9)
    #ax0.indicate_inset_zoom(subax,edgecolor="black")

    for ispecies in range(len(axes)):
        if (i < ncols): plt.setp(axes[i].get_xticklabels(), visible=False)

    B=4.0
    C=0.966
    wspec= f'cmor{B}-{C}'
    cwt_outfile = os.path.join(morlet_path,wspec+"_result.csv")
    mor = pd.read_csv(cwt_outfile,sep=",",index_col=0)
    adjphi = {1: 0.0096, 2:0.0, 3:0.0, 4: 0.0}
    adjphi = {1: 0.0, 2:0.0, 3:0.0, 4: 0.0}

    for ispecies in range(1,5):
        ax = axes[ispecies-offset_spec]
        lab=f"D{ispecies}_{item}"
        ax.set_title(f"D{ispecies}")
        #analy=process_phase(analytical[lab],cfreq[ispecies-1])
        analy=process_phase(analytical[lab],0.)
        
        if ispecies == 3:
            analy = analy - 6.2*2.*np.pi*t/4000 - 2*np.pi
            analy = analy.where(analy > (-6.5),analy+2*np.pi) 
            analy = analy + 2*np.pi
            
        analy = analy.where(analy>(2.*np.pi),analy - (2.*np.pi))
        if (ispecies ==1): 
            analy = analy - 6.2*2.*np.pi*t/4000
            analy = analy.where(analy > (-2.*np.pi),analy+2*np.pi)
            analy = analy.where(analy > (-2.*np.pi),analy+2*np.pi)
        if (ispecies ==2): analy = analy.where(analy > (-4),analy+2*np.pi)       
        if (ispecies ==4): analy = analy.where(analy > (-4),analy+2*np.pi)          
        print(ispecies,cfreq[5-ispecies])

        col = f"M{ispecies}_{item}_parm"
        sc23 = -scha23[col]
        sc23 = process_phase(sc23,cfreq[ispecies-1]+adjphi[ispecies])
        sc24 = -scha24[col]
        sc24 = process_phase(sc24,cfreq[ispecies-1]+adjphi[ispecies])
        if (ispecies == 3): 
            sc24 = sc24.where(sc24 < 7.,sc24-2*np.pi)   
        if (ispecies == 3): 
            sc23 = sc23 - 6.18*2.*np.pi*t/4000
            sc24 = sc24 - 6.18*2.*np.pi*t/4000            
            sc23 = sc23.where(sc23 < 3.75, sc23-2*np.pi)            
            sc24 = sc24.where(sc24-sc23 > -1., sc24+2*np.pi)
            #sc24 = sc24.where(sc24 < 6.,sc24-2*np.pi)
            #sc24 = sc24.where(sc24 < 6.,sc24-2*np.pi)    
            #sc24 = sc24.where(sc24 > -2.,sc24+2*np.pi)
            #sc24 = sc24.where(sc24 < 6,sc24-2*np.pi)
        
        if ispecies == 1: 
            sc23 = sc23 - 6.2*2.*np.pi*t/4000
            sc24 = sc24 - 6.2*2.*np.pi*t/4000            
            sc24 = sc24.where(sc23 - sc24 < 4,sc24+2*np.pi)
            sc24 = sc24.where(sc24 <6,sc24-2*np.pi)
            sc24 = sc24.where(sc24 <6,sc24-2*np.pi)
        if ispecies == 4: 
            sc24 = sc24.where(sc24 > -4,sc24+2*np.pi)
            sc23 = sc23.where(sc23<2.5,sc23-2*np.pi)
            
        
        vm = vmd[lab]
        vm = process_phase(vm,cfreq[ispecies-1])
        if (ispecies ==1):vm = vm - 6.2*2.*np.pi*t/4000
        if (ispecies ==3):vm = vm - 6.18*2.*np.pi*t/4000
        if (ispecies ==3): vm = vm.where(vm > -40.,vm+12*np.pi)   
        if (ispecies ==3): vm = vm.where(vm > -40.,vm+12*np.pi)         
        if (ispecies ==3): vm = vm.where(vm > -7.,vm+2*np.pi) 
        if (ispecies ==3): vm = vm.where(vm > -7.,vm+2*np.pi)  
        if (ispecies ==4): vm = vm.where(vm <2.5,vm-2*np.pi) 


        mo = mor[lab]

        print("adjphi",adjphi)
        #mo = process_phase(mo,cfreq[ispecies-1]+adjphi[ispecies])
        #mo = process_phase(mo,cfreq[ispecies-1])
        mo = process_phase(mo,cfreq[ispecies-1]+adjphi[ispecies])      
        if ispecies == 1: 
            mo = mo - 6.2*2.*np.pi*t/4000
        if ispecies == 3: 
            mo = mo + 0.05*2.*np.pi*t/4000
        
        if (ispecies == 3): mo = mo - 12.5*np.pi*t/4000        
        if (ispecies ==3): mo = mo.where(mo < 0.,mo-4*np.pi) 
        if (ispecies ==3): mo = mo.where(mo < 0.,mo-2*np.pi)         
        if (ispecies ==3): mo = mo.where(mo > -4.,mo+2*np.pi)  
        if (ispecies ==3): mo = mo.where(mo > -4.,mo+2*np.pi)  
        if (ispecies ==3): mo = mo.where(mo < 1 ,mo-2*np.pi)

        datum_adj = {1: 4., 2:0, 3:4.,4:0}

        for ndxitem in [sc24,sc23,vm,analy,mo]:
            print("Adjusting",ispecies,datum_adj[ispecies])
            ndxitem += datum_adj[ispecies]
            ndxitem.index = analytical.index/24.
            ndxitem.index.name = "Time (days)"
            ndxitem.iloc[:] = ndxitem*360/(2.*np.pi)
            print(ndxitem)
        

        #else: anyadj = analy
        analy.plot(ax=ax,color=anycol,label="Analytical",zorder=100)
        sc23.plot(ax=ax, label="SCHA(2,3)",color=col23,style=[':'])        
        sc24.plot(ax=ax, label="SCHA(2,4)",color=col24,style=['-'])
        if ispecies != 12: vm.plot(ax=ax,label="VMD",color=vmdcol,linestyle=vmdline)
        mo.plot(ax=ax,label="Morlet, B=4",color=colmo,style=['--'])
     
        print("lab",lab)
        if lab in ylims: ax.set_ylim(ylims[lab])
        
        
        if ispecies == 7:
            subax=ax.inset_axes([0.73,0.79,0.35,0.65])
            anyadj.plot(ax=subax,color=anycol)
            sc23.plot(ax=subax,color='0.2',style=[':'])
            sc24.plot(ax=subax,color="0.1")
            vm.plot(ax=subax,color=vmdcol,linestyle=vmdline,zorder=0)
            mo.plot(ax=subax,color="0.4",style=['--'])
    
            subax.set_xlim(71,85)
            subax.minorticks_off()
            subax.set_ylim(-4.0,4.0)
            subax.get_xaxis().set_visible(False)
            subax.get_yaxis().set_visible(False)
            ax.indicate_inset_zoom(subax,edgecolor="black")
        
        
    #axes[1].legend(bbox_to_anchor=(0.7,2),loc="lower left")
    axes[1].set_ylabel("Phase (degrees)")
    axes[1].set_xlim(0,3000./24)
    plt.tight_layout()
    
    plt.savefig(f'realistic_phase.png')
    plt.savefig(f'realistic_phase.pdf')       
    plt.show()
        
if __name__ == "__main__":
    do()        
        