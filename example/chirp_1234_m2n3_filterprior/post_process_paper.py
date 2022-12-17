#!/usr/bin/env python
# -*- coding: utf-8 -*-
import pandas as pd
import matplotlib.pyplot as plt
from vtools.datastore.read_ts import *
import numpy as np

"""
Index(['subtide', 'subtide_parm', 'subtide_postmean', 'subtide25', 'subtide50',
       'subtide75', 'D1_phase', 'D1_phase_parm', 'D1_phase_postmean',
       'D2_phase', 'D2_phase_parm', 'D2_phase_postmean', 'D3_phase',
       'D3_phase_parm', 'D3_phase_postmean', 'D4_phase', 'D4_phase_parm',
       'D4_phase_postmean', 'D1_amp', 'D1_amp_parm', 'D1_amp_postmean',
       'D1_amp_25', 'D1_amp_50', 'D1_amp_75', 'D2_amp', 'D2_amp_parm',
       'D2_amp_postmean', 'D2_amp_25', 'D2_amp_50', 'D2_amp_75', 'D3_amp',
       'D3_amp_parm', 'D3_amp_postmean', 'D3_amp_25', 'D3_amp_50', 'D3_amp_75',
       'D4_amp', 'D4_amp_parm', 'D4_amp_postmean', 'D4_amp_25', 'D4_amp_50',
       'D4_amp_75'],
      dtype='object')
"""
m=2
n=3
markerspace = 24
plt.style.use(['seaborn-paper','seaborn-colorblind'])
output=pd.read_csv(f"chirp{m}{n}.csv",index_col='times')
output.index.name="time"
print(output.columns)

fig,(ax0,ax1,ax2,ax3)=plt.subplots(4,1,sharex=True,figsize=(7,6))

output.subtide.plot(ax=ax0,label='Analytical',color='black')
output.subtide_parm.plot(ax=ax0,label='SCHA(2,3) Mean Param',color='0.4',markersize=5,markevery=markerspace,style=['-^'])
#ax0.fill_between(output.index,output.subtide25,output.subtide75,color='.8')
#output.subtide25(ax=ax0,label='Posterior 25 quantile')
#output.subtide50(ax=ax0,label='Posterior 50 quantile')
#output.subtide50(ax=ax0,label='Posterior 75 quantile')

output.D1_amp.plot(ax=ax1,label='Analytical',color='black')
output.D1_amp_parm.plot(ax=ax1,label='SCHA(2,3) Mean Param',
                        color='0.4',markersize=5,markevery=markerspace,style=['-^'])
ax1.fill_between(output.index,output.D1_amp_25,output.D1_amp_75,color='.8')

output.D2_amp.plot(ax=ax2,label='Analytical',color='black')
output.D2_amp_parm.plot(ax=ax2,label='SCHA(2,3) Mean Param',color='0.4',markersize=5,markevery=markerspace,style=['-^'])
ax2.fill_between(output.index,output.D2_amp_25,output.D2_amp_75,color='.8')

output.D3_amp.plot(ax=ax3,label='Analytical',color='black')
output.D3_amp_parm.plot(ax=ax3,label='SCHA(2,3) Mean Param',color='0.4',markersize=5,markevery=markerspace,style=['-^'])
ax3.fill_between(output.index,output.D3_amp_25,output.D3_amp_75,color='.8')

#[ax.set_ylabel("amplitude") for ax in (ax0,ax1,ax2,ax3)]
#ax0.legend(loc='lower right')
ax0.legend(loc='lower right',bbox_to_anchor=(0.5, 1.05, 0.5, 0.5),framealpha=0.9)

ax0.set_xlim(0.,960.)
ax0.annotate("Subtide",(0.88,0.8),xycoords='axes fraction')
ax1.annotate("D1",(0.9,0.8),xycoords='axes fraction')
ax2.annotate("D2",(0.9,0.8),xycoords='axes fraction')
ax3.annotate("D3",(0.9,0.8),xycoords='axes fraction')

ax1.set_ylabel("Amplitude")
#plt.tight_layout()
plt.show()


fig,(ax1,ax2,ax3)=plt.subplots(3,1,sharex=True,figsize=(7,5))
(180.*output.D1_phase/np.pi).plot(ax=ax1,label='Analytical',color='black')
(180.*output.D1_phase_parm/np.pi).plot(ax=ax1,label=f'SCHA({m},{n}) Mean Param',color='0.4',markersize=5,markevery=markerspace,style=['-^'])

(180.*output.D2_phase/np.pi).plot(ax=ax2,label='Analytical',color='black')
(180.*output.D2_phase_parm/np.pi).plot(ax=ax2,label=f'SCHA({m},{n}) Mean Param',color='0.4',markersize=5,markevery=markerspace,style=['-^'])

(180.*output.D3_phase/np.pi).plot(ax=ax3,label=f'Analytical',color='black')
(180.*output.D3_phase_parm/np.pi).plot(ax=ax3,label=f'SCHA({m},{n}) Mean Param',color='0.4',markersize=5,markevery=markerspace,style=['-^'])

ax1.legend(loc='lower right',bbox_to_anchor=(0.5, 1.05, 0.5, 0.5),framealpha=0.9)
ax1.set_xlim(0.,960.)
ax1.annotate("D1",(0.9,0.8),xycoords='axes fraction')
ax2.annotate("D2",(0.9,0.8),xycoords='axes fraction')
ax3.annotate("D3",(0.9,0.8),xycoords='axes fraction')

ax2.set_ylabel("Phase (degrees)")
#plt.tight_layout()
plt.show()







