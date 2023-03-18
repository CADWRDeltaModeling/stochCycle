#!/usr/bin/env python


"""
This script produces Figure 2

There is a dummy plot that may be useful for looking at filter gain if you want to check out some parameters.
"""

import matplotlib.gridspec as gridspec
from scipy.special import comb
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
from scipy.special import comb
import pandas as pd
import matplotlib.patches as patches
def freqresponse(lambdaf,m,q):
    elambdaf = -1.j*lambdaf
    elambdaf = np.exp(elambdaf)
    denom = 1. + np.power((1.-elambdaf)*(1.-1./elambdaf),m)/q
    fr = 1./denom
    gain = np.absolute(fr)
    return gain

def lowpass_gain_term(lambdaf,m,q):
    elambdaf = -1.j*lambdaf
    elambdaf = np.exp(elambdaf)
    return q/np.power((1.-elambdaf)*(1.-1./elambdaf),m)


def lowpass_gain(lambdaf,m,q):
    print("q=",q)
    out = np.power(2.-2.*np.cos(lambdaf),float(-m))*q
    return out

def lowpass_gain2(lambdaf,m,q):
    print("q=",q)
    out = 1.+np.power(2.-2.*np.cos(lambdaf),float(m))/q
    return 1./out


def spectral_density2(lambdaf,lambdac,rho,nn,q):
    assert q is not None
    rho2 = rho*rho
    clc = np.cos(lambdac)
    num = 1. + rho2*clc*clc - 2.*rho*clc*np.cos(lambdaf)
    denom = 1. + rho2*rho2 + 4.*rho2*clc*clc \
           - 4.*(rho + rho*rho2)*clc*np.cos(lambdaf) + 2.*rho2*np.cos(2.*lambdaf)
    return num/denom

def butter_c(L,rho,lambdac):
    coslambdacL=np.cos(lambdac*L)
    return (1-rho*coslambdacL)/(1.-2.*rho*coslambdacL+rho*rho*L*L)

def spectral_density_butter(lambdaf,lambdac,rho,nn,q=None):
    dlambdaf = np.exp(-1.j*lambdaf)
    return q*np.power(butter_c(dlambdaf,rho,lambdac)*butter_c(1./dlambdaf,rho,lambdac),nn)

def acgf_aj(kk,L,rho,nn):
    return np.power(L,nn-1)*comb(nn,kk)*np.power(-rho,kk)*np.power(L,kk)
    
def acgf_balanced(lambdaf,lambdac,rho,nn,q):
    #L = np.cos(-lambdaf) + 1j*np.sin(-lambdaf)
    L = np.exp(-1.j*lambdaf)
    #Lm1=np.cos(lambdaf) + 1j*np.sin(lambdaf)
    Lm1 = np.exp(1.j*lambdaf)
    L2 = L*L
    rho2 = rho*rho
    sum0=0.;sum1=0.;sum2=0.;sum3=0.
    for kk in np.arange(nn+1):
        kkf = float(kk)
        sum0 += np.cos(kkf*lambdac)*acgf_aj(kkf,L,rho,nn)
        sum1 += np.cos(kkf*lambdac)*acgf_aj(kkf,Lm1,rho,nn)
        sum2 += np.sin(kkf*lambdac)*acgf_aj(kkf,L,rho,nn)
        sum3 += np.sin(kkf*lambdac)*acgf_aj(kkf,Lm1,rho,nn)     
        print("run sum ",sum0*sum1,sum2*sum3)
        
    ttot0 = 0.
    ttot1 = 0.
    for jjj in np.arange(nn+1.):
        for kkk in np.arange(nn+1.):
            ttot0 += comb(nn,jjj)*comb(nn,kkk)*np.cos(float(jjj-kkk)*lambdaf)*np.cos(float(jjj)*lambdac)*np.cos(float(kkk)*lambdac)*np.power(-rho,jjj+kkk)
            ttot1 += comb(nn,jjj)*comb(nn,kkk)*np.cos(float(jjj-kkk)*lambdaf)*np.sin(float(jjj)*lambdac)*np.sin(float(kkk)*lambdac)*np.power(-rho,jjj+kkk)
            print("running ",jjj,kkk,ttot0,ttot1)
    print("here",lambdac,lambdaf)        
    print(sum0*sum1,sum2*sum3) 
    print("other way")    
    print(ttot0,ttot1)
    print(ttot0/(sum0*sum1),ttot1/(sum0*sum1))
    num = sum0*sum1 + sum2*sum3
    denom = np.power(1. - 2.*rho*np.cos(lambdac)*L + rho2*L2,float(nn))*np.power(1.-2.*rho*np.cos(lambdac)*Lm1 + rho2*Lm1*Lm1,float(nn))
    other = np.power(1.+4.*rho2*np.cos(lambdac)**2. + rho2*rho2 - 4.*rho*(1.+rho2)*np.cos(lambdac)*np.cos(lambdaf) + 2.*rho2*np.cos(2.*lambdaf),float(nn))
    print("denoms",denom)
    print("other",other)
    return (q/(2.*np.pi))*num/denom

    

def spectral_density_real(lambdaf,lambdac,rho,nn,q=None):
    n = float(nn)
    rho2 = rho*rho
    num = np.zeros_like(lambdaf)
    for j in range(nn+1):
        for k in range(nn+1):
            jkdiff = float(j-k)
            jksum = float(j+k)
            term = np.power(-1.,jksum)*comb(nn,j)*comb(nn,k)*np.cos(lambdac*jkdiff)*np.cos(lambdaf*jkdiff)*np.power(rho,jksum)
            num += term 
            #print("Term {} {} {}".format(j,k,term))
    denom= 1.+ 4.*rho2*np.cos(lambdac)**2. + rho2*rho2 - 4.*rho*(1.+rho2)*np.cos(lambdac)*np.cos(lambdaf)+2.*rho2*np.cos(2.*lambdaf)
    denom = np.power(denom,n)
    
    two_pi = np.pi*2.
    if q is not None:
        # return power spectrum (Trimbur 2006 pg 11, no eq number)
        g = (q/two_pi)*(num/denom)
        return g


    g_without_sigma_term = num/(two_pi*denom)

    #return spectral density
    # these are the numerator and denominator of the variance (16) ignoring 
    # q or sigma_kappa because it divides out    
    var_den = (1.-rho2)**(2.*nn-1.)
    var_num = 0.
    for i in range(nn):
        rho_to_2i = rho2**float(i)
        var_num += np.square(comb(nn-1,i)*rho_to_2i)
    
    var_without_sigma_term = var_num/var_den 
    return g_without_sigma_term/var_without_sigma_term



def plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax,linestyle,samples_tideday):
    if len(species) != len(sigma2kappa):
        raise ValueError("Length of species must match length of sigm2kappa")
    

    lambdac = 2.*np.pi/float(samples_tideday)

    #Frequencies in response plot
    nsample = 5000
    lambdaf = np.linspace(1./nsample,18*np.pi/float(samples_tideday),nsample)
    lambdaflab = np.linspace(0,9,nsample)


    qtrend= sigma2zeta/sigma2eps
    qkappa=sigma2kappa/sigma2eps

    print("freqs",lambdac,lambdaf)
    calc_meth = acgf_balanced #spectral_density_butter  #spectral_density_real #
    calc_meth = spectral_density_real
    sdense0 = lowpass_gain_term(lambdaf,m,q=qtrend)  *1  # subtidal response
    
    gains = [sdense0]
    totaldense = sdense0
    for i,(a,b) in enumerate(zip(qkappa,species)):
        print(a,b)
    
    ax.set_prop_cycle(None)
    for ispec,(qkap,spec) in enumerate(zip(qkappa,species)):
        gain = calc_meth(lambdaf,spec*lambdac,rho,n,q=qkap).flatten()
        gains.append(gain)    #diurnal
        totaldense = totaldense+gain
    print(totaldense.shape)
    
    one = 1. 
    sdenses = []
    labels = ["sub","D1","D2","D3","D4","D5","D6"]
    ax.set_prop_cycle(color=['black']+sns.color_palette("deep",n_colors=7))
    for lab,gain in zip(labels,gains):
        sdense=gain/(totaldense+one)
        sdenses.append(sdense)
        ax.plot(lambdaflab,sdense,label=lab,linestyle=linestyle)


plt.style.use(['seaborn-paper'])
plt.margins(x=0,y=0)
samples_tideday=99

##### This not used in paper, but it is a good dummy for plotting a new set of parameters and also 
# the figure is used to provide line styles that are used for Figure 2. 
fig,dummyax = plt.subplots(1)
linestyle = "-"
m = 2
n = 4
rho = .994
sigma2eps = 0.0022
sigma2zeta = 2.0e-7
species = np.array([1.,2.,3.,4.,5.,6.])
sigma2kappa = np.array([1.35e-13,1.35e-13,1.35e-13,1.35e-13,1.35e-13,1.35e-13])
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,dummyax,linestyle,samples_tideday)
dummyax.set_title("Ignore This Plot. Only used for line styles")

#######################

nrow=3
ncol=2
gs = gridspec.GridSpec(3, 2,
         wspace=0.48, hspace=0.74, 
         top=1.-0.5/(nrow+1), bottom=0.5/(nrow+1), 
         left=0.3/(ncol+1), right=1-0.3/(ncol+1))

#fig,ax=plt.subplots(3,2,figsize=(6,7)) #,layout='constrained')
fig = plt.figure(figsize=(6.5,7))
ax = []
ax.append([plt.subplot(gs[0,0]),plt.subplot(gs[0,1])])
ax.append([plt.subplot(gs[1,0]),plt.subplot(gs[1,1])])
ax.append([plt.subplot(gs[2,0]),plt.subplot(gs[2,1])])
ax = np.array(ax)

# Top left plot: (2,3) vs (2,4)
rho = .994 
sigma2eps = 0.0022
sigma2zeta = 2.0e-6 
species = np.array([1.,2.,3.,4.])
sigma2kappa = np.array([1.35e-13,1.35e-13,1.35e-13,1.35e-13])
m = 2
n = 4
linestyle = "-"
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[0,0],linestyle,samples_tideday)
# SCHA (2,3)
linestyle = (5,(10,3))
n = 3
plot_freq(m,n,0.996,sigma2eps,sigma2zeta/2,species,sigma2kappa,ax[0,0],linestyle,samples_tideday)

#top right plot: rho=0.984 vs 0.994
linestyle = "-"
m = 2
n = 4
rho = .994 
sigma2eps = 0.0022
sigma2zeta = 2.0e-6 
species = np.array([1.,2.,3.,4.])
sigma2kappa = np.array([1.35e-13,1.35e-13,1.35e-13,1.35e-13])
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[0,1],linestyle,samples_tideday)
lines=dummyax.get_lines()
print("# lines", len(lines))
freqlabels = ["Subtidal","D1","D2","D3","D4","D5","D6"]
linestyles=["-",(5,(10,3)),"-."]

# generate a black legend for line styles

dummy_lines = [] 
for b_idx, b in enumerate(linestyles):
    dummy_lines.append(ax[0,1].plot([],[], c="black", ls = linestyles[b_idx])[0])



linestyle = (5,(10,3))
rho = .984
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[0,1],linestyle,samples_tideday)

#middle left plot: sigma2zeta = 2e-7,2e-6,2e-5
linestyle = "-"
m = 2
n = 4
rho = .994
sigma2eps = 0.0022
sigma2zeta = 2.0e-7
species = np.array([1.,2.,3.,4.])
sigma2kappa = np.array([1.35e-13,1.35e-13,1.35e-13,1.35e-13])
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[1,0],linestyle,samples_tideday)

sigma2zeta = 2.0e-6
linestyle = (5,(10,3))
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[1,0],linestyle,samples_tideday)
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[1,0],linestyle,samples_tideday)

sigma2zeta = 2.0e-5
linestyle = "-."
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[1,0],linestyle,samples_tideday)

#middle right plot: Nspecies 4, 6
linestyle = "-"
m = 2
n = 4
rho = .994
sigma2eps = 0.0022
sigma2zeta = 2.0e-6
species = np.array([1.,2.,3.,4.])
sigma2kappa = np.array([1.35e-13,1.35e-13,1.35e-13,1.35e-13])
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[1,1],linestyle,samples_tideday)

linestyle = (5,(10,3))
species = np.array([1.,2.,3.,4.,5.,6.])
sigma2kappa = np.array([1.35e-13,1.35e-13,1.35e-13,1.35e-13,1.35e-13,1.35e-13])
plot_freq(m,n,rho,sigma2eps,sigma2zeta,species,sigma2kappa,ax[1,1],linestyle,samples_tideday)

modlabels=["SCHA(2,4), 4 species","SCHA(2,3), 4 species"]



leg00=ax[0,0].legend([dummy_lines[i] for i in range(len(linestyles)-1)],
               ['SCHA(2,4)','SCHA(2,3)'], 
               bbox_to_anchor=(0.62,0.57),
               loc='lower left',
               title="n", 
               ncol=1)
leg00.get_frame().set_alpha(0.92)

fig.legend(lines, freqlabels, bbox_to_anchor=(1.2,.6),title="Species")
leg11=ax[0,1].legend([dummy_lines[i] for i in range(len(linestyles)-1)],
               ['0.994','0.984'], 
               loc='lower left', 
               bbox_to_anchor=(0.68,.59),
               title=r'$\rho$')
leg00.get_frame().set_alpha(0.92)

legend101 = ax[1,0].legend([lines[i] for i in range(5)], freqlabels,
                           bbox_to_anchor=(1.2,.6),loc='lower left',title="Species")
leg102=ax[1,0].legend([dummy_lines[i] for i in range(len(linestyles))],
                      ['2.0e-7','2.0e-6','2.0e-5'],
                       bbox_to_anchor=(0.65,.37),
                       loc='lower left',
                       title=r'$\sigma ^2 _{\zeta}$')
leg102.get_frame().set_alpha(0.92)



leg11=ax[1,1].legend([dummy_lines[i] for i in range(len(linestyles)-1)], 
                     ['4','6'], 
                     bbox_to_anchor=(0.69,0.58),
                     loc='lower left', 
                     title="$N_c$")
leg11.get_frame().set_alpha(0.92)
ax[0,0].set_xlim(0,7)
ax[0,0].set_xticks([0,1,2,3,4,5,6,7])
ax[0,0].set_ylabel("G($\lambda$)")
ax[0,0].set_xlabel("Frequency (cycles/tidal day)")

ax[0,1].set_xlim(0,7)
ax[0,1].set_xticks([0,1,2,3,4,5,6,7])
ax[0,1].set_ylabel("G($\lambda$)")
ax[0,1].set_xlabel("Frequency (cycles/tidal day)")

ax[1,0].set_xlim(0,7)
ax[1,0].set_xticks([0,1,2,3,4,5,6,7])
ax[1,0].set_ylabel("G($\lambda$)")
ax[1,0].set_xlabel("Frequency (cycles/tidal day)")

ax[1,1].set_xticks([0,1,2,3,4,5,6,7,8])
ax[1,1].set_ylabel("G($\lambda$)")
ax[1,1].set_xlabel("Frequency (cycles/tidal day)")

#bottom right: color (mse error) & contour (filter quality)
from scipy.interpolate import interp2d
data=pd.read_csv("out_scha_cv.csv",skipinitialspace=True,comment="#",)
#data = data.apply(pd.to_numeric, errors='coerce')
#data = data.dropna()
data = data[data.rho < 0.996]
x = data["rho"]
y = data["sigma2_zeta"]
z = data["cv_mse"]
f = interp2d(x, y, z, kind='linear')
x2 = np.linspace(0.991, 0.996, 6)
y2 = np.power(10.,np.linspace(-7,-5,6))
z2 = f(x2,y2)
print(data[['rho','sigma2_zeta','cv_mse']])
X2, Y2 = np.meshgrid(x2, y2)
q = data["filter_quality"]/10000.
fq = interp2d(x, y, q, kind='linear')
q2 = fq(x2,y2)


d3 = data.sort_values(by=['rho','sigma2_zeta'])
xvals = d3['rho'].unique()
yvals = d3['sigma2_zeta'].unique()
zvals = d3['cv_mse'].values.reshape(len(xvals), len(yvals)).T
cm=ax[2,1].pcolormesh(xvals,yvals,zvals)
ax[2,1].set_yscale('log')
fig.colorbar(cm,label='rmse')
qvals = d3['filter_quality'].values.reshape(len(xvals), len(yvals)).T/1.6e5
cp = ax[2,1].contour(xvals, yvals, qvals,colors="0.9",vmin=180,vmax=280) #,levels = [255,258,260,265])

ax[2,1].clabel(cp,inline=1, fontsize=10)
ax[2,1].set_ylabel(r'$\sigma ^2 _{\zeta}$')
ax[2,1].set_xlabel(r'$\rho$')

# get data you will need to create a "background patch" to your plot
xmin, xmax = ax[2,1].get_xlim()
ymin, ymax = ax[2,1].get_ylim()
xy = (xmin,ymin)
width = xmax - xmin
height = ymax - ymin

# create the patch and place it in the back of countourf (zorder!)
p = patches.Rectangle(xy, width, height, hatch='///',color="0.7",fill=None, zorder=-10)
ax[2,1].add_patch(p)


# Plot 2(e), bottom left: Spectral density sample frequencies and target values
from itertools import repeat
import math
def plot_filter_gain_prior(species,x1,x0,x_pass_stop,ax,color):
    y1 = np.repeat(0.2,len(x1))
    y0 = np.repeat(0,len(x0))
    ax.plot(x1, y1+species, 'X', color=color,markersize=5)
    ax.plot(x0, y0+species, 'o', color=color,markersize=2.0)


    
#Subtidal
x1=[0.10,0.20,0.30,0.40]
x0=[0.94,0.96,0.98,1.00,1.04,1.08,1.45,1.50,1.55,1.85,1.95,2.00,2.5,3.0,3.5,4.0,
   4.50,5.0,5.50,6.0,6.5,7.0,7.5,8.0]
plot_filter_gain_prior(0,x1,x0,[(0.4,0.94)],ax[2,0],'black')
#D1
x1=[0.96,0.98,1.00,1.02,1.04,1.06,1.08]
x0=[0.35,0.40,0.45,1.75,1.85,2.00,3.00,4.00,3,5,6]
plot_filter_gain_prior(1,x1,x0,[(0.96,0.45),(1.08,1.75)],ax[2,0],sns.color_palette("deep",n_colors=7)[0])
#D2
x1=[1.85,2.00,2.15]
x0=[0.50,1.00,1.35,2.65,3.00,4.00,5.00,6.00]
plot_filter_gain_prior(2,x1,x0,[(2.15,2.65),(1.85,1.35)],ax[2,0],sns.color_palette("deep",n_colors=7)[1])
#D3
x1=[2.85,3.00,3.15]
x0=[0.50,1.00,2.0,2.35,3.65,4.00,5.00,6.00]
plot_filter_gain_prior(3,x1,x0,[(2.85,2.35),(3.15,3.65)],ax[2,0],sns.color_palette("deep",n_colors=7)[2])
#D4
x1=[3.85,4.00,4.15]
x0=[0.50,1.00,2.0,3.0,3.35,4.65,5.00,6.00]
plot_filter_gain_prior(4,x1,x0,[(3.85,3.35),(4.15,4.65)],ax[2,0],sns.color_palette("deep",n_colors=7)[3])
#D5
x1=[4.85,5.00,5.15]
x0=[0.50,1.00,2.0,3.0,4.0,4.35,5.65,6.00]
plot_filter_gain_prior(5,x1,x0,[(4.85,4.35),(5.15,5.65)],ax[2,0],sns.color_palette("deep",n_colors=7)[4])
#D6
x1=[5.85,6.00,6.15]
x0=[0.50,1.00,2.0,3.0,4.00,5.00,5.35,6.65]
plot_filter_gain_prior(6,x1,x0,[(5.85,5.35),(6.15,6.65)],ax[2,0],sns.color_palette("deep",n_colors=7)[5])
ax[2,0].set_xticks([0,1,2,3,4,5,6,7,8])
ax[2,0].set_ylabel("Filter Band")
ax[2,0].set_xlabel("Frequency (cycles/tidal day)")
ax[2,0].set_xlim(-.25,6.5)
ax[2,0].set_ylim(-0.6,6.9)
ax[2,0].set_yticks([0,1,2,3,4,5,6])
ax[2,0].set_yticklabels(["D0","D1","D2","D3","D4","D5","D6"])
lines2 = ax[2,0].get_lines()[0:2]
filterlabels=["Pass","Stop"]
legend111 = ax[2,0].legend(lines2, filterlabels,bbox_to_anchor=[0.15,0.95],ncol=2)

                       
fig.tight_layout()
fig.savefig("filter_gain.png",dpi=300)
fig.savefig("filter_gain.pdf",dpi=300)
plt.show()








