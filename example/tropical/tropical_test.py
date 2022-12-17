# -*- coding: utf-8 -*-
"""
Created on Tue May 31 14:26:50 2022

@author: eli
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pywt



def tropical_test_calc(B):
    wspec= f'cmor{B}-0.96'
    print(wspec)
    trop = 2.*13.66*24.   # tropical period in hours
    
    
    
    A12 = 0.75
    A22 = 0.5
    twopi = 2.*np.pi
    A = np.array([0.2,0.8,1.,0.2,0.2])  # amplitudes
    A1adj = np.sqrt((0.5*A12*A[1]*A[2])**2. + A[1]*A[1])   # adjusted D1 to include contribution from D1*D2
    phi1adj = np.arctan(0.5*A12*A[2])
    w = np.array([0.,1.,2.,3.,4.])*0.966137/24. # frequencies in c/hr
    t = np.arange(0.,800,0.25)  # time in hours
    mod = np.cos(twopi*t/trop)  # modulation over tropical cycle
    
    A[3] = A12*A[1]*A[2]*0.5
    A[4] = -A22*A[2]*A[2]*0.5
    
    D0 = mod*A[0]
    D1_orig = mod*A[1]*np.cos(twopi*w[1]*t)
    D1 = mod*A1adj*np.cos(twopi*w[1]*t - phi1adj)
    D2 = A[2]*np.cos(twopi*w[2]*t - np.pi/2.)
    D3 = mod*A[3]*np.cos(twopi*w[3]*t - np.pi/2.)
    D4 = A[4]*np.cos(twopi*w[4]*t )
    
    y0 = D1_orig + A12*D1_orig*D2
    y1 = D1+D3
    
    y0 = D0 + D1_orig + D2 + A12*D1_orig*D2 + A22*D2*D2
    D0_change = 0.5*A22*A[2]*A[2]
    y1 =   D0 + D0_change + D1 +D2 + D3 +  D4
    y1 =   D0 + D1 +D2 + D3 +  D4
    
    print(y0[0:15])
    print(y1[0:15])
    
    noise = np.random.normal(scale=0.01,size=len(t))
    Dtest=[]
    for i in range(5):
        Dtest.append(np.cos(twopi*w[i]*t))
    ynoisy = Dtest[4]
    #ynoisy = ynoisy*0
    #ynoisy[900]=1.
    ynoisy = y1
    test_tide = pd.DataFrame(index=t,data=y1)
    test_tide.to_csv("test_tide.csv",sep=",",header=True)
    test_tide.plot()
    plt.show()
    
    amplitudes = pd.DataFrame(index=t)
    amplitudes["D1"] = abs(mod*A1adj)
    amplitudes["D3"] = abs(mod*A[3])
    #df.plot()
    
    
    wav = pywt.ContinuousWavelet(wspec)
    
    # print the range over which the wavelet will be evaluated
    print("Continuo1s wavelet will be evaluated over the range [{}, {}]".format(
        wav.lower_bound, wav.upper_bound))
    
    width = wav.upper_bound - wav.lower_bound
    print(f"Width: {width}")
    scales = np.array([12,24,32,48,96])
    #scale = 1 corresponds to 4 hours
    # 2=8 hours  
    dt=1/96.  # 15 minute sample duration expressed in days
    sampling_freq = 1./dt  # sampling freq in samples/day
    f = pywt.scale2frequency(wav, scales)/dt
    print("first freq calc",f)
    
    
    cwtmatr, freqs = pywt.cwt(ynoisy,scales, wavelet=wspec,sampling_period=dt)
    
    print("freqs",freqs)
    print("periods",1/freqs)
    #cwtmatr = np.absolute(cwtmatr)
    norms = [1.,4.891,3.457,2.8222,2.440]
    for ifreq in [1,3]:
        #norms.reverse()
        normalize = norms[ifreq]
        cwtnorm = cwtmatr/normalize
        #print(np.absolute(cwtmatr).sum(axis=1))
        lo=400
        hi=2200
        matndx = 5-ifreq
        plt.plot(t[lo:hi],cwtnorm[matndx,lo:hi].real)
        amplitude = np.absolute(cwtnorm[matndx,:])
        components = {1:D1,2:D2,3:D3,4:D4}
        comp = components[ifreq]
        plt.plot(t[lo:hi],comp[lo:hi])
        plt.plot(t[lo:hi],amplitude[lo:hi])
        plt.legend(["wavelet","signal","amplitude"])
        plt.show()
        amplitudes[f"B{B}_D{ifreq}"]=amplitude
        
    amplitudes.plot()
    amplitudes.to_csv(f"cmor_{B}-0.96_result.csv")
    #widths = np.array([6,8,12,25,49,99])
    #cwtmatr = signal.cwt(res.values,signal.morelet,widths)
    #plt.imshow(cwtmatr, extent=[-1, 1000, 5, 69], cmap='PRGn', aspect='auto',
    #           vmax=abs(cwtmatr).max(), vmin=-abs(cwtmatr).max())  # doctest: +SKIP


for B in [2.5,4,8]: tropical_test_calc(B)