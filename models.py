#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct  5 04:51:40 2020

@author: Ayca
"""
import math
import scipy 
#%% dB-dBm
def mW2dBm(x):
    return 10*math.log(x,10)
def mW2dB(x):
    return x-30 #10log(10^3)
#%% Doppler Shift: Transmitter or Receiver movements in changein_t
# causes changein_d since TX signal need to travel to reach the receiver
def doppler(v,changein_t,angle,fc):
    changein_d = v*changein_t*math.cos(math.radians(angle)) #change in wireless path
    wave_lambda = scipy.constants.c/fc
    changein_ph = 2*math.pi*changein_d/wave_lambda #phase change due to path difference
    fd = changein_ph/(2*math.pi*changein_t) #Doppler Frequency
    print("Doppler freq alternative works",fd ==v*math.cos(math.radians(angle))/wave_lambda)
    print("Max Doppler freq =",v/wave_lambda)
    return changein_d,wave_lambda,changein_ph,fd

#%% Friis TX Formula
def FriisTX(pt,gt,gr,wave_length,d):
    pr = pt*gt*gr*math.pow(wave_length,2)/math.pow((4*math.pi*d),2)
    pr_dBm = mW2dBm(pt)+mW2dBm(gt)+mW2dBm(gr)+mW2dBm(math.pow(wave_length,2)) - mW2dBm(math.pow((4*math.pi*d),2))
    pt_dBm = mW2dBm(pt)
    pl = pt/pr
    pl_dBm = pt_dBm - pr_dBm
    print("Pr: {}, Pr[dBm]: {}".format(pr,pr_dBm))
    print("Free-Space Path Loss: {}, Free-Space Path Gain: {}".format( pl,-pl))
    print("Pl[dBm]: {}, Pg[dBm]: {}".format( pl_dBm,-(pl_dBm)))
    return pr,pr_dBm,pl,pl_dBm
def FSPL(pt,gt,gr,wave_length,d):
    return FriisTX(pt,gt,gr,wave_length,d)[2]
#%%Okumura Model
#fspl(fc,d) - free space loss at distance d, freq fc
#A(fc,d) - additional attenuation averaged across all environments
#gArea - environment specific gain
def Okumura(fspl,hr,ht,gArea,attenuation): #ht&hr in meters
    if ht<30 or ht>1000 or hr>10:
        return
    ght = 20*math.log10(ht/200)
    if hr<=3:
        ghr =10*math.log10(hr/3)
    elif hr>=3 and hr<=10:
        ghr = 20*math.log(hr/3)
    pl = fspl + attenuation -ght -ghr-gArea
    pl_dBm = mW2dBm(fspl)+mW2dBm(attenuation)-mW2dBm(ght)-mW2dBm(ghr)-mW2dBm(gArea)
    return pl,pl_dBm

#%%Simplified Path Loss Model
def SimplifiedPL(pt, fc,d0,d,a):
    k_dB = 10*math.log10((scipy.constants.c/fc)/(4*math.pi*d0))
    pr_dBm = mW2dBm(pt)+k_dB-10*a*math.log10(d/d0)
    return pr_dBm
    
    