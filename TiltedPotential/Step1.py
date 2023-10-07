# -*- coding: utf-8 -*-
"""
Created on Mon Jul  4 14:22:25 2022

@author: jain_
"""

# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
In this program, we are going to look at how to calculate the 
relationship between "stuck time" and the "potential depth" in a 
non-equilibrium situation. Here, a 2 micrometer sized particle is stuck in a 
circular trap of radius ~8um. On this trap, there is a tangential drive force 
and a cosine potential of depth of the order kbT created by an optical tweezer 
and inside water medium. The drive force is of the order of ~1 pN. 
"""
import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
# plt.style.use('ggplot')
# import sympy as sp
import math
import time
#%%
nums=int(1e6) # Simulation steps
dt =100     # Delta t of the simulation (Could be variable)
nump=20       # number of particles
F0 = 2.1e-15# Drive force
T=300         #T #Kelvin is 25degreeCelsius
eta=8.9e-4    #eta  #kg m^{-1}s^{-1} Dynamic Viscosity of water
pi=np.pi;    
a=2e-6        # Micrometers Size of particle
r=10e-6        #Radius of the outer circle
periodR = 2*pi*r;
k_b=1.38e-23
zeta=3*6*pi*eta*a  #A factor 3 is needed to account for the plate interaction
m=4.0*pi/3*a**3*1100 #density is 1100kg/m^3
kBT = k_b*T;
#%
#We calculate V_p as, the time it takes to cover a circle of radius 30 um in 1 second
V_p=20-6; #micrometers per second
# Initial Conditions: Potential well

rMat= np.linspace(0.0,periodR,1000);
thetaMat = np.linspace(0.0,2*pi,1000)
#Potential Function
Ufunc_ = 10*k_b*T*np.cos(thetaMat) 
Ufunc_ = Ufunc_- max(Ufunc_)
Uforce_ = -np.diff(Ufunc_)/np.diff(rMat)


fig, axs = plt.subplots(3)
fig.suptitle('Potential, force, and tilted potential')
axs[0].plot(rMat,Ufunc_)


rMat_ = np.linspace(0.0,periodR,999)
axs[1].plot(rMat_,Uforce_)

from scipy.interpolate import interp1d

Ufunc = interp1d(rMat,Ufunc_)
Uforce = interp1d(rMat_,Uforce_)

axs[2].plot(rMat/periodR*2*pi,Ufunc_-F0*rMat)
# We calculate deltaU and the minima position
