#!/usr/bin/python3

# The exact solution is graphed.  For the derivation of the solution,
# please see the handwritten notes.

from scipy.constants import *
from math import *

from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

# rewriting using Emma's notation (see second handwritten notes)

def spin_rot(k,alpha):
    sintheta=1/sqrt(1+k**2)
    costheta=k/sqrt(1+k**2)
    omegat=alpha*sqrt(1+k**2)
    sx=sintheta**2*np.cos(omegat)+costheta**2
    sy=sintheta*np.sin(omegat)
    sz=sintheta*np.costheta*(np.cos(omegat)-1)
    return sx,sy,sz

def spin(k,alpha): # non-rotating frame
    sintheta=1/sqrt(1+k**2)
    costheta=k/sqrt(1+k**2)
    omegat=alpha*sqrt(1+k**2)
    omega1t=alpha
    sx=(sintheta**2*np.cos(omegat)+costheta**2)*np.cos(omega1t)+sintheta*np.sin(omegat)*np.sin(omega1t)
    sy=-(sintheta**2*np.cos(omegat)+costheta**2)*np.sin(omega1t)+sintheta*np.sin(omegat)*np.cos(omega1t)
    sz=sintheta*costheta*(np.cos(omegat)-1)
    return sx,sy,sz

def spin_perfect(k,alpha): # non-rotating frame, but assuming perfect adiabatic
    sx=np.cos(alpha)
    sy=-np.sin(alpha)
    sz=0*alpha
    return sx,sy,sz

dim=3
fig,axs=plt.subplots(dim,sharex=True)

for d in range(dim):
    axs[d].set_xlim(0,pi/2)
    axs[d].set_ylim(-1,1)
    axs[d].set_xticks(np.arange(0,pi/2+pi/8,step=(pi/8)), ['0','$\pi/8$','$\pi/4$','$3\pi/8$','$\pi/2$'])

    
axs[0].set_ylabel('$S_x$')
axs[1].set_ylabel('$S_y$')
axs[2].set_ylabel('$S_z$')

N=100

alpha=np.linspace(0,pi/2,N)

colors=plt.cm.gist_yarg(np.linspace(0,1,20))


for k in range(2,14,2):
    sx,sy,sz=spin(k,alpha)
    axs[0].plot(alpha,sx,label='$k=%d$'%(k),color=colors[k])
    axs[1].plot(alpha,sy,label='$k=%d$'%(k),color=colors[k])
    axs[2].plot(alpha,sz,label='$k=%d$'%(k),color=colors[k])
    print(k,sy[-1])
    
k=0
sx_perfect,sy_perfect,sz_perfect=spin_perfect(k,alpha)
axs[0].plot(alpha,sx_perfect,label='$k=\infty$',color='red')
axs[1].plot(alpha,sy_perfect,label='$k=\infty$',color='red')
axs[2].plot(alpha,sz_perfect,label='$k=\infty$',color='red')


axs[0].legend(loc='lower left',ncol=4)
axs[1].legend(loc='upper left',ncol=4)
axs[2].legend(loc='upper left',ncol=4)

    
plt.show()
