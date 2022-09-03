#!/usr/bin/python3

from scipy.constants import *
from math import *

from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import animation

from matplotlib.patches import FancyArrowPatch
from mpl_toolkits.mplot3d import proj3d

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        FancyArrowPatch.__init__(self, (0,0), (0,0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def draw(self, renderer):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
        self.set_positions((xs[0],ys[0]),(xs[1],ys[1]))
        FancyArrowPatch.draw(self, renderer)

gamma=physical_constants['neutron gyromag. ratio'][0]
gamma=-gamma/1e6 # rad/s/uT
# gamma is positive in the physical constants library
# to make it negative, you have to add a minus sign
# I think I've corrected signs below so that either sign can be selected

# Set the field and rate at which it should change, here.
B1=1 # uT
T1=1 # s, time to go around once
omega1=2*pi/T1 # Hz
# At 1 uT, gamma*B1 = 2*pi*(30 Hz)
# So pick an omega1 that's a bit slower than that so we're fairly adiabatic
# I picked omega1=2*pi*(1 Hz) here.

Beff=sqrt(B1**2+(omega1/gamma)**2)
sintheta=omega1/abs(gamma*Beff)
costheta=B1/Beff

omega=gamma*Beff

# The 100 in the line below should make sure we get about 100 points
# in each of the revolutions in the rotating frame.
N=abs(int(omega/omega1*100))
print(N)

def spin_rot(t):
    sx=sintheta**2*np.cos(omega*t)+costheta**2
    sy=sintheta*np.sin(omega*t)
    sz=sintheta*costheta*(np.cos(omega*t)-1)
    return sx,sy,sz

def spin(t): # non-rotating frame
    sx=(sintheta**2*np.cos(omega*t)+costheta**2)*np.cos(omega1*t)+sintheta*np.sin(omega*t)*np.sin(omega1*t)
    sy=-(sintheta**2*np.cos(omega*t)+costheta**2)*np.sin(omega1*t)+sintheta*np.sin(omega*t)*np.cos(omega1*t)
    sz=sintheta*costheta*(np.cos(omega*t)-1)
    return sx,sy,sz


def spin_perfect(t): # non-rotating frame, but assuming perfect adiabatic
    sx=np.cos(omega1*t)
    sy=-np.sin(omega1*t)
    sz=0*t
    return sx,sy,sz

fig = plt.figure()
ax = fig.add_subplot(projection='3d')

t=np.linspace(0,T1,N)
sx,sy,sz=spin(t)
sx_perfect,sy_perfect,sz_perfect=spin_perfect(t)
ax.set_xlim(-1,1)
ax.set_ylim(-1,1)
ax.set_zlim(-1,1)

line,=ax.plot(sx[:1],sy[:1],sz[:1])
line2,=ax.plot(sx_perfect[:1],sy_perfect[:1],sz_perfect[:1])
a=Arrow3D([0,sx[0]],[0,sy[0]],[0,sz[0]],mutation_scale=20,arrowstyle="-|>",color="r")
a2=Arrow3D([0,sx_perfect[0]],[0,sy_perfect[0]],[0,sz_perfect[0]],mutation_scale=20,arrowstyle="-|>",color="b")
ax.add_artist(a)
ax.add_artist(a2)

def update(num):
    line.set_xdata(sx[:num])
    line.set_ydata(sy[:num])
    line.set_3d_properties(sz[:num])
    line2.set_xdata(sx_perfect[:num])
    line2.set_ydata(sy_perfect[:num])
    line2.set_3d_properties(sz_perfect[:num])
    a._verts3d=[0,sx[num]],[0,sy[num]],[0,sz[num]]
    a2._verts3d=[0,sx_perfect[num]],[0,sy_perfect[num]],[0,sz_perfect[num]]
    return [line,line2,a2,a]

ani=animation.FuncAnimation(fig,update,frames=N,interval=1,blit=False,repeat=True)
ani.save('asr.mp4')
plt.show()
