import os
import numpy as np
from scipy.io import FortranFile
from pylab import *
from astropy.io import ascii
from ast import literal_eval
from matplotlib import pyplot as plt
from scipy import stats
import matplotlib.ticker as mtick
import matplotlib.pylab as pylab
from astropy.table import Table, join
import array as arr
import warnings
import pickle
import time


y=np.arange(0.,2,0.005)
z=np.arange(-20,20,0.1)

th=np.arctan2(z,y)
#######################################


fig=figure()

plt.rc('text',usetex='True')

plt.subplots_adjust(wspace=0.,hspace=0.,bottom=0.,left=0.)

plt.subplot(1,2,1)

plt.axis('equal')
plt.axis('off')
plt.xlim(-0.3,3)
#plt.title('(a)')

plt.plot(np.cos(th), np.sin(th),'b--',linewidth=0.7)

# angle annotation and curve
plt.text(0.24,0.04,r'$\alpha$')
plt.plot(0.4*np.cos(th[200:206]), 0.4*np.sin(th[200:206]),'k')


plt.plot(y,z*0,'k-')
plt.text(y[-1],0.02,'P')

plt.plot(0.9*y,y/2,'k-')
plt.text(0.9*y[-1],y[-1]/2+0.02,'Q')

## lvs for p and q directions
plt.text(1.14,0.03,'$\mathit{proj. l.v.s}$',rotation=0,fontsize=5)
plt.text(1.01,0.62,'$\mathit{proj. l.v.s}$',rotation=31,fontsize=5)
#
# coordinate orientation x and y
# 
plt.arrow(.0,.0,0.5,0.,head_width=0.06,zorder=5)
plt.text(.05,-0.2,'$y$')

plt.arrow(.0,.0,0.,0.5,head_width=0.06)
plt.text(-.18,0.1,'$z$')

plt.text(-0.13,-0.13,'$O$')
plt.text(0.01,-1.3,'(a) Plane of the sky (POS) $zy$; $x=0$ ',fontsize=10)

#plt.text(1.2,1,'(a)')


#######################################
# second plot
##
plt.subplot(1,2,2)

#plt.title('(b)')

plt.axis('equal')
plt.axis('off')
plt.xlim(-1.2,1.2)

x=np.arange(-1,1)

# x and y axes
plt.plot([0,0.0001],[0.,0.0001],'w')
#plt.arrow(1,-0.3,-2,1,color='k',head_width=0.05)  ##leftarrow LVS
plt.arrow(-1.0,0.7,2,-1,color='k',head_width=0.05)  ##rightarrow LVS
#plt.text(-1,0.82,'To Sun center',rotation=-30)
plt.text(+0.47,-0.22,'$\leftarrow$to Sun center',rotation=-27,fontsize=6)
#plt.arrow(-0.74,0.62,-0.05,0.025,color='k',head_width=0.02,width=0.001) #Small arrow to sun center
# angle annotation and curve
xc=0.0
yc=0.2
plt.text(xc+0.25,yc-0.07,r'$\vartheta_B$') #r
plt.text(xc+0.18,yc-0.27,r'$\vartheta_B$')  #s
plt.text(xc-0.34,yc-0.01,r'$\vartheta_B$') #p
plt.text(xc-0.29,yc+0.15,r'$\vartheta_B$')  #q
#plt.plot(0.3*np.cos(th[206:]), 0.3*np.sin(th[206:]),'k')


#
# now line plots for other parts of cone
#
x=xc+0.85*np.array([0,0.2])
yy=yc+0.85*np.array([0,-.5])
plt.plot(x,yy,'k:')


x=xc+0.85*np.array([0,-0.2])
yy=yc+0.85*np.array([0.,.5])
plt.plot(x,yy,'k:')
#
#
cenx=-1
ceny=.2
# Original axis
# plt.arrow(cenx,ceny,0.3,0.,head_width=0.03)
# plt.text(cenx+.25,ceny+0.1,'y')

# plt.arrow(cenx,ceny,0.,-0.3,head_width=0.03)
# plt.text(cenx-.17,ceny-0.3,'x')
# axis
plt.text(cenx-.1,ceny+0.5,'$\mathit{O}$',fontsize=10)
plt.arrow(cenx,ceny+0.5,0.33,0.,head_width=0.03)
plt.text(cenx+.06,ceny+0.56,'$\mathit{y}$',fontsize=10)

plt.arrow(cenx,ceny+0.5,0.,-0.3,head_width=0.03)
plt.text(cenx-.1,ceny+0.4,'$\mathit{x}$',fontsize=10)
plt.text(cenx-.08,ceny+0.15,'$\hat\mathbf{k}$',fontsize=7)

plt.text(cenx-0.24,ceny-0.0,'$\leftarrow$ to observer',rotation=90,fontsize=7)


#
#######################################################################
# now plot angles around the b vector
#######################################################################
#  first one
#
alpha = np.arange(100)/99.
mn=81* np.pi/180
mx=154* np.pi/180
alpha = mn + (alpha)*(mx-mn)
plt.plot(xc+np.sin(alpha)/5.4,yc+np.cos(alpha)/5.4,'k',linewidth=.7)
#
#reflect this 
#
alpha+=np.pi
plt.plot(xc+np.sin(alpha)/5.9,yc+np.cos(alpha)/5.9,'k',linewidth=.7)
#plt.plot(xc+np.sin(alpha)/2+.65,yc+np.cos(alpha)/2-.27,'k')

plt.arrow(xc,yc+0.006,0.5,0.106,color='k',head_width=0.05)
plt.text(xc+0.48,yc+0.18,'$\hat\mathbf{b}$',fontsize=7)
plt.text(xc+1.01,yc-0.47,'$\hat\mathbf{r}$',rotation=-27,fontsize=7)
plt.text(xc+0.54,yc-0.45,'$l.v.s$',rotation=-27,fontsize=7)

x=xc+0.85*np.array([0,-0.5])
yy=yc+0.85*np.array([0,-0.1])
plt.plot(x,yy,'k:')


##sun angle
alpha = np.arange(100)/99.
mn=25* np.pi/180
mx=205* np.pi/180
alpha = mn + (alpha)*(mx-mn)
plt.plot(-0.99+np.sin(alpha)/1.8,0.66+np.cos(alpha)/1.8,'b--',linewidth=.7)
#
#######################################################################
# now plot angles around the b vector
#######################################################################
#  second one
#
alpha = np.arange(100)/99.
mn=90* np.pi/180
mx=138* np.pi/180
alpha = mn + (alpha)*(mx-mn)
#
#reflect this 
#
alpha+=np.pi
#plt.plot(xc+np.sin(alpha)/2,yc+np.cos(alpha)/2,'k')
alpha-=np.pi
#plt.plot(xc+np.sin(alpha)/2-.65,yc+np.cos(alpha)/2+.27,'k:',linewidth=.5)

plt.text(-1,-.48,'(b) Ecliptic plane $xy$; $z=0$ ',fontsize=10)


#
# define p q r s

plt.text(xc+0.02,yc+0.05,'o')
plt.text(xc-.42,yc-.17,'b$_3$',fontsize=6)
plt.text(xc-.1,yc+.37,'b$_4$',fontsize=6)
plt.text(xc+.5,yc+.01,'b$_1$',fontsize=6)
plt.text(xc+.04,yc-.4,'b$_2$',fontsize=6)
fig.tight_layout()
fig.set_size_inches(6.2,2.7) 
savefig('figg.pdf',transparent=True, bbox_inches='tight', pad_inches=0)

#plt.close('all')

#os.system("open figg.pdf")


