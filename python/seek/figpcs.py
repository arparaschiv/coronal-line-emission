import os
import numpy as np
from scipy.io import FortranFile
from pylab import *
from astropy.io import ascii
from matplotlib import pyplot as plt
from scipy import stats
import matplotlib.ticker as mtick
import matplotlib.pylab as pylab
from astropy.table import Table, join
import array as arr
import warnings
import pickle
import time
import cle

os.system('ls sdata*pkl')
print('There are lines, enter number before .pkl to process that number')
ans=input()
ans="1"
sline=ans




pfile='pcadata'+sline+'.pkl'
print(pfile, ' read...')

u,sig,coeff = pickle.load(open(pfile,"rb"))

sh=u.shape
for nc in range(0,sh[0]):
  if sig[nc] < sig[0]/1.e5: break
print("Need only ",nc," principal components, next evalue is ",sig[nc])


pfile='pcagrid'+sline+'.pkl'
ned,ngx,ngy,nbphi,nbtheta,\
            xed,gxmin,gxmax,gymin,gymax,\
            bphimin,bphimax,bthetamin,bthetamax,\
            nline,nq = pickle.load(open(pfile,"rb"))
#
# plotting
#
nline=int(ans)

fig=figure(figsize=[8,5])
plt.rc('text',usetex='True')

plt.subplots_adjust(wspace=.2,hspace=.4,bottom=0.1,left=0.1)

sig/=sig[0]

nx=3
ny=3

ss=['I','Q','U','V']
ss=ss+ss
sl=['$1074.7\longrightarrow$',' ',' ',' ','$1079.8\longrightarrow$','','','','']


if(nline <2):
  ny=2
  sl=sl[:4]
  ss=ss[:4]

plt.subplot(ny,nx,1)

plt.plot(np.log10(sig[:10]),'k.')
plt.title("$\log_{10}\Sigma_{i=1\ldots 10}$")




for j in range(0,ny*nx-1) :
    plt.subplot(ny,nx,j+2)
    plt.plot(u[:,j],'k')
    ddd="{:0.2f}".format(np.log10(sig[j]))
    plt.title(ddd)
    axes=plt.gca()
    ylim=axes.get_ylim()
    for k in range(0,nline*4):
      plt.plot([nq[0]*k,nq[0]*k],[ylim[0]+0.1*(ylim[1]-ylim[0]),ylim[1]],'k',linewidth=.5)
      if j ==0: plt.text(nq[0]*k+5, ylim[0] + (ylim[1]-ylim[0]),ss[k])
      if j == 0 & (k == 4 or k == 0): plt.text(nq[0]*k+5, ylim[0] - (ylim[1]-ylim[0])*0.1,sl[k])
    plt.axis('off')

savefig("pca"+sline+".pdf")
plt.close('all')

os.system("open pca"+sline+".pdf")
