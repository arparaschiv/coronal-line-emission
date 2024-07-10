#####################################################################
#  PROGRAM TO FIND BEST FIT PARAMETERS TO STOKES
#  VECTORS FROM M1 CORONAL LINES
#   P. Judge and A. Paraschiv 2020 December 18
######################################################################
#
# Needed libraries
#
import numpy as np
from scipy.io import FortranFile
from pylab import *
from matplotlib import pyplot as plt
from scipy import stats
import time
import sys
import os
import importlib
import cle as cle   # there is a cle.py that must be in the same
                    # directory as this file

np.set_printoptions(precision=3,suppress=False)
plt.rcParams["figure.figsize"] = [6.5,5]
matplotlib.rcParams.update({'font.size': 9})

######################################################################
#
# MAIN PROGRAM : some definitions
#
lab=['I','Q','U','V']
lab+=lab
s2d=180./np.pi
print("\n\n\n\n\n\n\n\n\n")
#
#######################################################################
# here are some variables that should be passed as parameters to a call
# to a routine
#
######################################################################
# LOOP OVER COUNTS
#  the number of counts which correspond to the
#  largest intensity in each of the nobs observations
#  for example: counts = 1.e8 gives rms fluctuations 10^-4 * max(intensity)
#
lmin=4.5
lmax=7.
lc = np.arange(lmin,lmax,1./3)
nl=len(lc)
fr=lc*0.
kount=0
#
start=time.time()
start1=start+0
yobs=1.37
#
#######################################################################
#
#  Read Observations: 
#   
#######################################################################
#
# sobs is the observed Stokes array (nline*4)
#
dbdir='dbcle/test/'
numsol=lc+0
kount=0
#
# here are the observations to be chosen in fake_obs
#
g=cle.readhdr(dbdir+'db.hdr')
print(g)
g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax, \
        bthetamin, bthetamax, nline, wavel  = cle.parsehdr(g)
oindex=randint(0,ned*ngx*nbphi*nbtheta) # pick one realization
#oindex=715101
countsm=10.**lmin

nline,sobs,orms,oindex = cle.obs_fake(dbdir,yobs,oindex,countsm)

start0=time.time()

for lcounts in lc:
    rms = orms * np.sqrt(10.**(lmin-lcounts) )
    counts = countsm*10.**(lcounts-lmin)
    print("\n")
    strc="{:0.1f}".format(np.log10(counts))
    print('Fake observations: log counts per state =',strc)
    print('S     ', [ "{:9.2e}".format(x) for x in sobs ])
    print('rms/S ',[ "{:9.2e}".format(x) for x in rms/sobs ])

    #
    #**********************************************************************
    #
    #  Find solutions with chisq < 1.4 or so.
    #
    #**********************************************************************
    
    index,sfound,chisq,physics=\
        cle.clematch(sobs,yobs,rms,iplot=0,verbose=0,maxchisq=1.4,dbdir=dbdir)
    numsol[kount]=size(index)
    print("\nFor observed fake index = ",oindex)
    #print(cle.phys(oindex,yobs,g))
    print("chisq =",chisq)
    print(size(index), "solutions found:",index)
    kount+=1
#
# end of example
#
dt= "{:4.3f}".format(time.time()-start0)
dtn= "{:4.3f}".format((time.time()-start0)/nl)
print(dt,' SECONDS FOR TOTAL ',dt, ' average per search ',dtn)
#
# Some examples of outputs as plots
#
fig1, ax1 = plt.subplots() # generates a Figure and Axis object
ax1.plot(lc,(numsol),'-')
ax1.set_title("")

i=20. # erg/cm2/s/sr
eff=0.04
px= 0.5 #" arcsec
f = (px*np.pi/180./60/60)**2 * i
hh=6.626e-27
cc=3.e10
lam=1.0747e-4
aperture=400.
countps = f * eff* (lam/hh/cc)*np.pi*(aperture/2.)**2
lgc=log10(countps)

plt.arrow(lgc,10,0,-5,head_width=.14,head_length=2)
lgc-=.4
plt.text(lgc,21,'CRYO-NIRSP')
plt.text(lgc,16,'counts/s/0.5"pixel for')        
plt.text(lgc,12,'I$_{1.0747}= 20$ millionths')

ax1.set_ylabel('Number of acceptable solutions ')
ax1.set_xlabel('Log$_{10}$ counts')
filen='example_stats'
print("Saving ", filen+".pdf")
savefig(filen+".pdf")
plt.close('')
os.system("open "+filen+".pdf")

alpha=0.3
filen='example_ambiguities'
fig, ax = plt.subplots(2) # generates a Figure and Axis object
plt.yscale("linear")
ax[0].set_xlabel('Physical parameter ')
ax[0].set_ylabel('Value ')
ax[0].set_title("Solution ambiguities")
#ax[0].set_yticks([-10,-3,0,3,10])
#plt.style.use('classic')
ax[0].autoscale(enable=True, axis='y', tight=None)
ax[0].xaxis.set_major_formatter(plt.FuncFormatter(cle.phys_axnames))



ax[1].set_xlabel('Stokes component $S_{i}$ ')
ax[1].set_ylabel('1-comp/obs (%)')
ax[1].set_title("")

nn=size(index)
if nn >= 2:
    plt.yscale("linear")
    
    for isol in np.arange(0,nn):
        strc="{:5.2f}".format(np.log10(chisq[isol]))
        ax[0].plot(physics[isol,:],'-o',alpha=alpha,markersize=(3 + isol*1.2),label="$\log_{10}\chi^2$="+strc)

    ax[0].legend(loc="lower left",fontsize=6)

    plt.yscale("linear")
    for isol in np.arange(0,nn):
        bvec=physics[isol,3:5]
        bmag=np.sqrt(dot(bvec,bvec))
        splot=sfound[isol,:]
        splot[[3,7]]*= bmag
        splot-=sobs
        splot[0]=0
        ax[1].plot(-splot*100,'o-',alpha=alpha,markersize=(3 + isol*1.2))
        ax[1].errorbar(np.arange(0,8),  -splot,rms*100,alpha=alpha)

    # recompute the ax.dataLim
    
    ax[1].autoscale(enable=True, axis='y', tight=None)
    plt.xlim(-1,8)
    fig.tight_layout()
    print("Saving ", filen+".pdf")
    savefig(filen+".pdf")
    plt.close()
    os.system("open -F "+filen+".pdf")

