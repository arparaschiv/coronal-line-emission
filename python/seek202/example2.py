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
# LOOP OVER DX
#  the number of counts which correspond to the
#  largest intensity in each of the nobs observations
#  for example: counts = 1.e8 gives rms fluctuations 10^-4 * max(intensity)
#
counts=1.e6
dmin=0
dmax=90
dc = np.arange(dmin,dmax,10)
nl=len(dc)
fr=dc*0.
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
# variables for plotting:
numsol=dc+0
chi=dc*0.
sol=dc*0.
#
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


start0=time.time()

kount=0
for dx in dc:
    print("\n")
    strc="{:0.1f}".format(np.log10(counts))
    nline,sobs,orms,oindex = cle.obs_fake(dbdir,yobs,oindex,counts,dx=dx)
    rms = orms
    i,j,k,l=cle.params(oindex,g)
    pp0=cle.physa(oindex,yobs,g,1)
    xindex=cle.invparams(i,j+dx,k,l,g)
    pp1=cle.physa(xindex,yobs,g,1)
    print("sobs", sobs[:4])
    #
    #**********************************************************************
    #
    #  Find solutions with chisq < 1.4 or so.
    #
    #**********************************************************************
    
    index,sfound,chisq,physics=\
        cle.clematch(sobs,yobs,rms,iplot=0,verbose=0,maxchisq=100,dbdir=dbdir)
    numsol[kount]=size(index)
    print("\n Observed fake index = ",index)
    print("sobs", sfound," found")
    pp =cle.physa( index,yobs,g,1)
    #print(pp.shape)
    p0="{:4.2f}".format(pp0[2])
    p1="{:4.2f}".format(pp1[2])

    print("x input",p0,p1,"dx=",dx,"physics ", pp[2])
    #print("chisq =",chisq)
    #print(size(index), "solutions found:",index)
    
    sol[kount]=pp[2]
    chi[kount]=chisq
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
fig1, ax1 = plt.subplots(2,1) # generates a Figure and Axis object
ax1.plot(dc,chi,'o')
ax1.set_title("")

ax1.set_ylabel('$\chi^2$  ')
ax1.set_xlabel('$\Delta x$ ')
filen='example_chi2_dx'
print("Saving ", filen+".pdf")
savefig(filen+".pdf")
plt.close('')
os.system("open "+filen+".pdf")


