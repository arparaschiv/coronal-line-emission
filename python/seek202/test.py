
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
plt.rcParams["figure.figsize"] = [5.,4.5]
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
lmax=8.5
lc = np.arange(lmin,lmax,1)
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
# here are the observations to be chose in fake_obs
#
g=cle.readhdr(dbdir+'db.hdr')
print(g)
g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax, \
        bthetamin, bthetamax, nline, wavel  = cle.parsehdr(g)
oindex=randint(0,ned*ngx*nbphi*nbtheta) # pick one realization
#oindex=715101
countsm=10.**lmin

for i in arange(0,5): 
    print('\n\n  I =  ',i, '<--------------------------------------------------\n')
    nline,sobs,orms,oindex = cle.obs_fake(dbdir,yobs,oindex,countsm)

