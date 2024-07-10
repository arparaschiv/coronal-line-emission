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

i=10. # erg/cm2/s/sr
eff=0.04
px= 0.5 #" arcsec

f = (px*np.pi/180./60/60)**2 * i

hh=6.626e-27
cc=3.e10
lam=1.0747e-4
aperture=400.

countps = f * eff* (lam/hh/cc)*np.pi*(aperture/2.)**2

print(i,f,countps/1.e6, ' counts per second per pixel')
