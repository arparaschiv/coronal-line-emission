import astropy as astro
import numpy as np
from pylab import *
import matplotlib.pyplot as plt
import statistics as stat
from matplotlib import gridspec
import scipy
import scipy.constants as const
from scipy.io.idl import readsav
import datetime
import matplotlib.dates as mdates
import csv
## to import the custom pa library in ipython
import sys
sys.path.insert(1, '/home/alin/Documents/physics_prog/cle/python/noise_model')
import pa

######################################################################
#
# NOISE MODEL FOR STOKES VECTORS IN THE CORONA
# 
#
######################################################################

np.set_printoptions(precision=3)

#
print(" ")
print(" ")
print("Order of magnitude Muller Matrix elements from Fig 13 of ")
print(" Harrington and Sueoka 2017 Journal of Astronomical Telescopes")
print(" Instruments, and Systems, Volume 3, id. 018002 (2017). ")
print(" The values chosen are read roughly from central points of each element")

mm= np.array([    [1.,.005,-0.01,0.02],\
                  [-0.01,-1.,1.,-0.1],\
                  [0.005,-1.,-1,0.04],\
                  [0.02,-0.1,0.05,1.] \
    ])


tel=pa.telescope('dkist',4.,mm,0.1,0.)

# cryo data from https://www.nso.edu/telescopes/dkist/instruments/cryo-nirsp/

aa=pa.analysis('n-s','si',10746.8,4,1,0,1)
dual=aa.dual
dual=1

sp=pa.spectrograph('cryo-nirsp',42, 'milli-Angstrom')   ## is this the spectral resolution bin? changed 134 to 42, which is corresponding to the 44A range spread onto 1024 pixels. 

dt=0.1
det=pa.detector('cryo-detector','n/a',dt,30,0.7,1,1,0.5,18.e-6,100000.) ##ARP changed maxdn 64000 to 100000 e- as per ipc; changed pxsize 4.0 to 1.8microns as this is the quoted pitch?

#####################################################################################
# fake data from solar corona
#####################################################################################
l0=10746.8e-10
doppsamp = sp.bin * 1.e-13 / l0 * const.c /1.e3 
x=np.arange(-50,50,doppsamp)*1.e3
dl=x/const.c*l0
shift=-4.e3
width=18.e3
xx= (x-shift)/width
int_disk=pa.iplanck(l0,6000.)  # w/m2/m/sr
i= 20.*exp(-xx**2)*int_disk*1.e-5
q=0.1*i
u=-0.01*i
v=np.gradient(i,x)*1.e1
#print ('bin in mA = %.1f ' % (l0*1.e13*doppsamp*1.e3/const.c))

i2cgs=1.e3*1.e-10  # erg/cm2/s/sr/a

######################################################################
#
# main integration loop
#
######################################################################
print("------------------------------------------------------------")
print("Stokes I emerging from the corona cgs")
print((i*i2cgs))

s=np.array([i,q,u,v]) # output from Sun
l0=1.07468e-6
sc=pa.sc1(l0,0)

print(" ")
print("------------------------------------------------------------")
print ("ATMOSPHERE ")
print ("Steady atmos scattering I  %0.2e " % sc[0])
tmp=sc[0]*i2cgs
print ("Steady atmos scattering I  %0.2e cgs " % tmp)
M,N = s.shape
sc=np.array(N*[sc]) # equivalent to rebin 
sc=np.swapaxes(sc,0,1)
################################################################
# sum of steady component and solar data
################################################################
ssun=s
s=s+sc

sperfect=s

######################################################################
#
#  Noise model
#
######################################################################
# 
# power spectrum and variances

p=pa.wlpspec(l0)  # power spectrum 

nu=p[0]
ps=p[1]
print ( "RMS all power above  0 Hz  %0.2e " % np.sqrt(pa.trapez(nu,ps)) )
idx=np.where(nu >= 1./det.readt)
wlrms=np.sqrt(pa.trapez(nu[idx],ps[idx]))
print ( "RMS power above 10 Hz      %0.2e " % wlrms )
print ( "RMS power above 10 Hz      %0.2e cgs " % (wlrms*i2cgs) )
print("------------------------------------------------------------")

ip=0
while ip == 1:
    f = figure(figsize=[9,7])
    plt.figure()
    plt.rc('text',usetex='True')
    plt.yscale('log')
    plt.xscale('linear')

    plt.text(1,1,"Variance $>$ 10 Hz " + "{:.2e}"%wlrms )
    
    plt.plot(p[0],p[1],color='blue')
    plt.title("Power spectrum model  SI units")
    
#    savefig("pspec.pdf")
    ip=0




######################################################################
#
#  Apply telescope matrix
#
######################################################################
# smod is the S array at the entrance to the Cryo modulator
smod=tel.mm.dot(s)
#
print("Telescope Mueller Matrix")
print(mm)
invm=np.linalg.inv(mm)
print("Inverse")
print(invm)


s=s*i2cgs
ssun=ssun*i2cgs
smod=smod*i2cgs

ip=0
while ip == 1:

    f = figure(figsize=[16,13])

    plt.figure()
    plt.rc('text',usetex='True')
    
    plt.subplots_adjust(wspace=0.3,hspace=0.3)
    
    plt.subplot(2,2,1)
    plt.plot(dl,s[0],color='blue')
    plt.plot(dl,ssun[0],color='red')
    plt.plot(dl,smod[0],color='grey')
    
    plt.title("I")
    
    plt.text(-2.2e-10,26,"Steady light",color='black')
    plt.text(.7e-10,26,"sun only",color='red')
    plt.text(.7e-10,23,"sun plus sky",color='blue')
    plt.text(.7e-10,20,"at Cryo mod",color='grey')
    
    
    plt.subplot(2,2,2)
    plt.ylim(-3,3)
    plt.plot(dl,s[1],color='blue')
    plt.plot(dl,ssun[1],color='red')
    plt.plot(dl,smod[1],color='grey')
    plt.title("Q")
    
    plt.subplot(2,2,3)
    plt.ylim(-3,3)
    plt.title("U")
    plt.xlabel("$\Delta\lambda$ m")
    #plt.ylabel("Brightness (W/m2/s/sr) ")
    plt.ylabel("Brightness (erg/cm2/s/A/sr) ")
    plt.plot(dl,s[2],color='blue')
    plt.plot(dl,ssun[2],color='red')
    plt.plot(dl,smod[2],color='grey')
    plt.title("U")
    
    plt.subplot(2,2,4)
    plt.plot(dl,s[3],color='blue')
    plt.plot(dl,ssun[3],color='red')
    plt.plot(dl,smod[3],color='grey')
    plt.title("V")
#    savefig("steady.pdf")
    ip=0

#
# Back to SI units
s=s/i2cgs
ssun=ssun/i2cgs
smod=smod/i2cgs



#
######################################################################
#
#  Now start integrations and include noise from two sources
#  1. photon counting at detector of Cryo-nirsp
#  2. random polarized fluctuations from Elmore's time series
#     of coronal wl fluctuations at Maunda Loa 2001 day 282 and 384
#
######################################################################
#
# integrate over ns seconds
#
ns=300
print(" ")
print("------------------------------------------------------------")
print("integration over %0.1f seconds " % ns )
nstates=aa.nstates
ncycles=round(ns/det.readt/nstates)
print("Number of full modulation cycles to integrate = %4i " % ncycles)
print("------------------------------------------------------------")

print(" ")
      
nw=len(dl) # number of wavelengths

# we have a 4 times n array as input, this needs
# noise adding in each time integration interval
# i.e. for each state

epp= const.h*const.c/l0
b=sp.bin * 1.e-10 * 1.e-3  # wavelength bin in metres
photons = tel.syseff * tel.area() * det.omega * b / epp 

#
jdet = pa.iplanck(l0,77.)*det.pixsz**2*4.*const.pi # thermal photons detector
cdet = jdet *det.qeff /epp
print(" read noise counts =%0.1f thermal noise =%0.1e " % (det.readn ,  cdet) )
#

np.set_printoptions(precision=2)

state = np.zeros([(dual+1)*nstates,nw])  # polarization modulation states

wt=pa.modwt4(nstates,dual)

#
#
wlamp = photons* wlrms * det.readt/2.

for ic in range(0,ncycles):
    #print("Cycle %4i " % ic)
    #
    # time-dependent noise is not dependent on wavelength
    # but only on time.
    #
    # deterministic  counts per read, /2. for each beam
    #  counts is an (nstate , nw) matrix.
    #  each of the counts is subject to random read and shot noise
    #  external noise is the same for each state.
    counts = photons *det.readt/2. * wt.dot(smod)
    #print(" shape of counts is " )
    #print(shape(counts))
    # 
    # Random (shot and read noise)
    #
    sigma = sqrt(det.readn+jdet) # standard deviation readnoise
    # 
    r = np.random.normal(0., sigma, [(dual+1)*nstates,nw])  # readnoise
    counts = counts+det.readn+r
    sigma=sqrt(counts) # shot noise
    r = np.random.normal(0., sigma,shape(counts))  # shot noise
    counts = counts+r
    #
    # noise from measured power spectra at Mauna Loa
    # first compute noise for each read state
    #
    countsn= np.random.normal(0., wlamp, nstates)  # MLSO noise, just for 1 cycle
    if(dual == 1):
        countsn = np.append(countsn, countsn)
    M,N = shape(counts)
    countsn=np.array(N*[countsn]) # equivalent to rebin 
    countsn=countsn.transpose()
    counts=counts+countsn

    state = state + counts

#print("STATE %4i cycle" % ic)
#print(state)

import matplotlib.ticker as mtick




ip=0
while ip == 1:
    f = figure(figsize=[9,7])

    fig = plt.figure()

    plt.rc('text',usetex='True')

    plt.plot(state[0],color='black')
    plt.plot(state[1],color='blue')
    plt.plot(state[2],color='yellow')
    plt.plot(state[3],color='red')
    if(dual == 1):
        plt.plot(state[4],'k.')
        plt.plot(state[5],'b.')
        plt.plot(state[6],'y.')
        plt.plot(state[7],'r.')

    plt.title("States")

    plt.plot()
#    savefig("states.pdf")
    ip=0

np.set_printoptions(precision=2)
#
#  Now demodulation
#

wtt=wt.transpose()
m=wtt.dot(wt)
y=wtt.dot(state)
output = np.linalg.solve(m,y)/ncycles
#
# Knowledge of Telescope matrix is given here
# accuracy. here is one realization

invm = invm #- 0.1*invm

output = invm.dot(output)
diff=(s-output)


print("------------------------------------------------------------")
print(" RMS (S-measured) / max (abs(S)) ")
names=['I','Q','U','V']
for j in range(0,4):
    mx=np.amax(s[j])
    sd=np.std(diff[j]) / np.fabs(mx)
    print (" %s     %0.2e " % (names[j],sd) )
print("------------------------------------------------------------")


import matplotlib.ticker as mtick

import matplotlib.pylab as pylab
params = {'legend.fontsize': 'x-large',
         'axes.labelsize': 'x-large',
         'axes.titlesize':'x-large',
         'xtick.labelsize':'large',
         'ytick.labelsize':'large'}
pylab.rcParams.update(params)

ip=1
while ip == 1:

    s1=np.sum(s[3])
    output=output / ( photons *det.readt/2.)
    f = figure(figsize=[8,6])
    f.text(0.12,0.96,"Integration time: "+str(ns)+" seconds",color='black')
    f.text(0.12,0.93,"Integrated V (pure):"+format(np.sum(np.abs(s[3])),'.2e')+" and V(obs) "+format(np.sum(np.abs(output[3])),'.2e'),color='black')
    f.text(0.44,0.03,"Wavelength [A]",color='black',fontsize=12)
    f.text(0.01,0.35,"Brightness [W/m2/m/sr]",color='black',fontsize=12,rotation=90)
    for j in range(1,5):

        ax = f.add_subplot(2,2,j)
        plt.subplots_adjust(wspace=0.3,hspace=0.3)
        ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.0e'))
        ax.xaxis.set_major_formatter(mtick.FormatStrFormatter('%.0f'))
        ax.plot(10746.8+dl*1e10,s[j-1]/10,linewidth=2)  
        ax.plot(10746.8+dl*1e10,output[j-1]/10,linestyle='dashed')  
#        ax.set_xlabel("Wavelength [m]")
        for axis in ['top','bottom','left','right']:
            ax.spines[axis].set_linewidth(2)
        plt.title(names[j-1])

    savefig("final_300s.pdf")
    ip=0

#import os
#os.system("open states.pdf  final.pdf")
#if ip == 1: os.system("open steady.pdf pspec.pdf")
