import numpy as np
from pylab import *
import scipy
import scipy.constants as const
from scipy.io.idl import readsav
import datetime
import matplotlib.dates as mdates
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.modeling.blackbody import blackbody_lambda, blackbody_nu
import csv
#import pa

#####################################################################################
#  FUNCTIONS
#####################################################################################

def area(r):
    return 3.141592654 * r**2


def iplanck(lam,t):
    """Planck intensity in SI units
    lam := wavelength m (required)
    t   := temperature K (required)
    output := W/m^2/m/sr
    """
    x= const.h*const.c/(lam*const.k*t)
    i= 2.*const.h*const.c**2/lam/lam/lam/lam/lam / (exp(x)-1.)
    return i




def modwt4(rph,dual):

    x=[[1.,1.,1.,1.], \
        [1.,1.,0.,0.], \
       [1.,0.,1.,0.], \
       [1.,0.,0.,1.]] 
    if(dual ==1):
        y=[[1.,1.,1.,1.], \
           [1.,-1.,0.,0.], \
           [1.,0.,-1.,0.], \
           [1.,0.,0.,-1.]]
        x=x+y # add lists
    
    x=np.array(x)
    print(" ")
    print('------------------------------------------------------------')
    print("Modulation weights  for %2i beam(s) "  % (dual+1) )
    print('------------------------------------------------------------')
    print(x.shape)
    print(x)
    return x




def mrot(theta):
    """ 
    Gets rotation matrix for Stokes vectors when 
    reference direction is changed by angle theta
    use: 
    m=mrot(theta)
    Snew = m.dot(Sold)
    """
    s = np.sin(2.*theta)
    c = np.cos(2.*theta)

    m = np.array([[1.,0.,0.,0.],[0.,c,s,0.],[0.,-s,c,0.],[0.,0.,0.,1.]])

    return m


def sc1(lam,perfect):
    """ Steady state scattering from Earth's atmosphere
    SI units.   Use 10 millionths of disk centre
    """
    if perfect ==1:
        return np.zeros(4)
    else:
        i1=10.e-6*iplanck(lam,5900.) 
        p1=i1/100.  # small angle scattering
        x = 5000.e-10/lam
        y=np.array([i1,p1,0.,0.])
        m=mrot(1.)
        y = m.dot(y)
        # 4th power (rayleigh) and 1st power (particulates)
        sc1 = y*x*x*x*x + 0.1*y*x
        return sc1
    
def sc2(lam,index,perfect):
    """ 
    Scattering of solar light from Earth horizon light
    and then off the atmosphere. This is NOT FINISHED
    """
    return np.zeros(4)  # set to zero
    if perfect ==1:
        return np.zeros(4)
    else:
        ealbedo=0.28
        thetah = 1.
        sc2 = sc1(lam,perfect)*ealbedo
        m=mrot(thetah)
        sc2 = m.dot(sc2)
        return sc2

def  stel(lam):
    """
    telescope scattered light in brightness of QS
    lam := wavelength in metres
    just a complete guess
    """
    stel = [.1,.005,.001,0.]*1.e-6 * (5.e-7 / lam)
    return stel

def trapez(x,y):
    """
    trapez(x,y)

    performs trapezoidal integration.
    """
    n=len(x)
    integrand=(y[:n-1]+y[1:n])*(x[1:n]-x[:n-1])*0.5
    return sum(integrand)

def wlpspec(lam):
    """

    Approximate noise power spectrum as a function of frequency 
    measured by D. Elmore at MLSO
    normalized to 1
    lam := wavelength in m (req)
    output 
    f := frequency in Hz
    
    """
    #print ("Power spectrum in millionths of disk brightness")
    fmax=500.
    f=np.arange(0.01,fmax,0.05)
    
    p= 10**(-2.0-1.1*np.log10(f) ) + 3.e-4
    #
    # SI units use 
    #
    i0=1.e-6*iplanck(7000.e-10,5900.) 
    #    print("Millionth of disk = %0.2e " % i0 )
    p=p*i0*i0  #  power is intensity squared
    return [f,p]

#####################################################################################
#  CLASSES
#####################################################################################

class telescope:

    """ Defines telescope class

    name :=telescope name (required)
    diam :=telescope diameter metres (required)
    mm   :=telescope Mueller Matrix (required)
    syseff:= telescope system efficiency (required)
    scatt := scattering intensity (millionths of disk, required)

    """
    
    def __init__(self,name, diam,mm,syseff,scat):
        self.name = name
        self.diam = diam
        self.mm = mm
        self.syseff=  syseff
        self.scat=  scat

        description = 'Telescope parameters '
        print()
        print('------------------------------------------------------------')
        print ('TELESCOPE')
        print('------------------------------------------------------------')
        print ('  name     =     %s  ' % self.name)
        print ('  diameter = %0.1f m ' % self.diam)
        print ('  area     = %0.2e m2' % self.area())
        print ('  sys effic= %0.2e   ' % self.syseff)
        print ('  scat     = %0.2e   ' % self.scat)
        print ('  Mueller matrix:')
        print (self.mm)
        print('------------------------------------------------------------')
        print()
        
    def area(self):
            return area(self.diam/2.)


class spectrograph:

    def __init__(self,name, bin, binunit):
        self.name = name
        self.bin = bin
        self.binunit=binunit
        
        description = 'Spectrograph parameters '
        print()
        print('------------------------------------------------------------')
        print ('SPECTROGRAPH')
        print('------------------------------------------------------------')
        print ('  name =     %s ' % self.name)
        print ('  bin  = %0.0f   binunit  = %s ' % (self.bin,self.binunit) )
        print('------------------------------------------------------------')
        print()


class detector:

        def __init__(self,name, chip, readt, readn, qeff, sigx, sigy, \
                     pixel, pixsz,maxdn):

            self.name = name
            self.chip = chip
            self.readt = readt
            self.readn = readn
            self.qeff = qeff
            self.sigx=sigx
            self.sigy=sigy
            self.pixel=pixel
            self.pixsz=pixsz
            self.maxdn=maxdn
            self.omega=(pixel*const.pi/180./60./60.)**2.
                
            description = 'Detector parameters '
            print()
            print('------------------------------------------------------------')
            print ('DETECTOR')
            print('------------------------------------------------------------')
            print ('  name =     %s ' % self.name)
            print ('  chip  =    %s ' % self.chip)
            print ('  readt    = %s ' % self.readt)
            print ('  readn    = %s ' % self.readn)
            print ('  qeff     = %0.3f ' % self.qeff)
            print ('  sigx     = %0.3f' % self.sigx)
            print ('  sigy     = %0.3f' % self.sigy)
            print ('  pixel    = %0.3f arcsec' % self.pixel)
            print ('  pixsz    = %0.3e m     ' % self.pixsz)
            print ('  maxdn    = %0.3e ' % self.maxdn)
            print ('  omega    = %0.3e sterad' % self.omega)
            print('------------------------------------------------------------')
            print()
        

class analysis:

    def __init__(self,refd,iunit,lambda0,nstates,dt, \
                 perfect,dual):
        self.refd=refd
        self.iunit=iunit
        self.lambda0=lambda0
        self.nstates=nstates
        self.perfect=perfect
        self.dual=dual
        
        description = 'Analysis parameters '

        print()
        print('------------------------------------------------------------')
        print ('ANALYSIS')
        print('------------------------------------------------------------')
        print ('  refd =     %s ' % self.refd)
        print (' iunit =     %s ' % self.iunit)
        print ('  nstates=   %4i      ' % self.nstates)
        print ('  lambda0=   %0.1f AA' % self.lambda0)
        print ('  perfect=   %4i ' % self.perfect)
        print ('  dual=      %4i' % self.dual)
        print('------------------------------------------------------------')
        print()

####################################################################
