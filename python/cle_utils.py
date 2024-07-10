##collection of classes that cna be imported in scripts
## used to read and manipulate cle and database output.
## to import the custom ibrary in ipython
#sys.path.insert(1, '/home/alin/Documents/physics_prog/cle/python/')
#import cle_utils

## the following functions are included:

## out_params;   for reading the output cle stokes profiles/arrays
## atmos_params; for reading the cle atmosphere used to generate stokes profiles.
## get_db_entry; select the correct database entry corresponding to a radial height above the limb.
## readhdr;      reads the geometrical and phyics parameter spaces from a DB header file. 
## physa;        Returns the physics and geometry information for a selected (fitted) DB entry.
## atmos2vtk;    saves the nne te bx by bz contents as vtk files so they can be  rendered with a 3d viewer.
## scale_stokes  rescales all the quantities from both the db and observation to have all quantities between 0.1 and 1.

import numpy as np
import sys
import os.path
import math
import re
import struct
from pylab import *
from vtk.util import numpy_support
import vtk

#######################################################################
## OUT_PARAMS
## Created by:    Alin Paraschiv arparaschiv@nso.edu
## Original IDL version: Philip Judge
## For output from cle (OUT file) queries contact: judge@ucar.edu
##
## Change log
## 20190903 ARP ported the IDL version to python
## 20200929 ARP changed the terminal script to a class function that can be easily called in python
#######################################################################

## Read in the CLE output for python 


## Write the output map using a class
class out_params:
    file     = None    ## file:     file containing atmospheric parameters
    model    = None    ## model:    name of fortran routine used to determine coronal parameters
    x        = None    ## x,y       positions in solar radii',
    y        = None       
    xc       = None    ## position of central pixels in range
    yc       = None
    dx       = None    ## pixel sizes
    dy       = None
    punit    = None    ## position of limb in arcsec or solar radii
    unit     = None    ## string for Punit unit.
    wave     = None    ## Central wavelength of data
    dunits   = None    ## Phyical units
    data_q   = None    ## Wavelength Array to store full spectroscopy of profiles
    data_full= None    ## Data Array to store the full spectroscopy of profiles
    data1    = None    ## Stokes I array
    data2    = None    ## Stokes I array
    data3    = None    ## Stokes Q array
    data4    = None    ## Stokes U array
    data5    = None    ## Stokes V(magnetograph) array
    data1_id = None    ## String IDs for all 5 (IQUVV) stokes profiles calculated
    data2_id = None
    data3_id = None
    data4_id = None
    data5_id = None
    
    ## define a function that computes the spectroscopic properties from one pixel.
    def __init__(self,simpath_out,kr,no_rot):       ### filepath, line index, rotation from cle frame

        ## check if input arggument string containing path exists. 
        ## We need to define nline so that we know the range for kr

        if isinstance(simpath_out, str) and os.path.exists(simpath_out):
            file = simpath_out
        else:
            raise ValueError('OUT file path not found!')
        
        with open(file) as f:                   # Reads to "1 STOKES DATA OUTPUT"
            for line in f:
                if 'OUTPUTTED' in re.findall(r'\w+', line): 
                    nline=int(re.findall(r'\w+\.?\w*', line)[re.findall(r'\w+\.?\w*', line).index('OUTPUTTED')+1] )
                    break
                    
        if no_rot in [0,1] and kr in range(0,nline):   
            kr = kr
            no_rot =no_rot
        else:
            raise ValueError('input arguments not correct type, value or not in right order.')
        ##-------------------------------------------------------------------------------------------------

        with open(file) as f:                   # Reads to "1 INPUT PARAMETERS"
            for line in f:
                if 'QNORM' in re.findall(r'\w+', line): 
                    qnorm=float(re.findall(r'\w+\.?\w*', line)[re.findall(r'\w+\.?\w*', line).index('QNORM')+1] )
                    break
            for line in f:
                if 'IWLINE' in re.findall(r'\w+', line):
                    iwline=float(re.findall(r'\w+\.?\w*', line)[re.findall(r'\w+\.?\w*', line).index('IWLINE')+1] )
                    break
            for line in f:
                if 'CRTN' in re.findall(r'\w+', line):
                    crtn=re.findall(r'\w+\.?\w*', line)[re.findall(r'\w+\.?\w*', line).index('CRTN')+1] 
                    break

        with open(file) as f:                   # Reads to "1 GRID"
            for line in f:
                if 'X (LINE OF SIGHT)   :' in re.findall(r'\w.+:', line):
                    gxmin=float(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('SIGHT')+1] )
                    gxmax=float(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('SIGHT')+2] )
                    ngx  =int(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('SIGHT')+3] )
                if 'Y (PLANE OF SKY E-W):' in re.findall(r'\w.+:', line):
                    gymin=float(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('E-W')+1] )
                    gymax=float(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('E-W')+2] )
                    ngy  =int(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('E-W')+3] )
                if 'Z (PLANE OF SKY N-S):' in re.findall(r'\w.+:', line):
                    gzmin=float(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('N-S')+1] )
                    gzmax=float(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('N-S')+2] )
                    ngz  =int(re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line)[re.findall(r'\-?\w+\.?\w*\-?\+?\w*', line).index('N-S')+3] )

        ## coordinate space and results array definitions
        gx = np.linspace(gxmin,gxmax,ngx)
        gy = np.linspace(gymin,gymax,ngy)
        gz = np.linspace(gzmin,gzmax,ngz)

        ## Add extra information if a specific type of output is found
        ## (I need an example to fully write this part)
        #if crtn == 'FFL':
        #    cornotes = ''
        #    with open('OUT') as f:
        #    for line in f:
        #        if '1 CORONAL ROUTINE' in re.findall(r'\w+', line):

        alamb = np.zeros((nline,), dtype=np.float)
        nq =    np.zeros((nline,), dtype=np.int)
        q =     np.zeros((nline,100), dtype=np.float)

        with open(file) as f:                   # Reads to "1 STOKES DATA OUTPUT"
            i=0
            for line in f:
                if 'DELTA-WAVELENGTHS' in re.findall(r'\w+\-?\w+', line): 
                    nq[i]=float(re.findall(r'\w+\.?\w*', line)[-1] )
                    alamb[i]=float(re.findall(r'\w+\.?\w*', prevline)[-1] )
                    for j in range(0,nq[i]//10):
                        q[i,j*10:(j+1)*10]=[float(k) for k in re.findall(r'\-?\w+\.?\w+',next(f))]
                    i+=1
                prevline=line


        ## prepare to read stokes data

        emerg = np.zeros((nline,5,ngy,ngz), dtype=np.float, order='F')

        if iwline > 0: 
            full =  np.zeros((nline,5,max(nq),ngy,ngz), dtype=np.float)
            with open(file) as f:          # Reads to "1 STOKES DATA OUTPUT - FULL"
                for line in f:
                    if 'FULL EMERGENT STOKES PROFILES' in re.findall(r'\w.+', line): 
                        for i in range(0,ngy):
                            for j in range(0,ngz):
                                for k in range(0,nline):
                                    for l in range(0,6):
                                        if l==5:
                                            emerg[k,:,i,j]=[float(m) for m in re.findall(r'\-?\w+\.?\w+\+?\-?\w*',next(f))]
                                            if len(sys.argv) < 4 or no_rot == 0: 
                                                if emerg[k,1,i,j] != 0 or emerg[k,2,i,j] != 0:
                                                    alpha = 0.5*math.atan2(emerg[k,1,i,j],emerg[k,2,i,j])+math.pi/2
                                                    p = math.sqrt(emerg[k,1,i,j]**2+emerg[k,2,i,j]**2)
                                                    emerg[k,1,i,j]=p*math.cos(alpha)
                                                    emerg[k,2,i,j]=p*math.sin(alpha)
                                        else:
                                            full[k,l,0:nq[k],i,j]=[float(m) for m in re.findall(r'\-?\w+\.?\w+\+?\-?\w*',next(f))]
        else:
            with open(file) as f:                   # Reads to "1 STOKES DATA OUTPUT - FREQUENCY"
                for line in f:
                    if 'FREQUENCY-INTEGRATED' in re.findall(r'\w+\-?\w*', line):
                        for i in range(0,ngy):
                            for j in range(0,ngz):
                                for k in range(0,nline): 
                                    emerg[k,:,i,j]=[float(m) for m in re.findall(r'\-?\w+\.?\w+\+?\-?\w*',next(f))]
                                    if no_rot == 0: 
                                        if emerg[k,1,i,j] != 0 or emerg[k,2,i,j] != 0:
                                            alpha = 0.5*math.atan2(emerg[k,1,i,j],emerg[k,2,i,j])+math.pi/2
                                            p = math.sqrt(emerg[k,1,i,j]**2+emerg[k,2,i,j]**2)
                                            emerg[k,1,i,j]=p*math.cos(alpha)
                                            emerg[k,2,i,j]=p*math.sin(alpha)

        ## prepare map putputs
        print(type(kr))
        si    = emerg[kr,0,:,:].reshape(emerg.shape[2],emerg.shape[3])
        px    = emerg[kr,1,:,:].reshape(emerg.shape[2],emerg.shape[3])
        py    = emerg[kr,2,:,:].reshape(emerg.shape[2],emerg.shape[3])
        v     = emerg[kr,3,:,:].reshape(emerg.shape[2],emerg.shape[3])
        vmag  = emerg[kr,4,:,:].reshape(emerg.shape[2],emerg.shape[3])
        p     = np.sqrt(px**2+py**2)
        
        ## simple test of computing blos
        #cc=2.99792458e10 #cm s-1 
        #me= 9.1095e-28   #g
        #ee=4.8032e-10    #esu  g1/2 cm3/2 s−1
        #kb=1.38064852e-16  #erg/K
        #te= 10**6.22
        #mion=55.845*1.672621E-24  #  g    for FeXIII 
        #omegat=math.sqrt(2*kb*te/mion)     #thermal velocity without microturbulence
        #dlmd_d=omegat*alamb[kr]*1e-8/cc               #Doppler width in frequency units!   1e-8 -->conversion A to cm
        #dlmd_d=omegat/(1e-8* alamb[kr])               #Doppler width in wavelength units!  1e-8 -->conversion A to cm
        #blos= -8*math.pi*me*cc*v*dlmd_d/(3*(si+p)*ee)   #eq 14 Plowman, 2014: LOS magnetic field |B|cos(theta); units are cgs 
        ##not finished here; will do separately

        ## output to maps

        ## units 
        ##can be either arcsec or solar radii; default is solar radii.
        punit = 955.                     #arcsec
        unit = 'arcseconds'
        punit = 1.                       #solar radii
        unit = 'solar radii'

        xc = (gymax+gymin)/2*punit 
        dx = (gymax-gymin)/ngy*punit 
        yc = (gzmax+gzmin)/2*punit 
        dy = (gzmax-gzmin)/ngz*punit 
        x  = (gymin + np.arange(ngy)*dx)*punit 
        y  = (gzmin + np.arange(ngz)*dy)*punit 

        if no_rot !=0:
            sx='Q (ref=z-axis)'
            sy='U (ref=z-axis)'
        else:
            sx='Px'
            sy='Py'


        
        self.x            = x
        self.xc           = xc
        self.dx           = dx
        self.y            = y
        self.yc           = yc
        self.dy           = dy
        self.punit        = punit
        self.unit         = unit
        self.wave         = alamb[kr]
        self.dunits       = 'erg/cm^2/s'
        self.data1        = si
        self.data1_id     = 'Stokes I'
        self.data2        = px
        self.data2_id     = 'Stokes '+sx
        self.data3        = py
        self.data3_id     = 'Stokes '+sy
        self.data4        = v
        self.data4_id     = 'Stokes V'
        self.data5        = vmag
        self.data5_id     = 'Stokes V (magnetograph formula)'
        if iwline > 0: 
            self.data_q   = q[kr,0:nq[kr]]
            self.data_full= full[kr,:,:,:,:].reshape(full.shape[1],full.shape[2],full.shape[3],full.shape[4])
            print('Full Stokes profiles were also outputted, stored in s.data_q, s.data_full')
        else:
            print('Wavelength integrated Stokes were written')

#######################################################################
## End
#######################################################################



#######################################################################
# ATMOS_PARAMS
# Created by:    Alin Paraschiv arparaschiv@nso.edu
# Original IDL version: Philip Judge
#
# Change log
# 20190827 ARP ported the IDL version to python
# 20200930 ARP changed the terminal script to a class function that can be easily called in python
#######################################################################
#
#
# Read in atmospheric model for python
#

## put the data in a nice class so it mimics IDL structures
class atmos_params:
    file  = None        # file:     file containing atmospheric parameters
    model = None        # model:    name of fortran routine used to determine coronal parameters
    x     = None        # x,y,z:    positions in solar radii',$
    y     = None        #           (x is along LOS, y is E-W, z N-S in observer''s frame)
    z     = None
    bx    = None        # bx,by,bz: x,y,z-components of magnetic field in G
    by    = None
    bz    = None
    nne   = None        # nne:      electron density in /cm3
    te    = None        # te:       electron temperature in K
    v     = None        # v:        LOS velocity in km/s (+ve is a red-shift)
    vt    = None        # vt:       "turbulent" velocity in km/s


    ## define a function that reads the ATMOS content.
    def __init__(self,simpath_atmos):       ### filepath

        ## check if input arggument string containing path exists
        if isinstance(simpath_atmos, str) and os.path.exists(simpath_atmos):
            file = simpath_atmos
        else:
            raise ValueError('input ATMOS file not found!') # if no input argument assumes atmos is in current dir

        ## read the atmoshpere header to get the data dimensions.
        with open(file, 'br') as f3:
            buffer = f3.read()
            print ("Length of buffer is %d" % len(buffer))
            data_h = struct.unpack_from("<6s1d2f1l1d2f1l1d2f1l",buffer, offset=4)
        ## The format of the file contains 1 empty space, 6 string characters followed by a series of three sets (for x, y, and z directions) 
        ## of 1 double (empty), two floats (g_min, g_max), and long (ng_) 
        ## Format example: "<1f6s1d2f1l1d2f1l1d2f1l1d3l1d7f + repeat 3l1d7f as many times as dimensions allow"

        cname=data_h[0]
        gxmin,gxmax,ngx = float('%.8f'%data_h[2]), float('%.8f'%data_h[3]),data_h[4]  ##truncated so there are no insignificant trailing digits
        gymin,gymax,ngy = float('%.7f'%data_h[6]), float('%.8f'%data_h[7]),data_h[8]
        gzmin,gzmax,ngz = float('%.8f'%data_h[10]),float('%.8f'%data_h[11]),data_h[12]

        if ngx+ngy+ngz == 0:                    ##optional manual input for dimensions                  
           print ("dimensions not written correctly")
           ngx,ngy,ngz=map(float,input("Enter NGX, NGY, NGZ:").split() )

        ## coordinate space and results array definitions
        x = np.linspace(gxmin,gxmax,ngx)
        y = np.linspace(gymin,gymax,ngy)
        z = np.linspace(gzmin,gzmax,ngz)
        bx,by,bz,nne,te,v,vt = [np.zeros((int(ngx),int(ngy),int(ngz)), dtype=float,) for _ in range(7)]

        ## Read the data 
        ## The data starts at byte position 78 after header information ends
        ## Configuration is 3 longs (array positions) followed by 7 floats (bx, by,bz,nne,te,v,vt) and 1 doulbe space
        with open(file, 'br') as f3:
            buffer = f3.read()
            for i in range(78,len(buffer),56):       
                data_d = struct.unpack_from("<3l1d7f", buffer, offset=i)
                bx [data_d[0]-1,data_d[1]-1,data_d[2]-1] = data_d[4]
                by [data_d[0]-1,data_d[1]-1,data_d[2]-1] = data_d[5]
                bz [data_d[0]-1,data_d[1]-1,data_d[2]-1] = data_d[6]
                nne[data_d[0]-1,data_d[1]-1,data_d[2]-1] = data_d[7]
                te [data_d[0]-1,data_d[1]-1,data_d[2]-1] = data_d[8]
                v  [data_d[0]-1,data_d[1]-1,data_d[2]-1] = data_d[9]
                vt [data_d[0]-1,data_d[1]-1,data_d[2]-1] = data_d[10]

        self.x=x
        self.y=y
        self.z=z
        self.bx=bx
        self.by=by
        self.bz=bz
        self.nne=nne
        self.te=te
        self.v=v
        self.vt=vt
        
#######################################################################
## End
#######################################################################



# Small script for selecting the correct database entry corresponding to a radial height above the limb.
def get_db_entry(yobs,dbdir):
    
    ## get the distance in the db format
    iy=round((yobs-1.0)*1000)
    so='%0*d' % (4, iy)
    fil1=int(so)
    ##get the database entries
    dblist=os.listdir(dbdir)  
    norm=np.zeros(len(dblist))
    for i in range(0,len(dblist)):
        ## if file in list is not database generated, ignore it from the norm.
        if dblist[i].find(".DAT") == -1:
            norm[i]=1000    
        elif dblist[i].find("DB") == -1:
            norm[i]=1000 
        else:
            ## compute the norm between the observed distance and db file distance
            norm[i]=np.absolute(fil1- int(dblist[i].split('DB')[1].split(".DAT")[0]))     
        
    ##select the closest match to the inputed radius
    db_file =dblist[np.where(norm == norm.min())[0][0]]    
    
    return db_file
#######################################################################
## End
#######################################################################



##reads the geometrical and phyics parameter spaces from a DB header file. 
def readhdr(file):
#
    g=np.fromfile(file,dtype=float,sep=' ')
    ned=int(g[0]) 
    ngx=int(g[1])
    nbphi=int(g[2])
    nbtheta=int(g[3])
    xed=np.empty(ned)
    for k in np.arange(0,ned): xed[k]=g[4+k]
    gxmin=g[ned+4]
    gxmax=g[ned+5]
    bphimin=g[ned+6]
    bphimax=g[ned+7]
    bthetamin=g[ned+8]
    bthetamax=g[ned+9]
    bfield=g[ned+10]
    nline=int(g[ned+11])
    nq=np.ones(nline,dtype=int32)
    return [g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
            bthetamin, bthetamax, bfield, nline, nq]
#######################################################################
## End
#######################################################################

# 2.0.1 NEEDED FUNCTIONS BEGIN
def readhdr201(file):
    #
    # just read the header of the DB*DAT data cubes
    #
    g=np.fromfile(file,dtype=np.float,sep=' ')
    #return g

    #def parsehdr(g):
    #
    # parse the header of the DB*DAT data cubes
    #
    # two kinds of data are returned
    # 1. the linear coefficents of the form min, max, n
    #    - x,theta,phi
    # 2. logarithmically spaced parameters
    #    - xed  (electron density array)
    #    - bfield (magnetic field strength array)
    #   for these also the min and max and number are also returned.
    #  
    ned=int(g[0]) 
    ngx=int(g[1])
    nbphi=int(g[2])
    nbtheta=int(g[3])
    nb=int(g[4])
    # log of Ne
    emin=g[5]
    emax=g[6]
    xed=lgrid(emin,emax,ned)
    # log of |B|
    bmin=g[13]
    bmax=g[14]
    xb=lgrid(bmin,bmax,nb)
    #
    gxmin=g[7]
    gxmax=g[8]
    bphimin=g[9]
    bphimax=g[10]
    bthetamin=g[11]
    bthetamax=g[12]
    #
    nline=int(g[15])
    wavel=np.empty(nline)
    for k in np.arange(0,nline): wavel[k]=g[16+k]    
    return [g, ned, ngx, nbphi, nbtheta, nb, xed, gxmin,gxmax, bphimin, bphimax,\
            bthetamin, bthetamax, xb, nline, wavel]

def lgrid(mn,mx,n):
    #
    if n == 1:
        ff=10.** (mx-mn)
    else:    
        ff=10.** ((mx-mn)/(n-1))
    x=np.empty(n)
    x[0]=10.**mn
    for k in np.arange(1,n): x[k]=x[k-1]*ff
    return x
#######################################################################
## End
#######################################################################

##Returns the physics and geometry information for a selected (fitted) DB entry.
def physa(index,gy,grd):
    i,j,k,l = params(index,grd)
    #
    ned=int(grd[0])
    xed=np.empty(ned)
    for kl in np.arange(0,ned): 
        xed[kl]=grd[5+kl]
    #
    n=4+ned
    gx = grd[n] + (grd[n+1]-grd[n])*j/grd[1]
    ed = xed[i]* electron_density(np.sqrt(gy*gy+gx*gx))
    #
    n+=2
    bphi  = grd[n] + (grd[n+1]-grd[n])*k/grd[2]
    n+=2
    btheta= grd[n] + (grd[n+1]-grd[n])*l/grd[3]
    b=btheta*0.+1.
    out=np.array([np.log10(ed),gy,gx,bphi*180./np.pi,btheta*180./np.pi,b])
    
    return out

    ## index helper for physa
def params(index,grd):
    # for index, get i,j,k,l   indices
    ned = int(grd[0])
    ngx = int(grd[1])
    nbphi = int(grd[2])
    nbtheta = int(grd[3])
    n4=(ngx*nbphi*nbtheta)
    n3=(        nbphi*nbtheta)
    n2=(              nbtheta)
    #
    i = index   //        n4

    j = index   -     i * n4
    j = j      //         n3

    k = index   -     i * n4 - j * n3
    k = k      //         n2

    l=  index   -     i * n4 - j * n3 - k * n2

    return np.array([i,j,k,l])


    ##Gets a estimate of the local density based on the radial distance from the limb used in physa
    ##Baumbach formulation; see Allen, 1973 
def electron_density(r):
    baumbach = 1.e8*(0.036/r**1.5 + 1.55/r**6.)
    hscale=   7.18401074e-02  # 50 Mm
    n0=3.e8*np.exp(- (r-1.)/hscale) + baumbach
    return n0
    

#######################################################################
## End
#######################################################################


## write the ATMOS file with vtk for paraview viewing
## writes the te, nne, v,vt,bb,by,bz 3d arrays in vtk format.
## Call as: atmos2vtk(atmos_file,path_to_save,file_name_start(gets autocompleted with quantity),te=0,nne=0,bx=0,by=0,bz=0,v=0,vt=0 (set which atmos variables to write as individual files))

def atmos2vtk(atm,pathn,filen,te=0,nne=0,bx=0,by=0,bz=0,v=0,vt=0):
    
    if te != 0:
        data = atm.te
        print(data.shape)
        # vtkImageData is the vtk image volume type
        imdata = vtk.vtkImageData()
        # this is where the conversion happens
        depthArray = numpy_support.numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        # fill the vtk image data object
        imdata.SetDimensions(data.shape)
        imdata.SetSpacing([1,1,1])
        imdata.SetOrigin([0,0,0])
        imdata.GetPointData().SetScalars(depthArray)
        # f.ex. save it as mhd file
        writer = vtk.vtkMetaImageWriter()
        writer.SetFileName(pathn+filen+"_te.mhd")
        writer.SetInputData(imdata)
        writer.Write()
    #####------------------------------
    
    if nne != 0:
        data = atm.nne
        print(data.shape)
        # vtkImageData is the vtk image volume type
        imdata = vtk.vtkImageData()
        # this is where the conversion happens
        depthArray = numpy_support.numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        # fill the vtk image data object
        imdata.SetDimensions(data.shape)
        imdata.SetSpacing([1,1,1])
        imdata.SetOrigin([0,0,0])
        imdata.GetPointData().SetScalars(depthArray)
        # f.ex. save it as mhd file
        writer = vtk.vtkMetaImageWriter()
        writer.SetFileName(pathn+filen+"_nne.mhd")
        writer.SetInputData(imdata)
        writer.Write()
    #####------------------------------
    
    if bx != 0:
        data = atm.bx
        print(data.shape)
        # vtkImageData is the vtk image volume type
        imdata = vtk.vtkImageData()
        # this is where the conversion happens
        depthArray = numpy_support.numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        # fill the vtk image data object
        imdata.SetDimensions(data.shape)
        imdata.SetSpacing([1,1,1])
        imdata.SetOrigin([0,0,0])
        imdata.GetPointData().SetScalars(depthArray)
        # f.ex. save it as mhd file
        writer = vtk.vtkMetaImageWriter()
        writer.SetFileName(pathn+filen+"_bx.mhd")
        writer.SetInputData(imdata)
        writer.Write()
    #####------------------------------

    if by != 0:
        data = atm.by
        print(data.shape)
        # vtkImageData is the vtk image volume type
        imdata = vtk.vtkImageData()
        # this is where the conversion happens
        depthArray = numpy_support.numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        # fill the vtk image data object
        imdata.SetDimensions(data.shape)
        imdata.SetSpacing([1,1,1])
        imdata.SetOrigin([0,0,0])
        imdata.GetPointData().SetScalars(depthArray)
        # f.ex. save it as mhd file
        writer = vtk.vtkMetaImageWriter()
        writer.SetFileName(pathn+filen+"_by.mhd")
        writer.SetInputData(imdata)
        writer.Write()
    #####------------------------------
    
    if bz != 0:    
        data = atm.bz
        print(data.shape)
        # vtkImageData is the vtk image volume type
        imdata = vtk.vtkImageData()
        # this is where the conversion happens
        depthArray = numpy_support.numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        # fill the vtk image data object
        imdata.SetDimensions(data.shape)
        imdata.SetSpacing([1,1,1])
        imdata.SetOrigin([0,0,0])
        imdata.GetPointData().SetScalars(depthArray)
        # f.ex. save it as mhd file
        writer = vtk.vtkMetaImageWriter()
        writer.SetFileName(pathn+filen+"_bz.mhd")
        writer.SetInputData(imdata)
        writer.Write()   
    #####------------------------------

    if v != 0:     
        data = atm.v
        print(data.shape)
        # vtkImageData is the vtk image volume type
        imdata = vtk.vtkImageData()
        # this is where the conversion happens
        depthArray = numpy_support.numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        # fill the vtk image data object
        imdata.SetDimensions(data.shape)
        imdata.SetSpacing([1,1,1])
        imdata.SetOrigin([0,0,0])
        imdata.GetPointData().SetScalars(depthArray)
        # f.ex. save it as mhd file
        writer = vtk.vtkMetaImageWriter()
        writer.SetFileName(pathn+filen+"_v.mhd")
        writer.SetInputData(imdata)
        writer.Write()    
    #####------------------------------
    if vt != 0:    
        data = atm.vt
        print(data.shape)
        # vtkImageData is the vtk image volume type
        imdata = vtk.vtkImageData()
        # this is where the conversion happens
        depthArray = numpy_support.numpy_to_vtk(data.ravel(order='F'), deep=True, array_type=vtk.VTK_DOUBLE)
        # fill the vtk image data object
        imdata.SetDimensions(data.shape)
        imdata.SetSpacing([1,1,1])
        imdata.SetOrigin([0,0,0])
        imdata.GetPointData().SetScalars(depthArray)
        # f.ex. save it as mhd file
        writer = vtk.vtkMetaImageWriter()
        writer.SetFileName(pathn+filen+"_vt.mhd")
        writer.SetInputData(imdata)
        writer.Write()
    #####------------------------------
    
    return "ATMOS file converted to VTK: "+pathn
# #######################################################################
## End
#######################################################################


##  Re-scale q, u and v by fq and fv to make s roughly (1,0-1,0-1,0-1)
##  accordingly, the noise magnitudes must be scaled by the same factors:
##  Note: this tends to change only the chi2 values by a single multiplier
##  thereby not affecting the ordering of the chi2
##  Note 2: First note is not valid as this would not give cvasi equal weighting to the chi^2 of each components. 
##  e.g. v~10^-5-> very small chi^2 to be added for this component
##  This leads to fitting only I to become dominant in the sorting.

##input is the observed array, the db array, the number of lines and the counts used to compute the noise influence.
def scale_stokes(sobs_in,sdb_in,nuse,counts):
    
    ## scale the observation such that all 8 components are between -1 and 1
    mx1=np.max(sobs_in)
    sobs_temp = (sobs_in/mx1)
    plc = np.where(sobs_temp == 1)[0]   ## so we know which element to use as a pivot.

    # use the same scaling factor for both observations and DB entries
    if nuse == 1 or nuse == 2:
        f = np.zeros(4*nuse)
        for ii in range(0,4*nuse):
            if ii != plc:
                f[ii]=10.**np.fix(np.log10(np.abs(sobs_temp[plc]/sobs_temp[ii])))/sobs_temp[plc] #stokes*10^order of magnitude from bigest value/
            else:
                f[ii]=1
    else:
        print("Not a correct line number; Only 1 or 2 observation combinations are accepted.")
        
        
#     if nuse == 2: 
#         fq1=10.
#         fu1=10.
#         fv1=100.
#         fq2=10.
#         fu2=10.
#         fv2=100.
        
#         f1=np.array([1,fq1,fu1,fv1])
#         f2=np.array([1,fq2,fu2,fv2])
#         f=np.append(f1,f2)
    
#     elif nuse == 1:
#         fq1 = sobs_temp[1]*10**np.fix(np.log10(np.abs(sobs[0]/sobs[1])))/sobs[0]
#         fu1 =
#         fv1 = 
#         f=np.array([1,fq1,fu1,fv1])
    
#     else:
#         print("Not a correct line number; Only 1 or 2 observation combinations are accepted.")



#     # First, here is variance of each state with counts electrons
#     #
#     variance = counts
#     #
#     sobs_temp*=variance           # number of counts in each sobs

#     # Here is variance vars in ((I+S) - (I-S))/2 = S, in counts
#     #
#     vars = (variance +variance)/2.
#     #
#     # here is the variance in y = S/I:  var(y)/y^2 = var(S)/S^2 + var(I)/I^2
#     #
#     vary =  vars/sobs_temp**2 + variance/sobs_temp[plc]**2
#     vary *= (sobs_temp/sobs_temp[plc])**2
#     #
#     # Now add noise to the fake observations
#     # one realization 
#     #sobs *= (1. +rms*[0. , -1.73,  0.05, -0.13,  1.44,  1.29, -1.11, -1.02])
#     sobs_temp/=variance
#     sobs_temp[plc]=1.
#     rms=sqrt(vary)

    # Now fix the rms values for the denominator in the chi2 calculation
    # this is for each stokes value
    rms=f/np.sqrt(counts)     

    sobs_out=sobs_temp * f*( 1. + rms*np.random.normal(0.,1,f.shape) )#*f

    # apply form factor f to database as well as obs
    mx2=np.max(sdb_in,axis=1)
    sdb_out = (sdb_in.T/mx2.T).T*f # apply form factor f to database as well as obs


    
    return [sobs_out,sdb_out,rms]


# #######################################################################
## End
#######################################################################




