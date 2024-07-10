#
# Needed libraries
#
import numpy as np
from scipy.io import FortranFile
from pylab import *
from matplotlib import pyplot as plt
from scipy import stats
import time
import glob
from random import seed
from random import randint
import sys
import os
import numexpr as ne
import timeit
import pdb
import gc # garbage collection
from numba import jit,njit, prange
from numba.typed import List  ## standard reflected python lists will be deprecated in numba

#global database
#global r2d,nsearch,npartition,reduction, elongation

np.set_printoptions(precision=4,suppress=False)

#
######################################################################
#
def dcompress(i):
# Constants here must correspond to those in dbe.f in the CLE main directory
# 32767 is the limit of 4-byte ints
# 15 is the number of orders of magnitude for the range of intensities
#
    verbose=0
    cnst=-2.302585092994046*15./32767.
    start_dc=time.time()
    strt=time.time()
    if verbose ==1:print(np.int(sys.getsizeof(i)/1.e6) , 'MB in database file')
    negv=np.flatnonzero(i  < 0 )
    dt= "{:4.3f}".format(time.time()-strt)
    if verbose ==1:print('DCOMPRESS:', dt,' SECONDS FOR WHERE')
    
    strt=time.time()
    #f=np.abs(i.astype(float))
    f=np.abs(i)*cnst

    #t = timeit.Timer(lambda: np.exp(f))
    #print('--------------------TIMER np.exp ',t.timeit(1))

    f=ne.evaluate("exp(f)")
    #import pdb; pdb.set_trace()

    dt= "{:4.3f}".format(time.time()-strt)
    if verbose ==1:print('DCOMPRESS:', dt,' SECONDS FOR EXP')

    dt= "{:4.3f}".format(time.time()-start_dc)
    if size(negv) > 0:f[negv]=-f[negv]
    if verbose ==1:print('DCOMPRESS:', dt,' SECONDS TOTAL')
    del negv
    return f
#
######################################################################
#
def electron_density(r):
  baumbach = 1.e8*(0.036/r**1.5 + 1.55/r**6.)
  hscale=   7.18401074e-02  # 50 Mm
  n0=3.e8*np.exp(- (r-1.)/hscale) + baumbach
  return n0

#
######################################################################
#
def phys(index,gy,grd,b):
  #
  i,j,k,l = params(index,grd)
  phs=physa(index,gy,grd,b)
  bphi=phs[3]
  btheta=phs[4]
  bx=b*np.sin(btheta)*np.cos(bphi)
  by=b*np.sin(btheta)*np.sin(bphi)
  bz=b*np.cos(btheta)
  out=np.array([phs[0],phs[1],phs[2],bx,by,bz])
  return out
#
######################################################################
#
def physa(index,gy,grd,b):
  #
  i,j,k,l = params(index,grd)
  #
  grd,ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax,\
    bthetamin, bthetamax, nline,wavel  = parsehdr(grd)

  gx = gxmin + j*(gxmax-gxmin)/(ngx-1)
  #
  bphi  = bphimin + k*(bphimax-bphimin)/(nbphi-1)
  btheta  = bthetamin + l*(bthetamax-bthetamin)/(nbtheta-1)
  #
  # log of Ne
  #
  ed = xed[i]* electron_density(np.sqrt(gy*gy+gx*gx))
  out=np.array([np.log10(ed),gy,gx,bphi,btheta,b])
  return out
#
######################################################################
#
def params(index,grd):
  #
  # for index, get i,j,k,l   indices
  #
  ned = int(grd[0])
  ngx = int(grd[1])
  nbphi = int(grd[2])
  nbtheta = int(grd[3])
  
  n5=( ngx*nbphi*nbtheta)
  n4=(     nbphi*nbtheta)
  n3=(           nbtheta)
  #
  i = index   //        n5
  j = index   -     i * n5
  j = j      //         n4
  k = index   -     i * n5 - j * n4
  k = k      //         n3
  l=  index   -     i * n5 - j * n4 - k * n3

  return np.array([i,j,k,l])
#
######################################################################
#
def invparams(i,j,k,l,grd):
  #
  # for i,j,k,l get index
  #
  ned = int(grd[0])
  ngx = int(grd[1])
  nbphi = int(grd[2])
  nbtheta = int(grd[3])
  return i*ngx*nbphi*nbtheta + j*nbphi*nbtheta + k*nbtheta + l 
#
######################################################################
#
def readhdr(file):
#
# just read the header of the DB*DAT data cubes
#
  g=np.fromfile(file,dtype=np.float32,sep=' ')
  return g

#
######################################################################
#
def parsehdr(g):
#
# parse the header of the DB*DAT data cubes
#
# two kinds of data are returned
# 1. the linear coefficents of the form min, max, n
#    - x,theta,phi
# 2. logarithmically spaced parameters
#    - xed  (electron density array)
#   for these also the min and max and number are also returned.
#  
  ned=int(g[0]) 
  ngx=int(g[1])
  nbphi=int(g[2])
  nbtheta=int(g[3])
  #
  # log of Ne
  #
  emin=g[4]
  emax=g[5]
  xed=lcgrid(emin,emax,ned)
  gxmin=g[6]
  gxmax=g[7]
  bphimin=g[8]
  bphimax=g[9]
  bthetamin=g[10]
  bthetamax=g[11]
  nline=int(g[12])
  wavel=np.empty(nline)
  for k in np.arange(0,nline): wavel[k]=g[13+k]    
  return [g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
          bthetamin, bthetamax, nline, wavel]
#
######################################################################
#
def lcgrid(mn,mx,n):
  #
  ff=10.** ((mx-mn)/(n-1))
  xf=np.empty(n)
  xf[0]=10.**mn
  for k in np.arange(1,n): xf[k]=xf[k-1]*ff
  return xf
#
######################################################################
#
def phys_axnames(value, tick_number):
    # use:
    #  ax.xaxis.set_major_formatter(plt.FuncFormatter(phys_axnames))
    # 
    if value == 0:
        return "$\log_{10} N_e}$"
    elif value == 1:
        return r"$y_0$"
    elif value == 2:
        return r"$x$"
    elif value == 3:
        return r"$b_x$"
    elif value == 4:
        return r"$b_y$"
    elif value == 5:
        return r"$b_z$"
    else:
        return ""
#
######################################################################
#
def find_elongation(dbdir,y):
#
# returns filename and index of elongation y in database
#
  names=glob.glob(dbdir+"DB*.DAT")
  nm=len(names)
  numbers=np.zeros(nm,dtype=np.int32)
  for i in range(0,nm):
    nn=str.find(names[i],'DB')
    string=names[i]
    string=string[nn+2:nn+6]
    numbers[i]=int(string)
  #
  #Use naming convention
  #
  y_of_files= 1.+ numbers/1000.
  yi=np.argmin(np.abs(y_of_files-y))
  file=names[yi]
  iy=numbers[yi]  # set it to the number ID of the file
  strc="{:0.3f}".format(y_of_files[yi]) 
  #print("Elongation =",y,"nearest DB file is",file)
  return file,iy

      
def obs_fake(dbdir,yobs,oindex,counts,**keys):
    #
    # INPUTS  yobs, counts
    # 
    # Normally one would read from an observation file.
    # Instead in this case read from files cgridobs, obs.dat, obs.hdr
    # Reads files 'obs.hdr' and filen (e.g. OB0100.DAT) for observations
    # from filename, i.e. 0100/1000 + 1 solar radii = 1.1 solar radii
    # The idea is to return data along a certain elongation y from Sun
    # 
    #
    # Set all the keywords to values
    #
    nkeys=len(keys)
    values=["" for x in range(nkeys)]
    kwords = ["" for x in range(nkeys)]
    kount=0
    for arg in keys.values():
      values[kount] = arg
      kount+=1
    kount=0
    verbose=0
    if verbose ==1:print("Keywords...")
    for arg in keys:
      kwords[kount] = arg
      if kwords[kount] == 'dx': dx=values[kount]
      if verbose ==1:print(kwords[kount],'=',values[kount],end=" ")
      kount+=1
    if verbose ==1:print("\n")
    #
    #  end of keyword decoding
    #
    # The data returned are just
    #  sobs    -    an array of size [4*nuse, nobs], normalized to maximum = 1.0
    #  rms     -    an array of size (4*nuse]  uncertainties
    #               nuse = number of lines to use = nline by default
    #
    # header file

    g=readhdr(dbdir+'db.hdr')
    g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax, \
        bthetamin, bthetamax, nline, wavel  = parsehdr(g)
    cgrid=g
    #
    #  set number of lines to use 
    #
    nuse=nline
    #
    #  get the elongation file name
    #
    start=time.time()
    fildb,iy=find_elongation(dbdir,yobs)
    filn = open(fildb, 'rb')
    #
    print('FAKE_OBS: READING DATABASE FILE & DECOMPRESSING DATA...')
    database=dcompress( np.fromfile(file=filn, dtype=np.int16) )
    filn.close()
    # reshape on the fly
    i,j,k,l= params(oindex,g)
    shap=[ned,ngx,nbphi,nbtheta,nline,4]
    tmp=np.reshape(database,shap)

# with option dx
    if nkeys > 0:
        j=9   # select something small to begin addition
        print("dx=",dx)
    s0=tmp[i,j,k,l,0,:4]
    s1=tmp[i,j,k,l,1,:4]
    if nkeys > 0:
        #print('sobs bef',s0)
        s0+=.1*tmp[i,j+dx,k,l,0,:4]
        #print('sobs aft',s0/2.)
        s1+=.1*tmp[i,j+dx,k,l,1,:4]
    oindex=invparams(i,j,k,l,g)
# end of option    
    sobs=np.append(s0,s1)  
    sobs/=sobs[0]          # normalized
    sobs*=counts           # number of counts in each sobs
    del s0,s1,tmp
    # First, here is variance of each state with counts electrons
    #
    variance = counts
    #
    # Here is variance vars in ((I+S) - (I-S))/2 = S, in counts
    #
    vars = (variance +variance)/2.
    #
    # here is the variance in y = S/I:  var(y)/y^2 = var(S)/S^2 + var(I)/I^2
    #
    var =  vars/sobs**2 + variance/sobs[0]**2
    var *= (sobs/sobs[0])**2
    rms=np.sqrt(var)
    #
    # Now add noise to the fake observations
    #
    sobs = sobs * ( 1. + rms*np.random.normal(0.,1,sobs.shape) )
    sobs/=counts
    sobs[0]=1. 
    #
    #  |B|
    #
    r=-0.5+2*np.random.uniform(0.,1)
    b=10.**r
    print('Observed B is ',b, 'G , r=',r)
    sobs[[3,7]] *= b
    rms[[3,7]] *= np.sqrt(b)
    ######################################################################
    # finished observations, which consist only of
    #
    #  nuse, sobs, rms
    #
    # yobs is the elongation in Rsun of the observed y-position 
    # index is the index of the calculation used to make fake obs
    #       (not needed for real observations)
    #
    ######################################################################
    #print("\n")
    #strc="{:0.1f}".format(np.log10(counts))
    #print('Fake observations: log counts per state =',strc)
    #print('S     ', [ "{:9.2e}".format(x) for x in sobs ])
    #print('rms/S ',[ "{:9.2e}".format(x) for x in rms/sobs ])
    return [nuse,sobs,rms,oindex]


#############################################################################  
#############################################################################  
#
#  MAIN SEARCH ALGORITHM
#
#############################################################################  
#############################################################################  

#@njit(parallel=True)
def clematch(sobs,yobs,rms,**keys):
#
# keep these variables in memory, analogous to common fortran
#
    global database
    global r2d,nsearch,npartition,reduction, elongation
    #
    # check number of elements of sobs
    #
    numlines=int(len(sobs)/4)
    if numlines < 2:
      print('clematch: must be more than 1 observed line to continue')
      sys.exit()
    timeentry=time.time()
    #
    # Default values replaced with any keyword values that are set:
    dbdir='dbcle/test/' ; outdir='results/fe13/' ; maxchisq=1.4 ; iplot=0 ; verbose=0
    #
    # Set all the keywords to values
    #
    nkeys=len(keys)
    values=["" for x in range(nkeys)]
    kwords = ["" for x in range(nkeys)]
    kount=0
    for arg in keys.values():
      values[kount] = arg
      kount+=1
    kount=0
    if verbose ==1:print("Keywords...")
    for arg in keys:
      kwords[kount] = arg
      if kwords[kount] == 'dbdir': dbdir=values[kount]
      if kwords[kount] == 'outdir': outdir=values[kount]
      if kwords[kount] == 'maxchisq': maxchisq=float(values[kount])
      if kwords[kount] == 'iplot': iplot=int(values[kount])
      if kwords[kount] == 'verbose': verbose=int(values[kount])
      if verbose ==1:print(kwords[kount],'=',values[kount],end=" ")
      kount+=1
    if verbose ==1:print("\n")
    #
    #  end of keyword decoding
    #
    # Here are parameters needed but not keywords
    #
    nsearch=4
    try: nsearch
    except:
      nsearch=4             #  number of theta angles to keep in reduced search
      print('DEFAULT NSEARCH = ',nsearch)
    npartition=8         # initial guess as to number of solutions with chi^2 < maxchisq
    reduction=1     # select database using U/Q (very fast when nsearch << nbtheta)
    r2d=180./np.pi     #; r2d=1.
    #
    # all parameters are now set.
    #    
    #######################################################################
    # Now identify the DB file to be read defined by yobs.
    #
    #
    #  get the elongation file name
    #
    fildb,iy=find_elongation(dbdir,yobs)
    #
    # Next read cgrid parameters from the accompanying db.hdr file in database
    #
    start=time.time()
    dbcgrid=readhdr(dbdir+'db.hdr')
    cgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
        bthetamin, bthetamax, nline, wavel  = parsehdr(dbcgrid)
    if numlines != nline:
      print("Must be the same number of lines in observations and database")
      print('numlines=',numlines,'nline=',nline)
      sys.exit()
    #
    dthb=(bthetamax-bthetamin)/ nbtheta # required for using U/Q to reduce datatbase search
    #
    #  decision- read a new data file?
    #
    try:elongation
    except:elongation=0
    delta = iy-elongation
    exist_data = 1
    try:database
    except: exist_data=0
    #
    #  here reading data is done and stored in variable database of size (n,2*nline)
    #  First check to see if the elongation has changed
    #
    #print('Delta ',delta, 'exist_data',exist_data)
    if delta !=0 or exist_data == 0:
        start=time.time()
        if verbose ==1:print('CLEMATCH: READING DATABASE...')
        elongation=iy
        #
        filn = open(fildb, 'rb')
        #
        start=time.time()
        database=dcompress(np.fromfile(file=filn, dtype=np.int16))
        filn.close()
        shape=[ned*ngx*nbphi*nbtheta,nline,4]
        database=np.reshape(database,shape)
        dt= "{:4.3f}".format(time.time()-start)
        if verbose == 1:
          print(dt,' SECONDS FOR DCOMPRESS')
          print("DB file parameters\n",dbdir,"\n",outdir,"\n",yobs,"\n", fildb,"\n")
          dt= "{:4.3f}".format(time.time()-start)
          print(dt,' SECONDS FOR DB READ/RESHAPE')
          #
        start=time.time()
        #
        #  Ok so the data are read.        
        #  THIS IS THE START OF THE SEARCH ALGORITHM
        #  get the cgrid associated with the data
        #
        dt= "{:4.3f}".format(time.time()-start)
        if verbose ==1:print(dt,' SECONDS FOR DB READ')
        start=time.time()
        #
        sdb=np.append(database[:,0,:4],database[:,1,:4],axis=1)
        dt= "{:4.3f}".format(time.time()-start)
        if verbose ==1:print(dt,' SECONDS FOR SDB CALC (append)')
    #######################################################################
    #  call routine to deliver the reduced dataset
    #
    sdb,outindex=get_subset(sobs,dbcgrid)
    # ndof = number of degree of freedom in model: ne, x, bphi, btheta, b
    ndata = nline*4 -1
    ndof=5
    denom=nline*4-1 - ndof  # used below in chi^2
    if denom < 1:
      print("Number of observables is insufficient to determine model parameters")
      print("N(data points)=",ndata," N(degrees of freedom)=",ndof)
      sys.exit()
    # So far so good.. 
    #
    #######################################################################
    ######################################################################
    #
    #  Now find the matched data
    #
    ######################################################################
    #
    #  This is the most time consuming 
    #
    start=time.time()
    # 
    # We do not use Stokes V to define the geometry
    # Therefore we use a subset of Stokes parameters
    # for  QU (line 1) and IQU (line 2)
    #
    subset=[1,2,4,5,6]
    
    diff= (sdb[:,subset]-sobs[subset]) / rms[subset] 
    diff*=diff
    dt= "{:4.3f}".format(time.time()-start)
    #
    # 
    # 
    chisq = np.sum( diff, axis=1)/denom
    dt= "{:4.3f}".format(time.time()-start)
    if verbose ==1:print(dt,' SECONDS FOR CALC CHISQ')
    #
    start=time.time()
    asrt = np.argpartition(chisq, range(nsearch))[:nsearch]
    dt= "{:4.3f}".format(time.time()-start)
    if verbose ==1:print(dt,' SECONDS FOR SORT  CHISQ')
    ix=asrt[0]
    ixr=ix
    #
    #  return indices and chisq of nearest solutions
    #
    si=0
    ixr=asrt[si]
    chisqmin=chisq[ixr]
    while chisq[ixr] < maxchisq and si < nsearch:
      ixr=asrt[si]
      dd=chisq[ixr]
      if(dd < maxchisq):
        #
        # here if reduced we must get the original value of ix
        #
        if reduction  == 1:
          cgridnew=dbcgrid + 0.0
          cgridnew[3]=nsearch
          i,j,k,l=params(ixr,cgridnew)
          newl=np.int(outindex[k,l])
          ix=invparams(i,j,k,newl,dbcgrid) # = original value ix
        #
        # Magnetic field strength from V/I of first line
        # 
        bfield=sobs[3]/database[ixr,0,3]
        #
        if si == 0:
          index=[ix]
          sfound=sdb[ixr]  # use reduced data index
          chisquare=[dd]
          ph=[phys(ix,yobs,dbcgrid,bfield)]
        else:
          #
          index = np.append(index, [ix], 0)
          sfound=np.append(sfound,sdb[ixr])
          chisquare=np.append(chisquare,[dd])
          ph=np.append(ph,[phys(ix,yobs,dbcgrid,bfield)])
        #
      si+=1
    dt= "{:4.3f}".format(time.time()-timeentry)
    if verbose ==1:print(dt,' SECONDS FOR SOLUTION')
    try:index
    except: return ixr, [0.,0.,0.,0.,0.,0.,0.,0.],chisq[ixr],[0.,0.,0.,0.,0.,0.]
    #
    # reshape arrays so they are readable by a human
    #
    dims=len(index)
    sshape=[dims,nline*4]
    sfound=np.reshape(sfound,sshape)
    pshape=[dims,(ndof+1)]
    ph=np.reshape(ph,pshape)
    gc.collect()
    return index, sfound, chisquare, ph
#

######################################################################
def get_subset(sobs,dbcgrid):
  # 
  # returns a subset of the sdb array compatible with yobs
  #
  global database
  global r2d,nsearch,npartition,reduction, elongation
  verbose=0

  start=time.time()

  cgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    bthetamin, bthetamax,  nline, wavel  = parsehdr(dbcgrid)

  kk=np.arange(0,nbphi)
  bphir  = bphimin + kk*(bphimax-bphimin)/(nbphi-1)  # cle phi array

  ll=np.arange(0,nbtheta)
  bthetar=bthetamin + ll*(bthetamax-bthetamin)/(nbtheta-1) #cle theta array
  dthb=(bthetamax-bthetamin)/nbtheta
  tt = np.tan(bthetar)  #cle tan theta array

  shape=[ned,ngx,nbphi,nbtheta,nline,4]
  # use reshape on the fly to get datasel
  datasel=(np.reshape(database,shape) )[:,:,:,:nsearch,:,:]   # selection of subset
  outindex=np.zeros([nbphi,nsearch])
  #
  gnew=dbcgrid + 0.
  gnew[3]=nsearch
  #
  #   Below, the relation is used: tan Phi_B = sin phi * tan theta
  #   here is tan Phi_B:
  #
  phib_obs=-0.5*np.arctan2( sobs[2],sobs[1] )
  tphib_obs= np.tan(phib_obs)
  #
  # Find those indices compatible with phib observed
  #
  shapeseg=[ned,ngx,nbphi,nbtheta,nline,4]
  # use reshape on the fly to get datasel
  datasel=(np.reshape(database,shape) )[:,:,:,:nsearch,:,:]   # selection of subset
  for amb in np.arange(0,2):
    ambphib=phib_obs+amb*np.pi/2.
    for ir in np.arange(len(bphir)):   # loop over bphi in cle frame
        ttp= tt * np.sin(bphir[ir])
        diff=r2d*np.abs(ambphib-arctan(ttp)) # this is an array over btheta
        srt = np.argpartition(diff, range(nsearch))[:nsearch]
        srt=np.argsort(diff)[0:nsearch]
        #
        # important to avoid phi=0 and theta = 0 case
        #
        if ir + srt[0] > 0: 
            # again use reshape on the fly to get datasel
            datasel[:,:,ir,:,:,:]=(np.reshape(database,shapeseg)) [:,:,ir,srt,:,:]
            outindex[ir,:]=srt
            #
            # nearest angular separation
            # if larger than cgrid spacing, remove from the search as no angle will
            # work.
            #
            #dr=np.abs(ambphib-bthetar[srt[0]])
            #if dr > dthb: datasel[:,:,ir,:,:,:] *=99.
            #else: datasel[:,:,ir,:,:,:] *=99.
            #
    # now work with reduced dataset with nbtheta replaced by nsearch
    #
    # data3 / datasel is not a memory gobbler.
    #
    data3=np.reshape(datasel,[ned*ngx*nbphi*nsearch,nline,4])
    if verbose ==1:print('Search over theta reduced by a factor', int(nbtheta/nsearch))
    dt= "{:4.3f}".format(time.time()-start)
    if verbose ==1:print(dt,' SECONDS FOR REDUCE (loop phi) CALC')
    #
    # shape of the stokes data sdb ("s database")
    # define new smaller array
    sdb=np.append(data3[:,0,:4],data3[:,1,:4],axis=1)
    del data3
    #
    # Always normalize to the strongest of the raw stokes parameter
    # 
    #norm=sdb[:,0]
    #sdb = (sdb.T/norm.T).T                    # consumes time
    dt= "{:4.3f}".format(time.time()-start)
    if verbose ==1:print(dt,' SECONDS FOR SDB CALC')
    return sdb,outindex



########################################################
###############Parallel implementation #################
########################################################

## global variables can't really be used with parallel python functions unless you aim to keep everything ocnstant.
#@njit(parallel=True) #,forceobj=True)
@njit
def clematch_par(sobs,yobs,database,dbhdr,rms,outdir,maxchisq,nsearch,reduction):#,verbose):
#
# keep these variables in memory, analogous to common fortran
#
    ##PAR in numba, global variablea are defined as uneditable constants; defined all of them inside the optional **keys
    #global database
    #global r2d,nsearch,npartition,reduction, elongation
    #
    # check number of elements of sobs
    #
#ARP: This portion moved to the main solver     
#     numlines=np.int_(len(sobs)/4)
#     if numlines < 2:
#         print('clematch: must be more than 1 observed line to continue')
#         return
    
    #timeentry=time.time()
    #
    # Default values replaced with any keyword values that are set:
#     try: dbdir == ''
#     except:
#         dbdir='dbcle/test/'  
    
#     try: outdir == ''
#     except:    
#         outdir='results/fe13/'
    
#     try: maxchisq == 0
#     except:    
#         maxchisq=1.4 
    
#     try: verbose
#     except:  
#         verbose=0
    #
#     # Set all the keywords to values
#     #
#     nkeys=len(keys)
#     values=["" for x in range(nkeys)]
#     kwords = ["" for x in range(nkeys)]
#     kount=0
#     for arg in keys.values():
#         values[kount] = arg
#         kount+=1
#     kount=0
#     if verbose ==1:print("Keywords...")
#     for arg in keys:
#         kwords[kount] = arg
#         if kwords[kount] == 'dbdir': dbdir=values[kount]
#         if kwords[kount] == 'outdir': outdir=values[kount]
#         if kwords[kount] == 'maxchisq': maxchisq=float(values[kount])
#         if kwords[kount] == 'iplot': iplot=int(values[kount])
#         if kwords[kount] == 'verbose': verbose=int(values[kount])
#         if kwords[kount] == 'database': database=values[kount]
#         #if kwords[kount] == 'r2d' r2d= float(values[kount]) not really needed.
#         if kwords[kount] == 'nsearch': nsearch= int(values[kount])
#         if kwords[kount] == 'npartition': npartition= int(values[kount])
#         if kwords[kount] == 'reduction': reduction= int(values[kount])
#         if kwords[kount] == 'elongation': elongation=float(values[kount])
#         if verbose ==1:print(kwords[kount],'=',values[kount],end=" ")
#         kount+=1
#     if verbose ==1:print("\n")
    #
    #  end of keyword decoding
    #
    # Here are parameters needed but not keywordsabase
    #
#     try: nsearch  #  number of theta angles to keep in reduced search
#     except:
#         nsearch=4            
#         if verbose ==1:print('DEFAULT NSEARCH = ',nsearch)
            
## ARP npartition not really used in code. Don't understand its purpose.   
#     try: npartition # initial guess as to number of solutions with chi^2 < maxchisq
#     except:
#         npartition=8        
#         if verbose ==1:print('DEFAULT NPARTITION = ',npartition)
        
#     try: reduction
#     except:
#         reduction=1     # select database using U/Q (very fast when nsearch << nbtheta)
#         if verbose ==1:print('DEFAULT REDUCTION = ',reduction)    
    r2d=182./np.pi     #; r2d=1.#PAR Not needed as a keyword or parameter.
    #
    # all parameters are now set.
    #    
    #######################################################################
    #
    #
    # Next read dbcgrid parameters from the accompanying db.hdr file in database
    #
    #start=time.time()
    #dbcgrid=readhdr_par(dbdir+'db.hdr')
    #cgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    #    bthetamin, bthetamax, nline, wavel  = parsehdr_par(dbcgrid)
  
    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
        bthetamin, bthetamax, nline, wavel  = dbhdr 

##ARP: This portion moved to the main solver    
#     if numlines != nline:
#         print("Must be the same number of lines in observations and database")
#         print('numlines=',numlines,'nline=',nline)
#         #sys.exit()
        #
        #
        #  Ok so the data are read.        
        #  THIS IS THE START OF THE SEARCH ALGORITHM
        #  get the dbcgrid associated with the data
        #
    ####ARP: NO OBSERVATION NO RUN!!!!!!!!!!!!!    
    if np.isnan(sobs).all():
        #return np.empty((1),dtype=np.int16)[0], [0.,0.,0.,0.,0.,0.,0.,0.],0,[0.,0.,0.,0.,0.,0.]   
        #return np.zeros(1,dtype=np.int_)[0],np.zeros(8,dtype=np.float64)[0],np.zeros(1,dtype=np.float64)[0],np.zeros(6,dtype=np.float64)[0]
        return np.zeros((nsearch),dtype=np.int_),np.zeros((nsearch,8),dtype=np.float64),np.zeros((nsearch),dtype=np.float64),np.zeros((nsearch,8),dtype=np.float64)    
        #return np.zeros((nsearch,10),dtype=np.float64)
    #######################################################################
    #  call routine to deliver the reduced dataset
    #   
    #ARP: need ot be careful here with appending as it only works for 2 lines, not for 3 or n. Superfluous, but still important.
    #sdb1=np.empty((database.shape[0],8),dtype=np.float64) ##append 2 lines database size, 8 observations! 
    #outindex=np.zeros((nbphi,nsearch))
    #print(sdb.shape,outindex.shape)
    #print(database.shape)
    if reduction == 1:
        #sdb=np.empty((database.shape[0],8),dtype=np.float32) ##append 2 lines database size, 8 observations! 
        #outindex=np.zeros([nbphi,nsearch])
        sdb,outindex=get_subset_par(sobs,dbhdr,database,nsearch) 
        #print("cle-match-reduced:",sdb.shape,outindex.shape)
    else:
        #ARP:append is much slower than reshape! will have to see if its possible to switch to a reshape method
        #start=time.time()
        sdb=np.append(database[:,0,:4],database[:,1,:4],axis=1)
        #dt= "{:4.6f}".format(time.time()-start)
        #if verbose ==1:print(dt,' SECONDS FOR SDB APPEND')
        #start=time.time()
        #sdb=np.reshape(database,(database.shape[0],8))
        #dt= "{:4.6f}".format(time.time()-start)
        #if verbose ==1:print(dt,' SECONDS FOR SDB RESHAPE')
        outindex=np.empty((nbphi,nsearch)) ##PAR: Even if never used, numba does not compile if it does not have this array defined because it is named in the code!
        #print("cle-match-full:",sdb.shape,outindex.shape)
    #sdb=np.append(database[:,0,:4],database[:,1,:4],axis=1)
    #sdb,outindex=get_subset_par(sobs,dbhdr,sdb,nsearch) 
    
    #dt= "{:4.3f}".format(time.time()-start)
    #if verbose ==1:print(dt,' SECONDS FOR SDB APPEND')

    # ndof = number of degree of freedom in model: ne, x, bphi, btheta, b
    ndata = nline*4 -1
    ndof=5
    denom=nline*4-1 - ndof  # used below in chi^2
    
#ARP: This portion moved to the main solver and merged with the number of lines check     
#     if denom < 1:
#         print("Number of observables is insufficient to determine model parameters")
#         print("N(data points)=",ndata," N(degrees of freedom)=",ndof)
#         sys.exit()        
        
    # So far so good.. 
    #
    #######################################################################
    ######################################################################
    #
    #  Now find the matched data
    #
    ######################################################################
    #
    #  This is the most time consuming 
    #
    #start=time.time()
    # 
    # We do not use Stokes V to define the geometry
    # Therefore we use a subset of Stokes parameters
    # for  QU (line 1) and IQU (line 2)
    #
    #subset=(1,2,4,5,6)
    # 0 1 2 3 4 5 6 
    #print(sdb[:,subset],sobs[subset],rms[subset])
    
    #diff= (sdb[:,subset]-sobs[subset]) / rms[subset] 
    #dt= "{:4.3f}".format(time.time()-start)
    #print(dt,' SECONDS FOR ALLOC DIFF')
    
    #start=time.time()
    diff=np.empty((sdb.shape[0],6),dtype=np.float64)
    diff[:,0:2]= (sdb[:,1:3]-sobs[1:3]) / rms[1:3] 
    diff[:,2:5]= (sdb[:,4:7]-sobs[4:7]) / rms[4:7] 
    #dt= "{:4.3f}".format(time.time()-start)
    #print(dt,' SECONDS FOR ALLOC DIFF2')
     
    #diff2=np.empty((sdb.shape[0],5),dtype=np.float64)
    #start=time.time()
    #diff2=np.empty((sdb.shape[0],5),dtype=np.float64)
    #diff2=sdb[:,np.arange(sdb.shape[1]) != 7 or np.arange(diff2.shape[1]) != 3 or np.arange(diff2.shape[1]) != 0]
    #diff2=diff2[:,np.arange(diff2.shape[1]) != 3 ]
    #diff2=diff2[:,np.arange(diff2.shape[1]) != 0 ]

    #ARP fast and pomising, does not work
#    diff[:,0:6]= (sdb[:,1:7]-sobs[1:7]) / rms[1:7] 
#    tra=np.empty((sdb.shape[0]),dtype=np.int_)
#     kk=0
#     for jj in range(3,sdb.shape[0]*sdb.shape[1],8):
#         tra[kk]=jj
#         kk+=1
#    diff2=np.empty((sdb.shape[0]*5),dtype=np.float64)
#     diff2=np.reshape(np.delete(diff,tra),(sdb.shape[0],5))
#    tra =np.array([ i for i in range(0,diff.shape[0]*6) if i%6==0])+2
#    diff2=np.delete(diff,tra).reshape(sdb.shape[0],5)
#     tra=np.array((0,3,7),dtype=np.int_)
#     tra2=np.array((0),dtype=np.int_)
#      for jj in range(0,sdb.shape[0]*sdb.shape[1],8):
#          tra2=np.append(tra2,jj)
#          tra2=np.append(tra2,jj+3)
#          tra2=np.append(tra2,jj+7)
        
    #diff2=np.empty((sdb.shape[0],5))
    #diff2=(np.delete(sdb,tra,0))#-np.delete(sobs,tra))
    #print(sdb.shape,np.delete(sdb,tra,0).shape,np.delete(sobs,tra).shape)
    #print(sdb[0,:],diff2[0:8])
    #dt= "{:4.3f}".format(time.time()-start)
    #print(dt,' SECONDS FOR ALLOC DIFF3')
    #print(diff2.shape, diff[432543,:],diff2[432543,:])
    
    #print("true?:",np.array_equal(diff,diff2))
    
    #start=time.time()
    diff*=diff
    #dt= "{:4.3f}".format(time.time()-start)
    #print(dt,' SECONDS FOR PYTHON MULTYPLY')
    
    ##PAR: PYTHON multiply seems faster than numpy multiply
    #start=time.time()
    #diff=np.multiply(diff,diff)
    #dt= "{:4.3f}".format(time.time()-start)
    #print(dt,' SECONDS FOR NP MULTPIPLY')
    
    ##PAR: numpy power seems hideous 15Xtime slower; avoid!
    #start=time.time()
    #diff=np.power(diff,2)
    #dt= "{:4.3f}".format(time.time()-start)
    #print(dt,' SECONDS FOR NP POWER')
    
    
    #print(diff,denom)
    #dt= "{:4.3f}".format(time.time()-start)
    #if verbose ==1:print(dt,' SECONDS FOR CALC DIFF')
    
    # 
    #start=time.time()
    #chisq =np.round_( np.sum( diff, axis=1)/denom,decimals=15) ##PAR: works in object mode
    chisq=np.empty((diff.shape[0]))
    np.round_( np.sum( diff, axis=1)/denom,15,chisq)
    #dt= "{:4.3f}".format(time.time()-start)
    #if verbose ==1:print(dt,' SECONDS FOR CALC CHISQ')
    #print(chisq.shape,chisq)
    #
    #start=time.time()
    ##standard!!!!!!
    #asrt = np.argpartition(chisq, range(nsearch))[:nsearch]
    
    ##ARP: numba compatible loop as np.argpartition does not exist!
    ##numpy partition+where would be compatible to numba
#     bla = np.partition(chisq,range(nsearch))[:nsearch]
#     asrt1=np.empty(nsearch,dtype=np.int32)
#     asrt=np.zeros(nsearch,dtype=np.int_)
#     gigi=len(np.argwhere(chisq == bla[0] ))
#     print(gigi,np.where(chisq == bla[0] )[0:2][0])
#     print("{:2.16f}".format(bla[0]),"{:2.16f}".format(bla[1]),"{:2.16f}".format(bla[2]),"{:2.16f}".format(bla[3]))

#     for i in range(0,nsearch,gigi):
#         #print(np.where([5,3,1,9,14,7,6,0] == bla[i] )[0][0])
#         asrt1[i:i+gigi]=np.where(chisq == bla[i] )[0:gigi][0]
    ## too fussy to make work. Also slow
    ####----------------------------
    ## argmin+where approach
#     j=0
#     for i in range(0,nsearch):
#         #print(np.min(chisq),np.where(chisq == np.min(chisq))[0][0])
#         #a=np.where(chisq == np.min(chisq))[0][0]
#         a=np.argmin(chisq)
#         print(a,j)
#         if a > asrt2[i-1]:
#             asrt2[i]=a+j
#             chisq=np.delete(chisq,a)
#             j+=1
#         else: 
#             asrt2[i]=a+j
#             chisq=np.delete(chisq,a)    #
    ##too convoluted; also could not set up the indices to make it work.
    ####----------------------------
    ##ARP: chisq needs to be passed as np.copy regardless if the sorting is done locally or in the subfunction.
    ## There is no explanation for htis behaviour
    asrt=subst_sort_par(np.copy(chisq),nsearch)  #ARP: made a function to do the simple sort needed here.      
    
    #dt= "{:4.3f}".format(time.time()-start)
    #if verbose ==1:print(dt,' SECONDS FOR SORT CHISQ')
    #ix=np.zeros(1,dtype=np.int_)
    ix=asrt[0]
    #print(ix.type())
    #ixr=ix
    #
    #  return indices and chisq of nearest solutions
    #
    si=0
    #ixr=np.zeros(1,dtype=np.int_)
    #ixr=np.float_()
    ixr=asrt[si]
    chisqmin=chisq[ixr]
    #chisquare=np.float_("NaN")

    #dt= "{:4.3f}".format(time.time()-timeentry)
    #if verbose ==1:print(dt,' SECONDS FOR SOLUTION')
    
    #print("manele:",chisq[ixr], maxchisq,si,nsearch)
    while chisq[ixr] < maxchisq and si < nsearch:       
        ixr=asrt[si]
        dd=chisq[ixr]#
        #print("lat:",dd,maxchisq)
        if(dd < maxchisq): 
            #
            # here if reduced we must get the original value of ix
            if reduction  == 1:
                cgridnew=dbcgrid + 0.0
                cgridnew[3]=nsearch
                i,j,k,l=params_par(ixr,cgridnew)
                newl=np.int_(outindex[k,l])
                ix=invparams_par(i,j,k,newl,dbcgrid) # = original value ix
            #
            # Magnetic field strength from V/I of first line 
            bfield=sobs[3]/sdb[ixr,3]
            #
            if si == 0:                
                index=np.zeros(1,dtype=np.int_)
                index[0]=ix
                sfound=np.zeros((nline*4),dtype=np.float64)
                #print('si is:',si,'---',sfound.shape)
                sfound=sdb[ixr,:]  # use reduced data index
                chisquare=np.zeros(1,dtype=np.float64)
                chisquare[0]=dd
                ph=phys_par(ix,yobs,dbhdr,bfield)
                #print("hallo!")
            else:
                index = np.append(index, ix)
                sfound=np.append(sfound,sdb[ixr,:])
                #print('si is:',si,'---',sfound.shape,'...',sdb[ixr,:].shape)
                chisquare=np.append(chisquare,dd)
                ph=np.append(ph,phys_par(ix,yobs,dbhdr,bfield))
                #print("else hallo!",index.shape[0])
        #
        si+=1
    
#     if not index:
#         return ixr, [0.,0.,0.,0.,0.,0.,0.,0.],chisq[ixr],[0.,0.,0.,0.,0.,0.]
    #
    # reshape arrays so they are readable by a human
    #
    #print("I am here!")
    dims=np.int_(index.shape[0])
    if dims > 1:
        sshape=(dims,np.int_(nline*4))
        sfound1=np.reshape(sfound,sshape) ##
        pshape=(dims,8) #ARP: 8 is the numer of outputs resulted from phys_par subroutine
        ph1=np.reshape(ph,pshape)        
        #print("dims:",sshape,sfound.shape,sfound1.shape,pshape,ph.shape,ph1.shape)
        return index, sfound1, chisquare, ph1
        #print(np.concatenate((index.reshape(nsearch,1), chisquare.reshape(nsearch,1), ph1)).shape)
        #return np.concatenate((index, chisquare, ph1))
    #gc.collect()
    else:
        #print('nodims:',nline,sfound.shape,ph.shape)
        ##ARP: This needs to return the same number of nsearch leements regardless of whats inside
        sfound2=np.zeros((nsearch,8),dtype=np.float64)
        #print(dims,sfound2.shape,sfound.shape,np.int_(nline*4),nline)
        sfound2[0,:]=sfound
        ph2=np.zeros((nsearch,8),dtype=np.float64)
        ph2[0,:]=ph
        #return index, sfound.reshape((1,nline*4)), chisquare, ph.reshape((1,6)) ##ARP: the reshapes are due to the fact that all 3 return lynes are forced to have the veriables of the same type!
        return index, sfound2, chisquare, ph2
        #return np.concatenate((index.reshape(nsearch,1), chisquare.reshape(nsearch,1), ph2))

#################################################################################################################  
##refined allocation from newer implementation; but does not get rid of the multiple array allocations
#             if si == 0:               # Initialize arrays on first iteration             
#                 index=np.zeros(1,dtype=np.int_)
#                 index[0]=ix
#                 sfound=np.zeros(8,dtype=np.float64)   ##8--> fixed for a two line obs at this point
#                 sfound=sdb[ixr,:]     # use reduced data index
#                 chisquare=np.zeros(1,dtype=np.float64)
#                 chisquare[0]=dd
#                 ph=phys_par(ix,yobs,dbhdr,bfield)
#             else:                      # append infor for each iterationn 
#                 index = np.append(index, ix)
#                 sfound= np.append(sfound,sdb[ixr,:]) 
#                 chisquare=np.append(chisquare,dd)
#                 ph=np.append(ph,phys_par(ix,yobs,dbhdr,bfield))
       
#         si+=1
#     #print(index.shape)
#     #######################################################################     
#     ## reshape arrays to a human readable format and return

#     #Both instances need to return the same number of nsearch elements regardless of whats inside
#     dims=np.int_(index.shape[0])
#     #print("dims",dims,nsearch)
#     if dims == nsearch:
#         #sshape=(nsearch,8)
#         #sfound1=np.reshape(sfound,sshape) ##
#         #pshape=(nsearch,9) # 9 is the numer of outputs resulted from phys_par subroutine
#         #ph1=np.reshape(ph,pshape)        
        
#         if verbose >=2: print("Warning: returned solutions less than desired nsearch")
#         #if verbose >=1: print("Return dims:",index.shape,sfound1.shape,chisquare.shape,ph1.shape)
        
#         #print("case 2:",index.shape, sfound1.shape, chisquare.shape, ph1.shape)    
#         #return index, sfound1, chisquare, ph1 
#         return index,sfound.reshape(nsearch,8),chisquare,ph.reshape(nsearch,9)
#         #return index,sfound,chisquare,ph
#     else:
#         index2=np.zeros((nsearch),dtype=np.int_)
#         index2[0:dims]=index
                              
#         sfound2=np.zeros((nsearch,8),dtype=np.float64)
#         sfound2[0:dims,:]=sfound
                              
#         chisquare2=np.zeros((nsearch),dtype=np.float64)
#         chisquare2[0:dims]=chisquare             
        
#         ph2=np.zeros((nsearch,9),dtype=np.float64) # 9 is the numer of outputs resulted from phys_par subroutine
#         ph2[0:dims,:]=ph
#         #print("case 3:",index.shape, sfound2.shape, chisquare.shape, ph2.shape)
#         return index2, sfound2, chisquare2, ph2
#################################################################################################################        
        
        
        
# @jit(parallel=True,forceobj=True)
# def sdbselect_par(yobs,yobs_prev,dbdir):
#     #
#     verbose=1
#     #  decision- read a new data file?
#     # ARP: so, check if the previous iteration elongation is the same, in that case, do not read the database again.
#         #
 
#     #  get the elongation file name
#     #
#     fildb,iy=find_elongation_par(dbdir,yobs)
#     fildb_prev,iy=find_elongation_par(dbdir,yobs_prev)
#     return fildb,fildb_prev

######################################################################
## 
## Simplest possible substitution partial sort. It just returns the first nsearch sorted indexes similar to np.argpartition. 
## It is numba non-python compatible!      
#@njit(parallel=True)#,forceobj=True)
@njit
def subst_sort_par(arr,nsearch):
    asrt=np.zeros((nsearch),dtype=np.int_)
    for i in range(nsearch):
        a=np.int_(np.argmin(arr[i:])+i)
        asrt[i]=a 
        arr_temp=arr[i]
        arr[i]=arr[a]
        arr[a]=arr_temp
    return asrt  
    ##works like a charm
######################################################################
#
#@njit(parallel=True)#,forceobj=True)
@njit
def get_subset_par(sobs,dbhdr,database,nsearch):
    # 
    # returns a subset of the sdb array compatible with yobs
    #
    #global database
    #global r2d,nsearch,npartition,reduction, elongation
    verbose=0
    r2d=180./np.pi 
    #start=time.time()

    #cgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    #bthetamin, bthetamax,  nline, wavel  = parsehdr_par(dbcgrid)
    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    bthetamin, bthetamax, nline, wavel  = dbhdr 

    kk=np.arange(0,nbphi)
    bphir  = bphimin + kk*(bphimax-bphimin)/(nbphi-1)  # cle phi array

    ll=np.arange(0,nbtheta)
    bthetar=bthetamin + ll*(bthetamax-bthetamin)/(nbtheta-1) #cle theta array
    dthb=(bthetamax-bthetamin)/nbtheta
    tt = np.tan(bthetar)  #cle tan theta array
    #print(database.shape[0]*database.shape[1]*database.shape[2])
    shape=(ned,ngx,nbphi,nbtheta,nline,4)
    #shape=np.array((ned,ngx,nbphi,nbtheta,nline,4),dtype=np.float32)
    # use reshape on the fly to get datasel;;;;
    
    database=np.ascontiguousarray(database)
    #datasel=(np.reshape(np.ascontiguousarray(database),shape) )[:,:,:,:nsearch,:,:]   # selection of subset
    datasel=(np.reshape(database,shape) )[:,:,:,:nsearch,:,:]   # selection of subset#datasel=np.reshape(database,shape) # selection of subset
    outindex=np.zeros((nbphi,nsearch))
    #
    gnew=dbcgrid + 0.
    gnew[3]=nsearch
    #
    #   Below, the relation is used: tan Phi_B = sin phi * tan theta
    #   here is tan Phi_B:
    #
    phib_obs=-0.5*np.arctan2( sobs[2],sobs[1] )
    tphib_obs= np.tan(phib_obs)
    #
    # Find those indices compatible with phib observed
    #
    #shapeseg=(ned,ngx,nbphi,nbtheta,nline,4)
    #shapeseg=np.array((ned,ngx,nbphi,nbtheta,nline,4))
    # use reshape on the fly to get datasel
    #datasel=(np.reshape(database,shape) )[:,:,:,:nsearch,:,:]   # selection of subset
    for amb in prange(0,2):
        ambphib=phib_obs+amb*np.pi/2.
        for ir in prange(0,len(bphir)):   # loop over bphi in cle frame
            ttp= tt * np.sin(bphir[ir])
            diff=r2d*np.abs(ambphib-arctan(ttp)) # this is an array over btheta
            ###ARP: Sort manually to appease numba gods!
            #srt = np.argpartition(diff, range(nsearch))[:nsearch]
            #srt=np.argsort(diff)[0:nsearch] ##ARP: argsort works with numba non-python but it stil ltakes an entire array and slices the first nindices. 
##ARP: Manual sort replased by function subst_sort_par 
#             diff1=np.copy(diff)
#             srt=np.zeros(nsearch,dtype=np.int32)
#             for i in range(0,nsearch):
#                 a=np.argmin(diff1[i:])+i
#                 srt[i]=a 
#                 bu=diff1[i]
#                 diff1[i]=diff1[a]
#                 diff1[a]=bu
            srt=subst_sort_par(np.copy(diff),nsearch)            
            #
            # important to avoid phi=0 and theta = 0 case
            #
            if ir + srt[0] > 0: 
                # again use reshape on the fly to get datasel
                #datasel[:,:,ir,:,:,:]=(np.reshape(database,shapeseg)) [:,:,ir,srt,:,:]
                #outindex[ir,:]=srt
                #ARP:The for+enumerate comes from numba slicing
                
                for jj in prange(0,len(srt)): 
                    #tt1=datasel[:,:,ir,jj,:,:].copy()
                    #tt1=(np.reshape(np.ascontiguousarray(database),shape))[:,:,ir,srt[jj],:,:].copy()
                    tt1=(np.reshape(database,shape))[:,:,ir,srt[jj],:,:].copy()
                    datasel[:,:,ir,jj,:,:]=tt1
                    outindex[ir,jj]=srt[jj]
                
                #
                # nearest angular separation
                # if larger than cgrid spacing, remove from the search as no angle will
                # work.
                #
                #dr=np.abs(ambphib-bthetar[srt[0]])
                #if dr > dthb: datasel[:,:,ir,:,:,:] *=99.
                #else: datasel[:,:,ir,:,:,:] *=99.
                #
    # now work with reduced dataset with nbtheta replaced by nsearch
    # data3 / datasel is not a memory gobbler.
    #
    #gg=ned*ngx*nbphi*nsearch
    shapesig=(ned*ngx*nbphi*nsearch,nline,4)
    #shapesig=np.array((ned*ngx*nbphi*nsearch,nline,4),dtype=np.float32)
    #data3=np.zeros((gg,nline,4),dtype=np.float64)
    #print("internal shapes:",datasel.shape,gg*2*4)
    data3=np.reshape(np.ascontiguousarray(datasel),shapesig)
    if verbose ==1:print('Search over theta reduced by a factor', np.int_(nbtheta/nsearch))
    #dt= "{:4.3f}".format(time.time()-start)
    #if verbose ==1:print(dt,' SECONDS FOR REDUCE (loop phi) CALC')
    #
    # shape of the stokes data sdb ("s database")
    # define new smaller array
    sdb=np.append(data3[:,0,:4],data3[:,1,:4],axis=1)
    #del data3
    #
    # Always normalize to the strongest of the raw stokes parameter
    # 
    #norm=sdb[:,0]
    #sdb = (sdb.T/norm.T).T                    # consumes time
    #dt= "{:4.3f}".format(time.time()-start)
    #if verbose ==1:print(dt,' SECONDS FOR SDB CALC')
    #print(sdb.shape,outindex.shape)
    return sdb,outindex



######################################################################
#
# @jit
# def dcompress_par(i):
# # Constants here must correspond to those in dbe.f in the CLE main directory
# # 32767 is the limit of 4-byte ints
# # 15 is the number of orders of magnitude for the range of intensities
# #
#     verbose=0
#     cnst=-2.302585092994046*15./32767.
#     start_dc=time.time()
#     strt=time.time()
#     if verbose ==1:print(np.int(sys.getsizeof(i)/1.e6) , 'MB in database file')
#     negv=np.flatnonzero(i  < 0 )
#     dt= "{:4.3f}".format(time.time()-strt)
#     if verbose ==1:print('DCOMPRESS:', dt,' SECONDS FOR WHERE')
    
#     strt=time.time()
#     #f=np.abs(i.astype(float))
#     f=np.abs(i)*cnst

#     #t = timeit.Timer(lambda: np.exp(f))
#     #print('--------------------TIMER np.exp ',t.timeit(1))

#     f=ne.evaluate("exp(f)")
#     #import pdb; pdb.set_trace()

#     dt= "{:4.3f}".format(time.time()-strt)
#     if verbose ==1:print('DCOMPRESS:', dt,' SECONDS FOR EXP')

#     dt= "{:4.3f}".format(time.time()-start_dc)
#     if size(negv) > 0:f[negv]=-f[negv]
#     if verbose ==1:print('DCOMPRESS:', dt,' SECONDS TOTAL')
#     #del negv
#     return f
#
######################################################################
#@jit(parallel=True,forceobj=True,looplift=True)
def sdbread_par(fildb,dbhdr):
    #
    verbose=0
    start=time.time()    
    #dbcgrid=readhdr_par(dbdir+'db.hdr')
    #cgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    #    bthetamin, bthetamax, nline, wavel  = parsehdr_par(dbcgrid)   
    
    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    bthetamin, bthetamax, nline, wavel  = dbhdr 
    #
    #  here reading data is done and stored in variable database of size (n,2*nline)
    #
    if verbose ==1:print('CLEMATCH: READING DATABASE...')
    #
    #filn = open(fildb, 'rb')
    #
    start1=time.time()
    db=dcompress(np.fromfile(fildb, dtype=np.int16))
    #filn.close()
    if verbose == 1:
        dt= "{:4.3f}".format(time.time()-start1)
        print(dt,' SECONDS FOR DCOMPRESS')
        print("DB file location\n", fildb,"\n")    
    #start=time.time()
    #shape=(ned*ngx*nbphi*nbtheta,nline,4)
    #database=np.reshape(database,shape)
    if verbose == 1:
        dt= "{:4.3f}".format(time.time()-start)
        print(dt,' SECONDS FOR TOTAL DB READ/RESHAPE')
        #    
    #return database
    return np.reshape(db,(ned*ngx*nbphi*nbtheta,nline,4))
    
######################################################################
#ARP: slightly simplified version with par2. decreases runtime by about 10%. 
# @jit(parallel=True,forceobj=True)    #
# def find_elongation_par(dbdir,y):
# #
# # returns filename and index of elongation y in database
# #
#     names=glob.glob(dbdir+"DB*.DAT")
#     nm=len(names)
#     numbers=np.zeros(nm,dtype=np.int_)
#     for i in range(0,nm):
#         nn=str.find(names[i],'DB')
#         string=names[i]
#         string=string[nn+2:nn+6]
#         numbers[i]=int(string)
#     #
#     #Use naming convention
#     #
#     y_of_files= 1.+ numbers/1000.
#     yi=np.argmin(np.abs(y_of_files-y))
#     file=names[yi]
#     iy=numbers[yi]  # set it to the number ID of the file
#     strc="{:0.3f}".format(y_of_files[yi]) 
#     #print("Elongation =",y,"nearest DB file is",file)
#     return file

@jit(parallel=True,forceobj=True)    #
def find_elongation_par2(dbdir,y):
#
# returns filename and index of elongation y in database
#
    names=glob.glob(dbdir+"DB*.DAT")
    #nm=len(names)
    numbers=np.empty(len(names),dtype=np.int_)
    nn=str.find(names[0],'DB')
    for i in prange(0,len(names)):
        #nn=str.find(names[i],'DB')
        #string=names[i][nn+2:nn+6]
        #string=string[nn+2:nn+6]
        #numbers[i]=np.int_(string)
        numbers[i]=np.int_(names[i][nn+2:nn+6])
    #
    #Use naming convention
    #
    y_of_files= 1.+ numbers/1000.
    yi=np.argmin(np.abs(y_of_files-y))
    #file=names[yi]
    #iy=numbers[yi]  # set it to the number ID of the file
    #strc="{:0.3f}".format(y_of_files[yi]) 
    #print("Elongation =",y,"nearest DB file is",file)
    return names[yi]#file

@jit(parallel=True,forceobj=True,looplift=True)    #
def find_elongation_num_par(y,numbers):
#
# returns filename and index of elongation y in database
# 
    #names=glob.glob(dbdir+"DB*.DAT")
    #nm=len(names)
    #numbers=np.empty(len(names),dtype=np.int_)
    #nn=str.find(names[0],'DB')
    #for i in prange(0,len(names)):
        #nn=str.find(names[i],'DB')
        #string=names[i][nn+2:nn+6]
        #string=string[nn+2:nn+6]
        #numbers[i]=np.int_(string)
        #numbers[i]=np.int_(names[i][nn+2:nn+6])
        #numbers[i]=np.array((names[i][nn+2:nn+6])).astype(np.int_)
        #
    #Use naming convention
    #
    #y_of_files= 1.+ numbers/1000.
    #yi=np.argmin(np.abs(y_of_files-y))
    #file=names[yi]
    #iy=numbers[yi]  # set it to the number ID of the file
    #strc="{:0.3f}".format(y_of_files[yi]) 
    #print("Elongation =",y,"nearest DB file is",file)
    #return yi#names[yi]#file
    return np.argmin(np.abs(1.+ numbers/1000-y))
#@jit#(parallel=True)#,forceobj=True)    #
#@jit(forceobj=True)
#def find_elongation_str_par(dbdir,yi): ## Removed to get rid of glob inside fast functions
#
# returns filename and index of elongation y in database
#
#    #names=glob.glob(dbdir+"DB*.DAT")
#    #return names[yi]#file
#    return glob.glob(dbdir+"DB*.DAT")[yi]
def file_preprocess(dbdir):   ## function to prepare the database directory files outside of the numba jitted functions as string operations do not play nice with numba.
    names=glob.glob(dbdir+"DB*.DAT")
    nn=str.find(names[0],'DB')
    numbers=np.empty(len(names),dtype=np.int_)
    for i in range(0,len(names)):
        numbers[i]=np.int_(names[i][nn+2:nn+6])
    return names,numbers
######################################################################


# @jit(parallel=True,forceobj=True)    #
# #@jit(forceobj=True)#(parallel=True)
# def find_sel_db(dbdir,y,dbhdr,sx,sy):
#     fildb=np.empty((0),dtype='S'+np.int_(len(find_elongation_par2(dbdir,y[0,0]))).astype(str))
#     filenc=np.empty((sx,sy),dtype=np.int_)
#     print(fildb.shape,np.int_(len(find_elongation_par2(dbdir,y[0,0]))).astype(str))
#     for xx in range(0,sx):
#         for yy in prange(0,sy): 
#             #start0=time.time()
#             #filtemp=cle202_utils.find_elongation_par(dbdir,yobs_a[xx,yy])  
#             #dt= "{:4.6f}".format(time.time()-start0)
#             #print(dt,' SECONDS FOR elong1 ',dt) 
#             #start0=time.time()
#             #filtemp=find_elongation_par2(dbdir,y[xx,yy])  
#             #dt= "{:4.6f}".format(time.time()-start0)
#             #print(dt,' SECONDS FOR elong2 ',dt)
#             #print(filtemp)   
            
#             #start0=time.time()
#             if find_elongation_par2(dbdir,y[xx,yy]) not in fildb:
#                 fildb=np.append(fildb,filtemp)     
#             filenc[xx,yy]=np.argwhere(fildb==filtemp)
#             #print(filenc[xx,yy])
#             #dt= "{:4.6f}".format(time.time()-start0)
#             #print(dt,' SECONDS FOR append ',dt)  
#     #dt= "{:4.6f}".format(time.time()-start0)
#     #print(dt,' SECONDS FOR append ',dt)    
#     print(fildb)
#     #start0=time.time()
    
#     database=np.empty((dbhdr[1]*dbhdr[2]*dbhdr[3]*dbhdr[4],dbhdr[12],4,np.max(filenc)+1),dtype=np.float64)
#     print(database.shape)
#     for ii in prange(0,np.max(filenc)+1):
#         database[:,:,:,ii]=sdbread_par(fildb[ii],dbdir,dbhdr)
#     #dt= "{:4.3f}".format(time.time()-start0)
#     #print(dt,' SECONDS FOR dbread ',dt)
#     return filenc,database


@jit(parallel=True,forceobj=True,looplift=True)    #
def find_sel_db(y,dbhdr,names,numbers,sx,sy):
    #fildb=np.empty((0),dtype='S'+np.int_(len(find_elongation_par2(dbdir,y[0,0]))).astype(str))
    fenc=np.empty((sx,sy),dtype=np.int_)
    #filenc1=np.empty((sx,sy),dtype=np.int_)
    #print(fildb.shape,np.int_(len(find_elongation_par2(dbdir,y[0,0]))).astype(str))
    #start0=time.time()
    #start0=time.time() 
    #for xx in prange(sx):
    #    for yy in prange(sy): 
            #start0=time.time()
            #filtemp=find_elongation_par2(dbdir,y[xx,yy])  
            #dt= "{:4.6f}".format(time.time()-start0)
            #print(dt,' SECONDS FOR elong classic ',dt) 
            #start0=time.time()  
            #dt= "{:4.6f}".format(time.time()-start0)
            #print(dt,' SECONDS FOR elong2 ',dt)
            #print(filtemp)   
    #        filenc[xx,yy]=np.argmin(np.abs(1.+ numbers/1000-y[xx,yy]))
            #print(filenc[xx,yy])
    #dt= "{:4.6f}".format(time.time()-start0)
    #print(dt,' SECONDS FOR elong direct ',dt)  
    start0=time.time()
    for xx in range(sx):
        for yy in prange(sy):                     
            fenc[xx,yy]=find_elongation_num_par(y[xx,yy],numbers)
    dt= "{:4.6f}".format(time.time()-start0)
    print(dt,' SECONDS FOR elong num ',dt)  
    
    #print(np.array_equal(filenc,filenc1))
            
            
    #dt= "{:4.6f}".format(time.time()-start0)
    #print(dt,' SECONDS FOR filenc ',dt)    
    
    start0=time.time()
    #print(np.unique(filenc),np.unique(filenc).shape)
    sz=np.unique(fenc)
    #database=np.empty((dbhdr[1]*dbhdr[2]*dbhdr[3]*dbhdr[4],dbhdr[12],4,sz.shape[0]),dtype=np.float64)
    #print(database.shape)
    print("read actual db")
    #database=sdbread_par(names[sz[0]],dbhdr)
    ## numpy large array implementation does not parralelize properly leading to a 5x increase in runtime per 1024 calculations
    ## reverted to use a list to feed the database set to the calculation
    database0=[None]*sz.shape[0]
    for ii in prange(sz.shape[0]):
        database0[ii]=sdbread_par(names[sz[ii]],dbhdr)
        #database[:,:,:,ii]=sdbread_par(names[sz[ii]],dbhdr)
        #database=np.append(database,sdbread_par(names[sz[ii]],dbhdr)).reshape(dbhdr[1]*dbhdr[2]*dbhdr[3]*dbhdr[4],dbhdr[12],4,ii+1)
    ## standard reflacted lists will be deprecated in numba 0.53; currently running 0.49
    database = List()
    [database.append(x) for x in database0]
    
    dt= "{:4.3f}".format(time.time()-start0)
    print(dt,' SECONDS FOR dbread ',dt)
    
    
    
    start0=time.time()
    cnt=0
    for i in sz:
        fenc[np.where(fenc == i)] = cnt
        cnt+=1    
    dt= "{:4.3f}".format(time.time()-start0)
    print(dt,' SECONDS FOR encoding ',dt)
    return fenc,database




######################################################################
#
@jit(forceobj=True)
def readhdr_par(file):
#
# just read the header of the DB*DAT data cubes
#
    g=np.fromfile(file,dtype=np.float32,sep=' ')
    return g

#
######################################################################
#
#@njit(parallel=True)
@njit
def parsehdr_par(g):
#
# parse the header of the DB*DAT data cubes
#
# two kinds of data are returned
# 1. the linear coefficents of the form min, max, n
#    - x,theta,phi
# 2. logarithmically spaced parameters
#    - xed  (electron density array)
#   for these also the min and max and number are also returned.
    ned=np.int_(g[0]) 
    ngx=np.int_(g[1])
    nbphi=np.int_(g[2])
    nbtheta=np.int_(g[3])
    # log of Ne
    emin=np.float64(g[4])
    emax=np.float64(g[5])
    xed=lcgrid_par(emin,emax,ned)
    gxmin=np.float64(g[6])
    gxmax=np.float64(g[7])
    bphimin=np.float64(g[8])
    bphimax=np.float64(g[9])
    bthetamin=np.float64(g[10])
    bthetamax=np.float64(g[11])
    nline=np.int_(g[12])
    wavel=np.empty(nline,dtype=np.float64)
    for k in prange(0,nline): wavel[k]=g[13+k]    
    #return np.array([g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,bthetamin, bthetamax, nline, wavel]) 
    #hdrparams=np.array([ned, ngx, nbphi, nbtheta, gxmin, gxmax, bphimin, bphimax,bthetamin, bthetamax, nline, wavel ])#,])
    return g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
        bthetamin, bthetamax, nline, wavel 
    #return hdrparams
######################################################################
#
#@njit(parallel=True)       ### grid for ned densities from mn to mx logarithms of density.
@njit
def lcgrid_par(mn,mx,n):
    #
    ff=10.** ((mx-mn)/(n-1))
    xf=np.empty(n,np.float64)
    xf[0]=10.**mn
    for k in prange(1,n): xf[k]=xf[k-1]*ff
    return xf
######################################################################
#
@njit #(parallel=True) nothing to parralelize here
def params_par(index,dbcgrid):
    #
    # for index, get i,j,k,l   indices
    #
    ned=np.int_(dbcgrid[0]) 
    ngx=np.int_(dbcgrid[1])
    nbphi=np.int_(dbcgrid[2])
    nbtheta=np.int_(dbcgrid[3])
    n5=( ngx*nbphi*nbtheta)
    n4=(     nbphi*nbtheta)
    n3=(           nbtheta)
    #
    i = index   //        n5
    j = index   -     i * n5
    j = j      //         n4
    k = index   -     i * n5 - j * n4
    k = k      //         n3
    l=  index   -     i * n5 - j * n4 - k * n3

    #return np.array((i,j,k,l),dtype=np.int_)
    return i,j,k,l
    
######################################################################
#
@njit #(parallel=True) nothing to parralelize here
def invparams_par(i,j,k,l,dbcgrid):
    #
    # for i,j,k,l get index
    #
    ned=np.int_(dbcgrid[0]) 
    ngx=np.int_(dbcgrid[1])
    nbphi=np.int_(dbcgrid[2])
    nbtheta=np.int_(dbcgrid[3])


    return np.int_(i*ngx*nbphi*nbtheta + j*nbphi*nbtheta + k*nbtheta + l)

######################################################################
#
@njit #(parallel=True) nothing to parralelize here
def electron_density_par(r):
    baumbach = 1.e8*(0.036/r**1.5 + 1.55/r**6.)
    hscale=   7.18401074e-02  # 50 Mm
    #n0=np.float64(3.e8*np.exp(- (r-1.)/hscale) + baumbach)
    #return n0
    return np.float64(3.e8*np.exp(- (r-1.)/hscale) + baumbach)
#
######################################################################
#
@njit #(parallel=True,forceobj=True) nothing to parralelize here
def phys_par(index,gy,dbhdr,b):
    #
    #i,j,k,l = params_par(index,dbcgrid) not needed
    phs=physa_par(index,gy,dbhdr,b)
    bphi=phs[3]
    btheta=phs[4]
    bx=b*np.sin(btheta)*np.cos(bphi)
    by=b*np.sin(btheta)*np.sin(bphi)
    bz=b*np.cos(btheta)
    #out=np.array((phs[0],phs[1],phs[2],bx,by,bz),dtype=np.float64)
    #return out
    return np.array((phs[0],phs[1],phs[2],phs[3],phs[4],bx,by,bz),dtype=np.float64)
#

######################################################################
#
@njit #(parallel=True,forceobj=True)
def physa_par(index,gy,dbhdr,b):
    #

    #
    #grd,ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax,\
    #bthetamin, bthetamax, nline,wavel  = parsehdr_par(grd)
    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    bthetamin, bthetamax, nline, wavel  = dbhdr 
    #
    i,j,k,l = params_par(index,dbcgrid)
    #
    gx = gxmin + j*(gxmax-gxmin)/(ngx-1)
    #
    bphi  = bphimin + k*(bphimax-bphimin)/(nbphi-1)
    btheta  = bthetamin + l*(bthetamax-bthetamin)/(nbtheta-1)
    #
    # log of Ne
    #
    ed = np.float64(xed[i]* electron_density_par(np.sqrt(gy*gy+gx*gx)))
    #out=np.array((np.log10(ed),gy,gx,bphi,btheta,b),dtype=np.float64)
    #return out
    return np.array((np.log10(ed),gy,gx,bphi,btheta,b),dtype=np.float64)
#
