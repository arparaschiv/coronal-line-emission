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

global database
global r2d,nsearch,npartition,reduction, elongation

np.set_printoptions(precision=4,suppress=False)

#
######################################################################
#
def dcompress(i):
# Constants here must correspond to those in dbe.f in the CLE main directory
# 32767 is the limit of 4-byte ints
# 15 is the number of orders of magnitude for the range of intensities
#
    cnst=-2.302585092994046*15./32767.
    start_dc=time.time()
    strt=time.time()
    print(int(sys.getsizeof(i)/1.e6) , 'MB in database file')
    negv=np.flatnonzero(i  < 0 )
    dt= "{:4.3f}".format(time.time()-strt)
    print('DCOMPRESS:', dt,' SECONDS FOR WHERE')
    
    strt=time.time()
    #f=np.abs(i.astype(float))
    f=np.abs(i)*cnst

    #t = timeit.Timer(lambda: np.exp(f))
    #print('--------------------TIMER np.exp ',t.timeit(1))

    f=ne.evaluate("exp(f)")
    #import pdb; pdb.set_trace()

    dt= "{:4.3f}".format(time.time()-strt)
    print('DCOMPRESS:', dt,' SECONDS FOR EXP')

    dt= "{:4.3f}".format(time.time()-start_dc)
    if size(negv) > 0:f[negv]=-f[negv]
    print('DCOMPRESS:', dt,' SECONDS TOTAL')
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
  print("Elongation =",y,"nearest DB file is",file)
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
    print('Delta ',delta, 'exist_data',exist_data)
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
          newl=int(outindex[k,l])
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
    return index, sfound, chisquare,ph
#

######################################################################
def get_subset(sobs,dbcgrid):
  # 
  # returns a subset of the sdb array compatible with yobs
  #
  global database
  global r2d,nsearch,npartition,reduction, elongation
  verbose=1

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
