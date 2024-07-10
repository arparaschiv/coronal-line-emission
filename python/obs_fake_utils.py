import numpy as np
import glob
import time
import numexpr as ne
from pylab import *
from numba import jit, prange

@jit(forceobj=True)
def fake_array(dbdir,xxl,yyl,x0,y0,dx,dy,counts):
    
    yfake=np.empty((xxl,yyl),dtype=np.float_)
    for xx in range(xxl):
        for yy in prange(yyl):
            yfake[xx,yy]=np.sqrt((x0+(xx*dx))**2+(y0+(yy*dy))**2 )
    print("done yfake")
## create a fake sobs array
    sfake=np.empty((xxl,yyl,8),dtype=np.float_)       
    rms=np.empty((xxl,yyl,8),dtype=np.float_) 
    oindex=np.empty((xxl,yyl),dtype=np.int_) 


## THIS IS SLOW!!!!
    verbose=1
    for xx in range(xxl):
        print("Executing ext. loop:",xx," of ",xxl)
        for yy in prange(yyl):
            counts1=counts+(1e-4*counts*np.random.normal(0.,1))   ##add another small fluctuation to counts to account for pixel to pixel variation in rms       
            oindex[xx,yy]= np.random.randint(10*60*180*90)   ## random index from database; calculation size ##ned*ngx*nbphi*nbtheta
            sfake[xx,yy,:],rms[xx,yy,:],oindex[xx,yy]=obs_fake(dbdir,yfake[xx,yy],oindex[xx,yy],counts1,verbose)

    return sfake,yfake,rms,xxl,yyl,oindex
## all required variables(sobs_a,yobs_a,rms,xxl,yyl) are computed.

def obs_fake(dbdir,yobs,oindex,counts,verbose,**keys):
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
    if verbose == 1: start=time.time()
    nkeys=len(keys)
    values=["" for x in range(nkeys)]
    kwords = ["" for x in range(nkeys)]
    kount=0
    for arg in keys.values():
      values[kount] = arg
      kount+=1
    kount=0

    if (verbose == 1): print("Keywords...")
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
    if verbose ==1: print("{:4.6f}".format(time.time()-start),'FAKE_OBS: TOTAL SECONDS FOR PREAM ')
    
    fildb,iy=find_elongation(dbdir,yobs)
    filn = open(fildb, 'rb')
    #
    database=dcompress( np.fromfile(file=filn, dtype=np.int16) )
    filn.close()
    # reshape on the fly
    i,j,k,l= params(oindex,g)
    shap=[ned,ngx,nbphi,nbtheta,nline,4]
    tmp=np.reshape(database,shap)
    if verbose == 1: print("{:4.6f}".format(time.time()-start),'FAKE_OBS: READING DATABASE FILE & DECOMPRESSING DATA')
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
    if verbose == 1: print("{:4.6f}".format(time.time()-start),'FAKE_OBS: MAKE SOBS')
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
    if verbose ==1:print('Observed B is ',b, 'G , r=',r)
    sobs[[3,7]] *= b
    rms[[3,7]] *= np.sqrt(b)    
    if verbose ==1: print("{:4.6f}".format(time.time()-start),'FAKE_OBS: FIN')
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
    #print(sobs,rms,oindex)
    return sobs,rms,oindex



#####################HELPER ROUTINES FOR OBS_FAKE#######################
## These are not updated with numba or optimized. Just used for fake_obs
########################################################################
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