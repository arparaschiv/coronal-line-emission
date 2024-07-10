###Code utils for findmany_phil from phil to seek database entries and solve gfactors for the dima degeneracy issue
##download from gmail 20201206
import numpy as np

def electron_density(r):
  baumbach = 1.e8*(0.036/r**1.5 + 1.55/r**6.)
  hscale=   7.18401074e-02  # 50 Mm
  n0=3.e8*np.exp(- (r-1.)/hscale) + baumbach
  return n0

def phys(index,gy,grd):
  #
  i,j,k,l,m = params(index,grd)
  p=physa(index,gy,grd)
  bphi=p[3]
  btheta=p[4]
  b=p[5]
  bx=b*np.sin(btheta)*np.cos(bphi)
  by=b*np.sin(btheta)*np.sin(bphi)
  bz=b*np.cos(btheta)
  out=np.array([p[0],p[1],p[2],bx,by,bz])
  return out

def physa(index,gy,grd):
  #
  i,j,k,l,m = params(index,grd)
  #
  grd,ned, ngx, nbphi, nbtheta, nb, xed, gxmin,gxmax, bphimin, bphimax,\
    bthetamin, bthetamax, bfield, nline,wavel  = parsehdr(grd)

  gx = gxmin + j*(gxmax-gxmin)/(ngx-1)
  #
  bphi  = bphimin + k*(bphimax-bphimin)/(nbphi-1)
  btheta  = bthetamin + l*(bthetamax-bthetamin)/(nbtheta-1)
  #
  # log of Ne
  #
  ed = xed[i]* electron_density(np.sqrt(gy*gy+gx*gx))
  #
  # |B|
  #
  out=np.array([np.log10(ed),gy,gx,bphi,btheta,bfield[m]])
  return out

def params(index,grd):
  # for index, get i,j,k,l,m   indices


  ned = int(grd[0])
  ngx = int(grd[1])
  nbphi = int(grd[2])
  nbtheta = int(grd[3])
  nb = int(grd[4])
  
  n5=( ngx*nbphi*nbtheta*nb)
  n4=(     nbphi*nbtheta*nb)
  n3=(           nbtheta*nb)
  n2=(                   nb)
  #
  i = index   //        n5
  
  j = index   -     i * n5
  j = j      //         n4
  
  k = index   -     i * n5 - j * n4
  k = k      //         n3
  
  l=  index   -     i * n5 - j * n4 - k * n3
  l= l       //         n2

  m=  index   -     i * n5 - j * n4 - k * n3 - l*n2
  #print('params ',i,j,k,l,m)
  return np.array([i,j,k,l,m])

   
######################################################################
# NEEDED FUNCTIONS BEGIN
def readhdr(file):
#
# just read the header of the DB*DAT data cubes
#
  g=np.fromfile(file,dtype=np.float,sep=' ')
  return g

def parsehdr(g):
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
  #
  # log of Ne
  #
  emin=g[5]
  emax=g[6]
  xed=lgrid(emin,emax,ned)
  #
  # log of |B|
  #
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

    ff=10.** ((mx-mn)/(n-1))
    x=np.empty(n)
    x[0]=10.**mn
    for k in np.arange(1,n): x[k]=x[k-1]*ff
    return x