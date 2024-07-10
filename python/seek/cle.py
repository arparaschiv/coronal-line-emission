import numpy as np

def electron_density(r):
  baumbach = 1.e8*(0.036/r**1.5 + 1.55/r**6.)
  hscale=   7.18401074e-02  # 50 Mm
  n0=3.e8*np.exp(- (r-1.)/hscale) + baumbach
  return n0

def phys(index,gy,grd):
  #
  i,j,k,l = params(index,grd)
  #
  ned=int(grd[0])
  xed=np.empty(ned)
  for kl in np.arange(0,ned): xed[kl]=grd[5+kl]
  #
  n=4+ned
  gx = grd[n] + (grd[n+1]-grd[n])*j/grd[1]
  #
  ed = xed[i]* electron_density(np.sqrt(gy*gy+gx*gx))
  #
  n+=2
  bphi  = grd[n] + (grd[n+1]-grd[n])*l/grd[2]
  n+=2
  btheta= grd[n] + (grd[n+1]-grd[n])*m/grd[3]
  bx=np.sin(btheta)*np.cos(bphi)
  by=np.sin(btheta)*np.sin(bphi)
  bz=np.cos(btheta)
  out=np.array([ed,gy,gx,bx,by,bz])
  return out




def physa(index,gy,grd):
  #
  i,j,k,l = params(index,grd)
  #
  ned=int(grd[0])
  xed=np.empty(ned)
  for kl in np.arange(0,ned): xed[kl]=grd[5+kl]
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

def pcprep(d):
  sh=d.shape
  v=np.average(np.abs(d))
  if len(sh) > 1: v = np.average(np.abs(d),axis=0)
  ret=(d/v)
  
  return v,ret
  

   
