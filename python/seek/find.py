import os
import numpy as np
from scipy.io import FortranFile
from pylab import *
from astropy.io import ascii
from matplotlib import pyplot as plt
from scipy import stats
import matplotlib.ticker as mtick
import matplotlib.pylab as pylab
import warnings
import pickle
import time

import h5py
import difflib

## to import the custom pa library in ipython
import sys
sys.path.insert(1, '/home/alin/Documents/physics_prog/cle/python/seek')
import cle as cle

sys.path.insert(1, '/home/alin/Documents/physics_prog/python')
from alin_misc import fullprint


#
np.set_printoptions(precision=3,suppress=True)
np.set_printoptions(precision=2,suppress=False)
lab=['I','Q','U','V']
lab+=lab
lab+=lab
lab+=lab

def readhdr(file):
#
    g=np.fromfile(file,dtype=float,sep=' ')
    ned=int(g[0]) 
    ngx=int(g[1])
    nbphi=int(g[2])
    nbtheta=int(g[3])
    xed=np.empty(ned)
    for k in arange(0,ned): xed[k]=g[4+k]
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

def obs_fake(dirn,filen,yobs,counts,index):
    #
    # Normally one would read from an observation file.
    # Instead in this case read from files gridobs, obs.dat, obs.hdr
    # Reads files 'obs.hdr' and filen (e.g. OB0100.DAT) for observations
    # from filename, i.e. 0100/1000 + 1 solar radii = 1.1 solar radii
    # The idea is to return data along a certain elongation y from Sun
    # from another PCA calculation from CLE.
    # The data returned are just
    #  sobs    -    an array of size [4*nuse, nobs], normalized to maximum = 1.0
    #  rms     -    an array of size (4*nuse]  scaled uncertainties
    #               nuse = number of lines to use = nline by default
    #
    g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
        bthetamin, bthetamax, bfield, nline, nq = readhdr(dirn+'obs.hdr')
    #
    #  set number of lines to use 
    #
    nuse=nline
    #
    shap=[4*nq[0]*nline,ned*ngx*nbphi*nbtheta]
    s=np.empty(shap)
    #sigm=np.empty(nline,shap[1])
    #
    f = open(dirn+filen, 'rb')
    #
    d = np.fromfile(file=f, dtype=np.float32)
    shape=[ned,ngx,nbphi,nbtheta,nline,5]
    data=np.reshape(d,shape)
    #
    # get stokes data s and alignments
    # the indices 1,3,7,5   are entirely random
    #
    mx=ned*ngx*nbphi*nbtheta
    #change index to specific if inputted
    # if sys.argv[5]==-1:
    #     print('Using random index')
    #     index=int(np.random.uniform(0.,1)*mx)
    # elif 0 <= sys.argv[5] < mx:
    #     index=sys.argv[5]
    # else:
    #     print('input index is out of bounds! Using random index')
    #     index=int(np.random.uniform(0.,1)*mx)
    if index ==-1:
        print('Using random index')
        index=int(np.random.uniform(0.,1)*mx)
    elif 0 <= index < mx:
        print("using index:", index)
    else:
        print('input index is out of bounds! Using random index')
        index=int(np.random.uniform(0.,1)*mx)

    i,j,k,l = cle.params(index,g)
    s0=data[i,j,k,l,0,:4]
    s1=data[i,j,k,l,1,:4]
    #
    # here are the physical data
    #
    sobs=np.append(s0,s1)
    #
    # This index is an encoding of the [i,j,k,l,:4] stokes array
    # it can be converted back to j,j,k,l  using cle.params(index, grid)
    #
    #index= ngx*nbphi*nbtheta*1 + nbphi*nbtheta*3 +nbtheta*7 + 5
    print('Obs.  s =', s0/np.max(s0), '\n      s =',s1/np.max(s0), '\nPhysics =',cle.physa(index,yobs,g))
    #
    # normalize sobs to individual maxima
    #
    #for i in range(0,2):
    #    mx=np.max(sobs[i*4:(i+1)*4])
    #    sobs[i*4:(i+1)*4] = sobs[i*4:(i+1)*4]/mx
    mx=np.max(sobs)
    sobs = sobs/mx
    #
    #  here we add random fluctuations to sobs determined by counts above
    #
    err=1./sqrt(counts)
    #
    # Now, re-scale qu and v by fq and fv to make s roughly (1,1,1,1)
    # accordingly, the noise magnitudes must be scaled by the same factors:
    #  Note: this tends to change only the chi2 values by a single multiplier
    #  thereby not affecting the ordering of the chi2
    fqu=10.
    fv=100.
    f0=np.array([1,fqu,fqu,fv])
    for ii in range(0,nuse): f=np.append(f0,f0)
    
    sobs=f * sobs * ( 1. + err*np.random.normal(0.,1,f.shape) )
    #
    # Now fix the rms values for the denominator in the chi2 calculation
    # this is for each stokes value
    #
    rms=f*err
    #
    ######################################################################
    # finished observations, which consist only of
    #
    #  nuse, sobs, rms, and the factor f that scales QUV wrt I
    #
    # yobs is the elongation in Rsun of the observed y-position 
    # index is the index of the calculation used to make fake obs
    #       (not needed for real observations)
    #
    ######################################################################

    return [nuse,sobs,rms,f,index,mx]

########################################################################
# Small script for selecting the correct database entry.
def get_db_entry(yobs,dbdir):
    
    ## get the distance in the db format
    iy=round((yobs-1.0)*1000)
    so='%0*d' % (4, iy)
    fil1=int(so)

    ##get the database entries
    dblist=os.listdir(dbdir)  
    norm=np.zeros(len(dblist))
    for i in range(0,len(dblist)):
        if dblist[i].find(".DAT") == -1:
            norm[i]=1000    ## if file in list is not database generated, ignore it from the norm.
        else:
           norm[i]=np.absolute(fil1- int(dblist[i].split('DB')[1].split(".DAT")[0])) ## compute the norm between the observed distance and db file distance

    ##select the closest match to the inputed radius
    db_file =dblist[np.where(norm == norm.min())[0][0]]    
    
    return db_file


######################################################################
# some controlling parameters
#
iplot=1  # to make pdf files in directory out1 (can be time consuming)
#
#  the number of counts which correspond to the
#  largest intensity in each of the nobs observations
#  for example: counts = 1.e8 gives rms fluctuations 10^-4 * max(intensity)
#
counts=1.e9
#
######################################################################


#######################################################################
#
# 1. Read Observations: 
#   
#######################################################################
#
# yobs=1.427
# iy=round((yobs-1.0)*1000)
# so='%0*d' % (4, iy)
# dir1='./obs/fe13/'
# fil='DB'+so+'.DAT'   ##OB in current dir
# # sobs is the observed Stokes array (nline*4)
# # index of the cle output to use as obs. use -1 to randomize.
# nuse,sobs,rms,f,oindex,mx = obs_fake(dir1,fil,yobs,counts,6397)
# sobs=sobs[0:4]
# f=f[0:4]
# rms=rms[0:4]

nuse=3 ## number of lines to be used
yobs=1.228

dir1='/home/alin/Documents/physics_prog/cle/test_db_obs/fe13/'
dir2='/home/alin/Documents/physics_prog/cle/test_db_obs/si09/'
dir3='/home/alin/Documents/physics_prog/cle/test_db_obs/si10/'

##get the closest matching database entry to yobs
fil=get_db_entry(yobs,dir1)


## sobs is the observed Stokes array (nline*4)
nuse1,sobs1,rms1,f1,oindex1,mx1 = obs_fake(dir1,fil,yobs,counts,6959)

if nuse > 2:
    nuse2,sobs2,rms2,f2,oindex2,mx2 = obs_fake(dir2,fil,yobs,counts,6959)
    if mx1 > mx2:
        sobs2=sobs2*mx2/mx1
    else:
        sobs1=sobs1*mx1/mx2

    sobs=np.append(sobs1,sobs2[0:4])
    f=np.append(f1,f2[0:4])
    rms=np.append(rms1,rms2[0:4])
else:
    sobs=sobs1[:nuse*4]
    f=f1[:nuse*4]   ##here you need to compute the correct outputs for the number of lines
    rms=rms1[:nuse*4]

#######################################################################
#
# 1. Read Observations - MURAM
#   
#######################################################################
# #Read the muram data
# f1a = h5py.File('/home/alin/Desktop/cle-invert/hh2-muram/Fe13_10747_hires_spectra_050000.hdf5','r')
# f1b = h5py.File('/home/alin/Desktop/cle-invert/hh2-muram/Fe13_10747_hires_tot_IQU_050000.hdf5','r')
# f2a = h5py.File('/home/alin/Desktop/cle-invert/hh2-muram/Fe13_10798_hires_spectra_050000.hdf5','r')
# f2b = h5py.File('/home/alin/Desktop/cle-invert/hh2-muram/Fe13_10798_hires_tot_IQU_050000.hdf5','r')
# ## get the scales for the wavelength and spatial axes
# ## they should be the same for the two data sets 
# wvvec1 = f1a['stokes'].dims[1][0][:]
# wvvec2 = f2a['stokes'].dims[1][0][:]
# yvec = f1a['stokes'].dims[2][0][:]
# xvec = f1a['stokes'].dims[3][0][:]
# ##compute the total integrated signals.
# f1aa=np.sum(f1a['stokes'],axis=1)
# f2aa=np.sum(f2a['stokes'],axis=1)
# f1aaa=np.zeros((101,1024,1024))
# f1aaa=np.abs(f1a['stokes'][3,:,:,:])
# f2aaa=np.zeros((101,1024,1024))
# f2aaa=np.abs(f2a['stokes'][3,:,:,:])
# f1aa[3,:,:]=np.sum(f1aaa,axis=0)
# f2aa[3,:,:]=np.sum(f2aaa,axis=0)
# # sobs is the observed Stokes array (nline*4)
# nuse=2

# ## just take a random pixel along the central x.
# xx=int(np.random.uniform(0.,1024))
# yy=512#int(np.random.uniform(0.,1024))

# ##form the array
# sobs=[f1aa[0,xx,yy],f1aa[1,xx,yy],f1aa[2,xx,yy],f1aa[3,xx,yy],f2aa[0,xx,yy],f2aa[1,xx,yy],f2aa[2,xx,yy],f2aa[3,xx,yy]]
# sobs = sobs/max(sobs)
# #here we add random fluctuations to sobs determined by counts above (see subroutine above)
# err=1./sqrt(counts)
# fqu=10.
# fv=100.
# f0=np.array([1,fqu,fqu,fv])
# for ii in range(0,nuse): f=np.append(f0,f0)
    
# sobs=f * sobs * ( 1. + err*np.random.normal(0.,1,f.shape) )
# rms=f*err

# # Compute the radial distance from teh sun
# xvec[xx]/696.34+1
# yobs=xvec[xx]/696.34+1
# iy=round((yobs-1.0)*1000)
# so='%0*d' % (4, iy)

# print(' Muram Obs.  s =', sobs, '\n' )
#
#
#######################################################################
#
# 2. Read computed database for selected y-position 
#
#######################################################################
#
#

## the fe 13 atom calculations contain 2 lines, 
## the si9 and si10 atoms contain 1 line each.

##nuse=1 #nline defined earlier
d_dir='/home/alin/Documents/physics_prog/cle/test_db_large'
#d_dir='/home/alin/Documents/physics_prog/cle/test_db_obs'
#d_dir='/home/alin/Documents/physics_prog/cle/test_db_small'
g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
    bthetamin, bthetamax, bfield, nline, nq = readhdr(d_dir+'/fe13/db.hdr')
#
#  set number of lines to use 
#
#
#shap=[4*nq[0]*nline,ned*ngx*nbphi*nbtheta]
#s=np.empty(shap)
#

filen=get_db_entry(yobs,d_dir+'/fe13') ## get the closest matching database entry

file = open(d_dir+'/fe13/'+filen, 'rb')
#
d = np.fromfile(file=file, dtype=np.float32)
shape=[ned*ngx*nbphi*nbtheta,nline,5]
data=np.reshape(d,shape)
#
# get stokes data s and alignments
#
#align=data[:,:,4]


# shape of the stokes data sdb ("s database")
#  is (number of calcs, 4*nuse)     

s0=data[:,0,:4]

if nuse >1:
    s1=data[:,1,:4]

if nuse == 2:
    sdb=np.append(s0,s1,axis=1)
elif nuse == 1:
   sdb=s0 
                                                                


######if a 3rd line exists###########
if nuse >2:
    g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
        bthetamin, bthetamax, bfield, nline, nq = readhdr(d_dir+'/si09/db.hdr')
    #

#    shap=[4*nq[0]*nline,ned*ngx*nbphi*nbtheta]
#    s=np.empty(shap)
    #
    filen=get_db_entry(yobs,d_dir+'/si09') ## get the closest matching database entry

    file = open(d_dir+'/si09/'+filen, 'rb')
    #
    d = np.fromfile(file=file, dtype=np.float32)
    shape=[ned*ngx*nbphi*nbtheta,nline,5]
    data=np.reshape(d,shape)
    #
    # get stokes data s and alignments
    #
    #align=data[:,:,4]
    s2=data[:,0,:4]
    if nuse == 3:
        sdb_temp=sdb=np.append(s0,s1,axis=1)
        sdb=np.append(sdb_temp,s2,axis=1)                                                                     

######if a 4th line exists###########
if nuse == 4:
    g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
        bthetamin, bthetamax, bfield, nline, nq = readhdr(d_dir+'/si10/db.hdr')
    #

#    shap=[4*nq[0]*nline,ned*ngx*nbphi*nbtheta]
#    s=np.empty(shap)
    #
    filen=get_db_entry(yobs,d_dir+'/si10') ## get the closest matching database entry

    file = open(d_dir+'/si10/'+filen, 'rb')
    #
    d = np.fromfile(file=file, dtype=np.float32)
    shape=[ned*ngx*nbphi*nbtheta,nline,5]
    data=np.reshape(d,shape)
    #
    # get stokes data s and alignments
    #
    #align=data[:,:,4]
    s3=data[:,0,:4]
    sdb=np.append(s0,[s1,s2,s3],axis=1)                                                                     

# Always normalize to the strongest of the raw stokes parameter
# 
#for i in range(0,nuse):
#    mx=np.max(sdb[:,4*i:4*(i+1)],axis=1)
#    sdb[:,4*i:4*(i+1)] = (sdb[:,4*i:4*(i+1)].T/mx.T).T

mx=np.max(sdb,axis=1)
sdb = (sdb.T/mx.T).T
#
# apply form factor f to database as well as obs

sdb = sdb*f[:4*nuse]

#
#
######################################################################
#
#  Now find the matched data
#
######################################################################
#
start=time.time()
#
outdir=d_dir+'/fitted/'
if (sdb.shape[1] <=8 ):
    oindex1=oindex1  ##check why?
os.system("rm "+outdir+str(oindex1)+".txt")
Fobj = open(outdir+str(oindex1)+".txt","a")
Fobj.write("      O       C      chi2        Ne         y         x     phi_B     theta_B   |B|")
nn=15  # upper limit of number of nearest fits to permit, 
#
nc=sobs.shape[0]
x=sobs
dif=(sdb-x)
sdif=dif*dif#/rms/rms
#sdif=dif
chi2 = np.sum( sdif, axis=1)/(nc-1)
#
#  Here is a minimal search algorithm- a mere sort of chi2
#
#  fast:
#     a = np.array([9, 4, 4, 3, 3, 9, 0, 4, 6, 0])
ind = np.argpartition(chi2, 16)[:16]
asrt= ind[np.argsort(chi2[ind])]
#
#
#asrt=np.argsort(chi2)
#print(asrt-asrt1[:16])

firstchi2=chi2[asrt[0]]
si=0
ix=asrt[si]
print("best s",sdb[ix,:])
Fobj.write("\n")

if(iplot > 0): fig, (ax1, ax2) = plt.subplots(2, 1)

grid, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
    bthetamin, bthetamax, bfield, nline, nq = readhdr(d_dir+'/fe13/db.hdr')

while ((chi2[ix]/firstchi2) < 1.2 and (si <= nn)):
    ix=asrt[si]
    dd=chi2[ix]
    ddi= "{:0.2e}".format(dd)
    if(dd/firstchi2 < 1.2):
        np.savetxt(Fobj,[oindex1,ix], fmt='%7i', newline=" ")
        np.savetxt(Fobj,[chi2[ix]], fmt='%+.2e', newline=" ")
        np.savetxt(Fobj,[cle.physa(ix,yobs,grid)],fmt='%+9.2f', newline=" ")
        Fobj.write("\n")

    if(iplot > 0):
        text="fit="+str(ix)+" $\chi^2$="+ddi
        if(si ==0):
            #ax1.ticklabel_format(style='sci',axis='y',limits=(0,1))
            #ax1.set_yscale('log')
            #ax1.set_ylim(-1,1)
            text="Obs sp="+str(oindex1)#+" fit="+str(ix)+" $\chi^2$="+ddi
            x1=np.copy(x)
            if (sdb.shape[1] >8 ): 
                x1[4:8]=x1[4:8]*mx1/(x1[4]*mx1) 
                x1[8:12]=x1[8:12]*mx2/(x1[8]*mx2)
            ax1.plot(x1,'.',label=text)
        if(si == 1):   
            text= "fit="+str(ix)+" $\chi^2$="+ddi
            sdb1=np.copy(sdb[ix,:]) 
            if (sdb.shape[1] >8 ):        
                sdb1[4:8]=sdb1[4:8]*mx1/(sdb1[4]*mx1) 
                sdb1[8:12]=sdb1[8:12]*mx2/(sdb1[8]*mx2)
            ax1.plot(sdb1,'.',label=text)
            ax1.errorbar(np.arange(nc),x,yerr=rms,fmt='none')
            ax1.legend(fontsize=7)
            ax1.set_xticks(np.arange(sdb.shape[1]))
            ax1.set_xticklabels(lab[:sdb.shape[1]])
            ax1.axvline( x=3.5, ymin=-0.25, ymax=1)
            ax1.text(x=0.5,y=0.9,s="Fe ${XIII}$ 1074.7nm",size=6.5)
            ax1.text(x=4.5,y=0.9,s="Fe ${XIII}$ 1079.8nm",size=6.5)
            if (sdb.shape[1] >8 ):
                ax1.axvline( x=7.5, ymin=-0.25, ymax=1)
                ax1.text(x=8.5,y=0.9,s="Si ${IX}$ 3934.3nm",size=6.5)
        if(si ==0): 
            ax2.plot(dif[ix,:]/rms,'g.-',label="(O-N)/$\sigma$")
            ax2.legend(fontsize=7)
            ax2.set_xticks(np.arange(sdb.shape[1]))
            ax2.set_xticklabels(lab)
            ax2.axvline( x=3.5, ymin=-0.25, ymax=1)
            ax2.text(x=0.5,y=0.9,s="Fe ${XIII}$ 1074.7nm",size=6.5)
            ax2.text(x=4.5,y=0.9,s="Fe ${XIII}$ 1079.8nm",size=6.5)
            if (sdb.shape[1] >8 ):
                ax2.axvline( x=7.5, ymin=-0.25, ymax=1)
                ax2.text(x=8.5,y=0.9,s="Si ${IX}$ 3934.3nm",size=6.5)
    si=si+1
if(iplot > 0):
    filen=outdir+'fit'+str(oindex1)
    print("Saving ", filen+".pdf")
    savefig(filen+".pdf")
    plt.close()
#
Fobj.close()
dt= "{:3.2f}".format(time.time()-start)
print(dt,' seconds')
