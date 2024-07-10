######################################################################


#  PROGRAM TO FIND BEST FIT PARAMETERS TO STOKES VECTORS FROM M1 CORONAL LINES
#
#  1. GET OBSERVATIONS, ELONGATION AND STOKES VECTORS
#  2. IDENTIFY AND READ THE SPECIFIC ELONGATION FILE
#  3. COMPUTE CHI^2 


######################################################################
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
import cle as cle
import sys
#
np.set_printoptions(precision=2,suppress=True)

lab=['I','Q','U','V']
lab+=lab
######################################################################
# NEEDED FUNCTIONS BEGIN
def readhdr(file):
#
    g=np.fromfile(file,dtype=np.float,sep=' ')
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
    nq=np.ones(nline,dtype=np.int32)
    return [g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
            bthetamin, bthetamax, bfield, nline, nq]

def obs_fake(dbdir,yobs,counts):
    #
    # INPUTS  yobs, counts
    # 
    # Normally one would read from an observation file.
    # Instead in this case read from files gridobs, obs.dat, obs.hdr
    # Reads files 'obs.hdr' and filen (e.g. OB0100.DAT) for observations
    # from filename, i.e. 0100/1000 + 1 solar radii = 1.1 solar radii
    # The idea is to return data along a certain elongation y from Sun
    # from another PCA calculation from CLE.
    # The data returned are just
    #  sobs    -    an array of size [4*nuse, nobs], normalized to maximum = 1.0
    #  rms     -    an array of size (4*nuse]  uncertainties
    #               nuse = number of lines to use = nline by default
    #
    # header file


    hfile=dbdir+'db.hdr'
    g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax, \
        bthetamin, bthetamax, bfield, nline, nq = readhdr(hfile)
    #
    #  set number of lines to use 
    #
    nuse=nline
    shap=[4*nq[0]*nline,ned*ngx*nbphi*nbtheta]
    s=np.empty(shap)
    #
    iy=round((yobs-1.0)*1000)
    so='%0*d' % (4, iy)
    fildb=dbdir+'DB'+so+'.DAT'

    f = open(fildb, 'rb')
    #
    d = np.fromfile(file=f, dtype=np.float32)
    shape=[ned,ngx,nbphi,nbtheta,nline,4]
    data=np.reshape(d,shape)
    #
    # get stokes data s and alignments
    # the indices i,j,k,l are picked at random
    #
    mx=ned*ngx*nbphi*nbtheta
    index= int(np.random.uniform(0.,1)*mx)
    index=97713 +5030010
    #index=97713 +3080010
    oindex=index
    i,j,k,l = cle.params(index,g)
    s0=data[i,j,k,l,0,:4]
    s1=data[i,j,k,l,1,:4]
    sobs=np.append(s0,s1)  # physical units

    sobs/=sobs[0]          # normalized
    sobs*=counts           # number of counts in each sobs

    #
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
    vary =  vars/sobs**2 + variance/sobs[0]**2
    vary *= (sobs/sobs[0])**2
    rms=sqrt(vary)
    #
    # Now add noise to the fake observations
    #
    sobs = sobs * ( 1. + rms*np.random.normal(0.,1,sobs.shape) )
    # one realization 
    #sobs *= (1. +rms*[0. , -1.73,  0.05, -0.13,  1.44,  1.29, -1.11, -1.02])
    sobs/=counts
    sobs[0]=1.
    #
    # Now fix the rms values for the denominator in the chi2 calculation
    # this is for each stokes value
    #
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
    print("\n")
    strc="{:0.1f}".format(log10(counts))

    print('Fake observations: log counts per state =',strc)
    print('S     ', [ "{:9.2e}".format(x) for x in sobs ])
    #print('rms   ',[ "{:9.2e}".format(x) for x in rms ])
    print('rms/s ',[ "{:9.2e}".format(x) for x in rms/sobs ])
    print("\n")
    
    return [nuse,sobs,rms,index]



######################################################################
# NEEDED FUNCTIONS END


######################################################################
# MAIN PROGRAM : some controlling parameters
#
iplot=1  # to make pdf files in directory fitted (can be time consuming)
#
#  the number of counts which correspond to the
#  largest intensity in each of the nobs observations
#  for example: counts = 1.e8 gives rms fluctuations 10^-4 * max(intensity)
#
xtra='ls'
xtra=''
xtra='a'
chi2max=2.
nstokes=4
dbdir='dbcle/fe13'+xtra+'/'
outdir='results/fe13'+xtra+'/'
yobs=1.18
iy=round((yobs-1.0)*1000)
so='%0*d' % (4, iy)
fildb=dbdir+'DB'+so+'.DAT'

######################################################################
# LOOP OVER COUNTS

lmin=5.
lmax=8.0
l = np.arange(lmin,lmax,0.3)
fr=l*0.
kount=0
for lcounts in l:
    counts=10.**lcounts
    strc="{:0.1f}".format(lcounts)
    
    #
    #######################################################################
    #
    # 1. Read Observations: 
    #   
    #######################################################################
    #

    # sobs is the observed Stokes array (nline*4)
    nuse,sobs,rms,oindex = obs_fake(dbdir,yobs,counts)
    #
    #
    #######################################################################
    #
    # 2. Read computed database for selected y-position 
    #
    #######################################################################
    #
    #
    g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
    bthetamin, bthetamax, bfield, nline, nq = readhdr(dbdir+'db.hdr')
    #
    #  set number of lines to use 
    #
    nuse=nline
    #
    shap=[4*nq[0]*nline,ned*ngx*nbphi*nbtheta]
    s=np.empty(shap)
    #
    file = open(fildb, 'rb')
    #
    dbd = np.fromfile(file=file, dtype=np.float32)
    shape=[ned*ngx*nbphi*nbtheta,nline,nstokes]
    data=np.reshape(dbd,shape)
    #
    # get stokes data s and alignments
    # 
    #align=data[:,:,4]
    s0=data[:,0,:4]
    s1=data[:,1,:4]
    #
    # shape of the stokes data sdb ("s database")
    #  is (number of calcs, 4*nuse)
    #
    sdb=np.append(s0,s1,axis=1)                                                                     
    #
    # Always normalize to the strongest of the raw stokes parameter
    # 
    norm=sdb[:,0]
    sdb = (sdb.T/norm.T).T
    #
    ######################################################################
    #
    #  Now find the matched data
    #
    ######################################################################
    #
    start=time.time()
    #
    filn=outdir+so+str(oindex)+"C"+strc+".txt"
    os.system("rm "+filn)
    Fobj = open(filn,"a")
    nn=16  # upper limit of number of nearest fits to permit, 
    #
    nc=sobs.shape[0]
    dif=( (sdb-sobs) / rms )[:,1:]
    sdif=dif*dif
    #
    # nf = number of degree of freedom in model: ne, x, bphi, btheta
    ndata = nline*4 -1
    ndof=4
    denom=nline*4-1 - ndof
    print("Ndata = ",ndata," N degrees of freedom = ",ndof)
    #
    chi2 = np.sum( sdif, axis=1)/denom
    #
    #  fast:
    #     a = np.array([9, 4, 4, 3, 3, 9, 0, 4, 6, 0])
    #
    #nnc=32
    #ind = np.argpartition(chi2, nn)[:nnc]
    #asrt= ind[np.argsort(chi2[ind])]
    #
    #
    asrt=np.argsort(chi2)

    ix=asrt[0]
    firstchi2=chi2[ix]
    print("best s",sdb[ix,:], '\nbest chi2 = ', firstchi2)

    if(iplot > 0): fig, (ax1, ax2) = plt.subplots(2, 1)

    grid, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
        bthetamin, bthetamax, bfield, nline, nq = readhdr(dbdir+'db.hdr')

    s2d=180./np.pi


    Fobj.write("      Observation    chi2     Ne        y       x    Bphi    Bthet   |B|")
    Fobj.write("\n")

    np.savetxt(Fobj,[oindex,oindex], fmt='%7i', newline=" ")
    np.savetxt(Fobj,[0.], fmt='%+.2e', newline=" ")
    led,gy,gx,bphi,btheta,b = cle.physa(oindex,yobs,grid)
    ph=cle.physa(oindex,yobs,grid)
    print('Physcs',ph) # [ "{:+0.2e}".format(x) for x in ph ])
    np.savetxt(Fobj,[cle.physa(oindex,yobs,grid)],fmt='%+7.2f', newline=" ")
    Fobj.write("\n")
    Fobj.write("\n")
    si=0
    Fobj.write("      O       C      chi2     Ne        y       x    Bphi    Bthet   |B|   POS AZ  LOS THETA")
    Fobj.write("\n")
    
    while (si <= nn):
        ix=asrt[si]
        dd=chi2[ix]
        ddi= "{:0.2e}".format(dd)
        if(dd < chi2max):
            np.savetxt(Fobj,[oindex,ix], fmt='%7i', newline=" ")
            np.savetxt(Fobj,[chi2[ix]], fmt='%+.2e', newline=" ")
            led,gy,gx,bphi,btheta,b = cle.physa(ix,yobs,grid)
            np.savetxt(Fobj,[cle.physa(ix,yobs,grid)],fmt='%+7.2f', newline=" ")
            los=np.arccos(np.sin(btheta)*np.cos(bphi))*s2d
            pos=np.arctan2(np.sin(btheta)*np.sin(bphi),np.cos(btheta))*s2d
            np.savetxt(Fobj, [pos,los],fmt='%+7.2f', newline=" ")
            Fobj.write("\n")
        if(iplot > 0):
            text="fit="+str(ix)+" $\chi^2$="+ddi
            if(si ==0):
                ax1.plot(sobs,'.')
                text="Obs sp="+str(oindex)+" fit="+str(ix)+" $\chi^2$="+ddi
            if(si <= 2):
                    ax1.plot(sdb[ix,:],label=text)
                    ax1.errorbar(np.arange(nc),sobs,yerr=rms,fmt='bo')
                    ax1.legend(fontsize=7)
                    ax1.set_xticks(np.arange(8))
                    ax1.set_xticklabels(lab)
            if(si ==0): ax2.plot(dif[ix,:],'.-',label="(O-N)/sigma")
            ax2.legend(fontsize=7)
            ax2.set_xticks(np.arange(8))
            ax2.set_xticklabels(lab)
            si=si+1
            
    if(iplot > 0):
        filen=outdir+'fit'+so+str(oindex)+"C"+strc
        print("Saving ", filen+".pdf")
        savefig(filen+".pdf")
        plt.close()
        #
    Fobj.close()
    dt= "{:3.2f}".format(time.time()-start)
    print(dt,' seconds for chi2')
    start=time.time()
                
    ######################################################################
    # Now here we examine the statistics of these apparent "fits"
    ######################################################################
    ok=np.argwhere(chi2 < chi2max)
    nok=ok.size
    ntot=chi2.size
    fr[kount]=np.float32(nok)/chi2.size
    print(fr[kount],' <--- fraction of ok fits, log counts=',strc,'nok=',nok,' ntot=',ntot)
    #
    d=np.zeros([6,nok])
    pobs= cle.physa(ix,yobs,grid)
    print("assigning physical parameters..")
    for n in range(0,nok):
        ix=(ok[n])[0]
        a,b,c,dd,e,f = cle.physa(ix,yobs,grid)
        d[0,n] = a
        d[1,n] = b
        d[2,n] = c
        d[3,n] = dd
        d[4,n] = e
        #ratio=sobs[3]/sdb[ix,3]  # ratio of Vobs/ Vdb gives B
        d[5,n] = f 

    ttext=['log$_{10}N_e$','y (pos)','x (LOS)','$\Phi_B$','$\Theta_B$','B']
    fig, axs = plt.subplots(2,3,tight_layout=True)
    title=strc
    #
    plt.rcParams.update({'font.size': 7})

    pobs=cle.physa(oindex,yobs,g) 
    scounts='log counts = '+strc+' '

    f=np.array([1,1,1,1,1,1])

    alpha=.7
    for i in range(0,6):
        j=int(i/3)
        k=int(i-j*3)
        axs[j,k].hist(d[i,:]*f[i],label=text[i],alpha=alpha,density=0,bins=50,align='left')

        if i > 0: scounts=''
        axs[j,k].set_xlabel(ttext[i])
        axs[j,k].set_title(scounts)
        
        yr = axs[j,k].get_ylim()
        ij=0
        xr=[0.,0.]+pobs[i]*f[i]
        axs[j,k].plot(xr, yr)

#
    print('Physcs',pobs*f) # [ "{:+0.2e}".format(x) for x in ph ])

    filen=outdir+so+'hist'+xtra+str(oindex)+"C"+strc
    print("Saving ", filen+".pdf")
    savefig(filen+".pdf")
    plt.close()
#    os.system("open "+filen+".pdf")
    dt= "{:3.2f}".format(time.time()-start)
    print(dt,' seconds for search')
    kount+=1

print('GASP!')    
#
plt.plot(l,log10(fr*ntot),'-o')
plt.title(str(oindex) + " " +xtra+ str(pobs))
plt.ylabel('Log # acceptable solutions ')
plt.xlabel('Log counts')
filen=outdir+so+'stats'+xtra+str(oindex)
print("Saving ", filen+".pdf")
savefig(filen+".pdf")
plt.close()
os.system("open "+filen+".pdf")
#


