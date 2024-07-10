###############################################################################
### Parallel implementation of utilities for CLE DBINVERT coronal inversion ###
###############################################################################

#Contact: Alin Paraschiv, National Solar Observatory
#         arparaschiv@nso.edu
#         paraschiv.alinrazvan+cledb@gmail.com


### Tested Requirements for DBINVERT: ################################################
### CLEDB database as configured in V2.0.2; 
# python                    3.7.10 
# scipy                     1.6.2
# numpy                     1.20.1 
# numba                     0.53.1 
# llvmlite                  0.36.0 (as numba dependency probably soed not need separate installing)

### Other used packages of secondary importance that will probably not become code-breaking 
# pylab    ## Check if still needed
# numexpr
# glob

### Auxiliary, only used in verbose mode
# time
# os
# sys

### Crutches for data handling during testing, not utilized by the inversion
# importlib ## to recompile the auxiliary file before each run
# pickle    ## convenient python/numpy data storage and loading
###


### NOTES: ##############################################################
## global variables can't really be used with numba parallel python functions unless you aim to keep everything constant.
### main chi2 fitting and matching algorithm

## Functions in the time module are not understood/suported by numba. enabling them forces the compiler to go back to object mode and throw warnings. 
## Full Parallelization will generally not be possible while time functions are enabled.

## Array indexing by lists and advanced indexing are generally not supported by numba.
## Do not be fooled by the np. calls. All np. functions used are rewritten by numba.
## Revert to simple slicing, and array initialization by tuples.

## Generally, if a numba non-python  (njit, or jit without explicit object mode flag) 
## calls another function, the second should also be non-python compatible 

#########################################################################
### Needed imports ######################################################
#########################################################################

#
# Needed libraries
#
import numpy as np

from numba import jit,njit, prange
from numba.typed import List  ## numba list is needed ad standard reflected python lists will be deprecated in numba


from scipy.io import FortranFile
from pylab import *
#from matplotlib import pyplot as plt
#from scipy import stats
import time
import glob
import os
#from random import seed
#from random import randint
import sys
import numexpr as ne

#########################################################################
### Main solvers ########################################################
#########################################################################

@njit#(forceobj=True)
def clematch_par(sobs,yobs,database,dbhdr,rms,outdir,maxchisq,nsearch,reduction,verbose):
## main solver for the geometry and magnetic field strength 
## Returns matched database index, and double IQUV vector,  and chi^2 fitting residual
## Returns matched observation physics:
## Y(radial) position, X(LOS) position, Bphi and Btheta angles in L.V.S coordinates; and their transform to cartesian Bx, By, Bx in LOS geometry, and B field strength. 
    #start=time.time()
    #######################################################################
    ## unpack dbcgrid parameters from the database accompanying db.hdr file 
    ##
    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
        bthetamin, bthetamax, nline, wavel  = dbhdr 

    ####NO OBSERVATION --> NO RUN! ###########    
    ## e.g. a pixel inside the solar disk, an invalid pixel, etc.
    ##returns an array-like 0 vector of the dimensions of the requested output.
    ## the index is set to -1 to warn about hte missing entry
    if np.isnan(sobs).all():
        #print("exit1")
        nullout=np.zeros((nsearch,11),dtype=np.float64)
        nullout[:,0]=np.full(nsearch,-1)
        if verbose >=2: print("Warning: No observation in voxel!")
        #return np.zeros((nsearch),dtype=np.int_),np.zeros((nsearch,8),dtype=np.float64),np.zeros((nsearch),dtype=np.float64),np.zeros((nsearch,9),dtype=np.float64)    
        return nullout,np.zeros((nsearch,8),dtype=np.float64)

    #######################################################################
    ##  Read the database in full or reduced mode  
    #start0=time.time()
    if reduction == 1:
        database,outindex=get_subset_par(sobs,dbhdr,database,nsearch,verbose) 
        #if verbose >=1:print("clematch with reduced DB:",sdb.shape,outindex.shape)
        ## sdb is the reduced database in this cass
        ## outindex its its corresponding index, that is different from the full index.
        #print("{:4.6f}".format(time.time()-start),' SECONDS FOR SUBSET')
    else:
        #ARP:append is much slower than reshape! will have to see if its possible to switch to a reshape method
        #this will be revisited for single ion DBs.  
        ##start=time.time()                     
        #sdb=np.append(database[:,0,:4],database[:,1,:4],axis=1)
        #dt= "{:4.6f}".format(time.time()-start)
        #if verbose ==1:print(dt,' SECONDS FOR SDB APPEND')
        #start=time.time()
        #sdb=np.reshape(database,(database.shape[0],8))
        #dt= "{:4.6f}".format(time.time()-start)
        #if verbose ==1:print(dt,' SECONDS FOR SDB RESHAPE')
        outindex=np.empty((nbphi,nsearch)) 
        #sdb=database
        database=np.reshape(database,(ned*ngx*nbphi*nbtheta,8))
        #if verbose >=1:print("clematch with full DB:",sdb.shape,outindex.shape)                     
        ## In this branch, outindex is never used
        ## numba does not properly compile if it does not have this array defined because it is named in the code!
        #if verbose >=2:
        #    if (sdb.shape[0] > 1000000):
        #        print("Warning: Full database match with a large calculation will be significantly slower with no proven benefits")

    #######################################################################
    ##  Geometric solution is based on a reduced chi^2 measure.   

    ## Number of observables for geometry solver
    ndata = 4*2 - 2 -1  ## -2 comes from not using the two V components, as V is independent of the observation geometry.
    # ndof = number of degree of freedom in model: ne, x, bphi, btheta
    ndof=4  ##ne, x, bphi, btheta
    denom= ndata - ndof  # Denominator used below in reduced chi^2
    
    #######################################################################
    ##
    ## match sobs data with sdb using the reduced chi^2 method

    ## We do not use Stokes V to define the geometry. Therefore we use a subset (1,2,4,5,6) 
    ## of Stokes parameters QU (line 1) and IQU (line 2)
    #print("{:4.6f}".format(time.time()-start),' SECONDS FOR PREAMBLE')
    
    #start=time.time()
#     diff=np.empty((sdb.shape[0],6),dtype=np.float64)
#     diff[:,0:2]= (sdb[:,1:3]-sobs[1:3]) / rms[1:3] 
#     diff[:,2:5]= (sdb[:,4:7]-sobs[4:7]) / rms[4:7] 
    diff=np.empty((database.shape[0],6),dtype=np.float64)
    diff[:,0:2]= (database[:,1:3]-sobs[1:3]) / rms[1:3] 
    diff[:,2:5]= (database[:,4:7]-sobs[4:7]) / rms[4:7]     
    
    #print("{:4.6f}".format(time.time()-start),' SECONDS FOR ALLOC DIFF')
    
    #These two above lines are the slowest of all the code due to requiring two operations. 
    #Could not find a better numba compatible alternative
    #Trials are in the old notebook
        
    #start=time.time()
    diff*=diff ## pure python multiplication is faster than numpy or any power operation
    #print("{:4.6f}".format(time.time()-start),' SECONDS FOR SQUARE DIFF')
    
    #start=time.time()
    chisq=np.empty((diff.shape[0]))
    np.round_( np.sum( diff, axis=1)/denom,15,chisq)   
    ## chisq needs to be initialized as the output variable, as np.round_ is not supported as a addresable function
    ## Precision is up to 15 decimals.
    ## Significant truncation errors appear if this is not enforced. if no rms and counts are present, the truncation is too low without the decimal fix.
    ## if no rms and counts are present, the truncation is too low without the decimal fix.
    #print("{:4.6f}".format(time.time()-start),' SECONDS FOR CHI^2')
    
    ## Need to create a numba compatible fast sorting as np.argpartition does not exist!
    ## subst_sort_par is a simple parallel function that produces an output similar to np.argpartition with sort.
    ## it just does a standard sorting search, but sorts just the first nsearch elements, making it extremely fast.
    #start=time.time()
    asrt=subst_sort_par(np.copy(chisq),nsearch) 
    #print(asrt)
    ## chisq needs to be passed as np.copy regardless if the sorting is done locally or in the subfunction.
    ## There is no apparent explanation for this behaviour
    #print("{:4.6f}".format(time.time()-start),' SECONDS FOR SORT')
    #######################################################################    
    ##
    ## the while loop will compute the first nsearch solutions or
    ## solutions that have chi^2 less than maxchisq, whichever comes first. 
    #start=time.time()    
    ## initialize return indices and chisq of nearest solutions
    si=0   ## needs to start -1 so updates to 0 i nthe loop because we have no functional do-while in python with numba.
    ix=asrt[0]
    ixr=asrt[0]
    ixrchisq=chisq[ixr]
    
    ##returns an array-like 0 vector of the dimensions of the requested output.
    ## the index is set to -1 to warn about hte missing entry
#     if chisq[ixr] > maxchisq :
#         nullout=np.zeros((nsearch,11),dtype=np.float64)
#         nullout[:,0]=np.full(nsearch,-1)
#         if verbose >=2: print("Warning: No solutions compatible with the maximum Chi^2 constraint.")
#         return nullout,np.zeros((nsearch,8),dtype=np.float64) 

    out=np.zeros((nsearch,11),dtype=np.float64)                
    out[:,0]=np.full(nsearch,-1)
    sfound=np.zeros((nsearch,8),dtype=np.float64)
    for si in range(nsearch):       
        ix=asrt[si]
        ixr=asrt[si]
        ixrchisq=chisq[ixr]
        #print(ixr,ixrchisq,maxchisq)
        if ixrchisq <=  maxchisq:
        
            if reduction  == 1:       # here if reduced we must get the original value of ix
                cgridnew=dbcgrid + 0.0
                cgridnew[3]=nsearch
                i,j,k,l=params_par(ixr,cgridnew)
                newl=np.int_(outindex[k,l])
                ix=invparams_par(i,j,k,newl,dbcgrid) # = original value ix
            #print(sobs[3]/(database[ixr,3]+1e-8),sobs[7]/((database[ixr,7]+1e-8)))
            bfield=sobs[3]/(database[ixr,3]+1e-8) # Magnetic field strength from ratio of sdb  with V/I of first line 
            #bfield=(sobs[3]/(database[ixr,3]+1.e8)+
            
            out[si,0]=ix
            out[si,1]=ixrchisq
            out[si,2:]=phys_par(ix,yobs,dbhdr,bfield)
            sfound[si,:]=database[ixr,:] 

# ####### fixing bug for very bad solutions with enormous chi2!                
#     for si in range(nsearch):       
#         ix=asrt[si]
#         ixr=asrt[si]
#         ixrchisq=chisq[ixr]
#         print(ixrchisq,maxchisq)
#         if ixrchisq <=  maxchisq:
        
#             if reduction  == 1:       # here if reduced we must get the original value of ix
#                 cgridnew=dbcgrid + 0.0
#                 cgridnew[3]=nsearch
#                 i,j,k,l=params_par(ixr,cgridnew)
#                 newl=np.int_(outindex[k,l])
#                 ix=invparams_par(i,j,k,newl,dbcgrid) # = original value ix
#             #print(database[ixr,3])
#             bfield=sobs[3]#/database[ixr,3] # Magnetic field strength from ratio of sdb  with V/I of first line 

#             if si == 0:              # Initialize arrays on first iteration                      

#                 out=np.zeros((nsearch,11),dtype=np.float64)
#                 out[0,0]=ix
#                 out[0,1]=ixrchisq
#                 out[0,2:]=phys_par(ix,yobs,dbhdr,bfield)
#                 sfound=np.zeros((nsearch,8),dtype=np.float64)   ##8--> fixed for a two line obs at this point
#                 sfound[0,:]=database[ixr,:]     # use reduced data index
                                
#             else:                     # append infor for each iterationn

#                 out[si,0]=ix
#                 out[si,1]=ixrchisq
#                 out[si,2:]=phys_par(ix,yobs,dbhdr,bfield)
#                 sfound[si,:]=database[ixr,:]


#     while ixrchisq <= maxchisq and si < nsearch:       
#         ix=asrt[si]
#         ixr=asrt[si]
#         ixrchisq=chisq[ixr]
        
#         if reduction  == 1:       # here if reduced we must get the original value of ix
#             cgridnew=dbcgrid + 0.0
#             cgridnew[3]=nsearch
#             i,j,k,l=params_par(ixr,cgridnew)
#             newl=np.int_(outindex[k,l])
#             ix=invparams_par(i,j,k,newl,dbcgrid) # = original value ix

#         bfield=sobs[3]/sdb[ixr,3] # Magnetic field strength from ratio of sdb  with V/I of first line 

#         if si == 0:               # Initialize arrays on first iteration             
#             out=np.zeros((nsearch,11),dtype=np.float64)
#             out[0,0]=ix
#             out[0,1]=ixrchisq
#             out[0,2:]=phys_par(ix,yobs,dbhdr,bfield)
#             sfound=np.zeros((nsearch,8),dtype=np.float64)   ##8--> fixed for a two line obs at this point
#             sfound[0,:]=sdb[ixr,:]     # use reduced data index
#         else:                      # append infor for each iterationn 
#             out[si,0]=ix
#             out[si,1]=ixrchisq
#             out[si,2:]=phys_par(ix,yobs,dbhdr,bfield)
#             sfound[si,:]=sdb[ixr,:]
#         si+=1

    #print("{:4.6f}".format(time.time()-start),' SECONDS FOR OUT')        
#     while chisq[ixr] <= maxchisq and si < nsearch:       
#         ixr=asrt[si]
#         dd=chisq[ixr]

#         if(dd <= maxchisq): 

#             if reduction  == 1:       # here if reduced we must get the original value of ix
#                 cgridnew=dbcgrid + 0.0
#                 cgridnew[3]=nsearch
#                 i,j,k,l=params_par(ixr,cgridnew)
#                 newl=np.int_(outindex[k,l])
#                 ix=invparams_par(i,j,k,newl,dbcgrid) # = original value ix
            
#             bfield=sobs[3]/sdb[ixr,3] # Magnetic field strength from ratio of sdb  with V/I of first line 

            
#             if si == 0:               # Initialize arrays on first iteration             
#                 out=np.zeros((nsearch,11),dtype=np.float64)
#                 out[0,0]=ix
#                 out[0,1]=dd
#                 out[0,2:]=phys_par(ix,yobs,dbhdr,bfield)
#                 sfound=np.zeros((nsearch,8),dtype=np.float64)   ##8--> fixed for a two line obs at this point
#                 sfound[0,:]=sdb[ixr,:]     # use reduced data index
#             else:                      # append infor for each iterationn 
#                 out[si,0]=ix
#                 out[si,1]=dd
#                 out[si,2:]=phys_par(ix,yobs,dbhdr,bfield)
#                 sfound[si,:]=sdb[ixr,:]
       
#         si+=1
    #print(out.shape,sfound.shape)
    return out,sfound





###########################################################################    
###########################################################################     


@njit
def get_subset_par(sobs,dbhdr,database,nsearch,verbose): 
## returns a subset of the sdb array compatible with  height in yobs
##
    #start=time.time()  
    #######################################################################
    ## unpack dbcgrid parameters from the database accompanying db.hdr file 
    ##
    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    bthetamin, bthetamax, nline, wavel  = dbhdr 
    
    ##    
    #print(database.flags)
    #database=np.ascontiguousarray(database)                  # fix database memory allocation 
    #print(database.flags)
    outindex=np.zeros((nbphi,nsearch))
    
    #######################################################################
    ## reduced CLE array calculation
    ##
    kk=np.arange(0,nbphi)
    bphir  = bphimin + kk*(bphimax-bphimin)/(nbphi-1)        # cle phi array
    
    ll=np.arange(0,nbtheta)
    bthetar=bthetamin + ll*(bthetamax-bthetamin)/(nbtheta-1) #cle theta array
    
    dthb=(bthetamax-bthetamin)/nbtheta
    tt = np.tan(bthetar)                                     #cle tan theta array

    ##   Below, the relation is used: tan Phi_B = sin phi * tan theta
    phib_obs=-0.5*np.arctan2( sobs[2],sobs[1] )
    tphib_obs= np.tan(phib_obs)  ##   here is tan Phi_B:
    
    ## use reshape on the fly to get datasel.
    ## we don't need more than the desired subsets
    #shape=(ned,ngx,nbphi,nbtheta,8)
    #datasel=(np.reshape(database,shape) )[:,:,:,:nsearch,:]   # selection of nsearch subset
    #datatemp=np.reshape(database,shape)
    datasel=np.zeros((ned,ngx,nbphi,nsearch,8),dtype=np.float64)
    #print("{:4.6f}".format(time.time()-start),' SUBSET: SECONDS FOR PREAMB') 
    #print(datatemp.flags)

    #start=time.time()     
    ## Find those indices compatible with phib observed
    r2d=180./np.pi 
    for amb in range(2):
        #start1=time.time()
        ambphib=phib_obs+amb*np.pi/2.
        for ir in range(bphir.shape[0]):          # loop over bphi in cle frame
            ttp= tt * np.sin(bphir[ir])
            diff=r2d*np.abs(ambphib-np.arctan(ttp)) # this is an array over btheta
            srt=subst_sort_par(np.copy(diff),nsearch)            
            
            #start0=time.time()
            if ir + srt[0] > 0:                  # important to avoid phi=0 and theta = 0 case
                for jj in range(srt.shape[0]):     # advanced slicing is not available, the for+enumerate comes from numba requirements
                    #tt1=(np.reshape(database,shape))[:,:,ir,srt[jj],:].copy() 
                    ## no explanation here for copy(). 
                    ## Something weird happens memory-wise. tt1 might be a view instead of a separate object in that case?
                    #datasel[:,:,ir,jj,:]=tt1
                    #datasel[:,:,ir,jj,:]=np.reshape(database,shape)[:,:,ir,srt[jj],:]
                    datasel[:,:,ir,jj,:]=database[:,:,ir,srt[jj],:]
                    outindex[ir,jj]=srt[jj]
            #print("{:4.6f}".format(time.time()-start0),' SUBSET LOOP: if') 
        #print("{:4.6f}".format(time.time()-start1),' SUBSET LOOP: for')     

    ##updated version of the linear polarization match block above from tries in cledb_proc.
    ##_________________________________________________________________________________________

    # # Find those indices compatible with phib observed
    # for amb in range(2):                                   ## take all ambiguous solutions from tangent
    #     ambphib=phib_obs+amb*np.pi/2.
    #     for ir in range(bphir.shape[0]):                   ## loop over bphi in cle frame
    #         ttp = tt * np.sin(bphir[ir])
    #         difft = r2d*np.abs(ambphib - np.arctan(ttp))    ## this is an array over btheta
    #         srt=cledb_partsort(difft,nsearch)               ## no speed gain to move to np.sort here as the arrays are small.
    #         ## advanced slicing is not available, the for jj enumeration comes from numba requirements
    #         if ir + srt[0] > 0:                             ## important to avoid phi=0 AND theta = 0 case
    #             for jj in range(nsearch):              ## note nsearch = srt.shape[0]
    #                 datasel[:,:,ir,jj,:]=database_sel[:,:,ir,srt[jj],:]
    #                 outredindex[ir,jj]=srt[jj]    ## Find those indices compatible with phib observed

    # # Find those indices compatible with phib observed
    # for ir in range(bphir.shape[0]):                   ## loop over bphi in cle frame
    #     ttp = tt * np.sin(bphir[ir])
    #     difft = np.abs(tphib_obs - ttp)    ## this is an array over btheta
    #     srt=cledb_partsort(difft,2*nsearch)               ## no speed gain to move to np.sort here as the arrays are small.
    #     ## advanced slicing is not available, the for jj enumeration comes from numba requirements
    #     if ir + srt[0] > 0:                             ## important to avoid phi=0 AND theta = 0 case
    #         for jj in range(nsearch):              ## note nsearch = srt.shape[0]
    #             datasel[:,:,ir,jj,:]=database_sel[:,:,ir,srt[jj],:]
    #             outredindex[ir,jj]=srt[jj]    ## Find those indices compatible with phib observed
    #     #print("a:",srt)

    # # for ir in range(bphir.shape[0]):                   ## loop over bphi in cle frame
    #     ttp = tt * np.sin(bphir[ir])
    #     difft = np.abs(tphib_deg_obs - ttp)    ## this is an array over btheta
    #     srt=cledb_partsort(difft,2*nsearch)               ## no speed gain to move to np.sort here as the arrays are small.
    #     ## advanced slicing is not available, the for jj enumeration comes from numba requirements
    #     if ir + srt[0] > 0:                             ## important to avoid phi=0 AND theta = 0 case
    #         for jj in range(nsearch):              ## note nsearch = srt.shape[0]
    #             datasel[:,:,ir,jj,:]=database_sel[:,:,ir,srt[jj],:]
    #             outredindex[ir,jj]=srt[jj]    ## Find those indices compatible with phib observed
    #     #print("b:",srt)


    # # Find those indices compatible with phib observed
    # for ir in range(bphir.shape[0]):                   ## loop over bphi in cle frame
    #     ttp = tt * np.sin(bphir[ir])
    #     difft = np.abs(tphib_obs - ttp)    ## this is an array over btheta
    #     srta=cledb_partsort(difft,2*nsearch)               ## no speed gain to move to np.sort here as the arrays are small.
    #     difftt = np.abs(tphib_deg_obs - ttp)    ## this is an array over btheta
    #     srtb=cledb_partsort(difftt,2*nsearch)               ## no speed gain to move to np.sort here as the arrays are small.
    #     ## advanced slicing is not available, the for jj enumeration comes from numba requirements
    #     if ir + srta[0] > 0 or  ir + srtb[0] > 0:                             ## important to avoid phi=0 AND theta = 0 case
    #         for jj in range(2*nsearch):              ## note nsearch = srt.shape[0]
    #             datasel[:,:,ir,jj,:]=database_sel[:,:,ir,srta[jj],:]
    #             outredindex[ir,jj]=srta[jj]    ## Find those indices compatible with phib observed
    #             datasel[:,:,ir,jj+(2*nsearch),:]=database_sel[:,:,ir,srtb[jj],:]
    #             outredindex[ir,jj+(2*nsearch)]=srtb[jj]    ## Find those indices compatible with phib observed

    ##WORKING
#     # Find those indices compatible with phib observed
#     for ir in range(bphir.shape[0]):                   ## loop over bphi in cle frame
#         ttp = tt * np.sin(bphir[ir])
#         diffa = np.abs(tphib_obs - ttp)    ## this is an array over btheta
#         diffb = np.abs(tphib_deg_obs - ttp)    ## this is an array over btheta
#         diffab = np.append(diffa,diffb)
#         srt=cledb_partsort(diffab,4*nsearch)
#         ## advanced slicing is not available, the for jj enumeration comes from numba requirements
#         if ir + srt[0] > 0:                             ## important to avoid phi=0 AND theta = 0 case
#             for jj in range(4*nsearch):              ## note nsearch = srt.shape[0]
#                 datasel[:,:,ir,jj,:]=database_sel[:,:,ir,srt[jj],:]
#                 outredindex[ir,jj]=srt[jj]    ## Find those indices compatible with phib observed




    # # Find those indices compatible with phib observed
    # difft=np.zeros((1),dtype=np.float64)
    # for ir in range(bphir.shape[0]):                   ## loop over bphi in cle frame
    #     ttp = tt * np.sin(bphir[ir])
    #     for it in range(tphib_obs.shape[0]):
    #         difft = np.append(difft,(np.abs(tphib_obs[it] - ttp))**2)                    ## this is an array over btheta
    #     srt = cledb_partsort(difft,nsearch)               ## no speed gain to move to np.sort here as the arrays are small.
    #     #print("partsort:", srt[:nsearch])
    #     ## advanced slicing is not available, the for jj enumeration comes from numba requirements
    #     if ir + srt[0] > 0:                             ## important to avoid phi=0 AND theta = 0 case
    #       for jj in range(nsearch):              ## note nsearch = srt.shape[0]
    #             datasel[:,:,ir,jj,:]=database_sel[:,:,ir,srt[jj],:]
    #             outredindex[ir,jj]=srt[jj]    ## Find those indices compatible with phib observed




    # #ttp=np.zeros((bphir.shape[0]),dtype=np.float64)
    # difft=np.array([1e7],dtype=np.float64)
    # for ir in range(bphir.shape[0]):                   ## loop over bphi in cle frame
    #     ttp= tt * np.sin(bphir[ir])
    #     difft =np.append(difft,r2d*np.abs(phib_obs - np.arctan(ttp)))
    #     difft =np.append(difft,r2d*np.abs(phib_obs+np.pi/2. - np.arctan(ttp)))    ## this is an array over btheta
    # print(bphir.shape[0],difft.shape)
    # srt=cledb_partsort(difft,nsearch)               ## no speed gain to move to np.sort here as the arrays are small.
    # ## advanced slicing is not available, the for jj enumeration comes from numba requirements
    # #if ir + srt[0] > 0:                             ## important to avoid phi=0 AND theta = 0 case
    # for jj in range(nsearch):              ## note nsearch = srt.shape[0]
    #     datasel[:,:,ir,jj,:]=database_sel[:,:,ir,srt[jj],:]
    #     outredindex[ir,jj]=srt[jj]
    # print(difft[srt][0:nsearch])    
    
    ##_____________________________________________________________________________________
    
    #print("{:4.6f}".format(time.time()-start),' SUBSET: SECONDS FOR LOOP')     
    # now work with reduced dataset with nbtheta replaced by nsearch
    # data3 / datasel is not a memory gobbler.
    #start=time.time()     
    #shapesig=(ned*ngx*nbphi*nsearch,8)
    #data3=np.reshape(np.ascontiguousarray(datasel),shapesig)
    #sdb=np.reshape(np.ascontiguousarray(datasel),shapesig)
    
    # define new smaller array
    #sdb=np.append(data3[:,0,:4],data3[:,1,:4],axis=1)
    
    #print("{:4.6f}".format(time.time()-start),' SUBSET: SECONDS FOR ARANGE')      
    if verbose >=1:print('Search over theta reduced by a factor: ', np.int_(nbtheta/nsearch),". New db size: (",ned*ngx*nbphi*nsearch,",8)")
    #dt= "{:4.3f}".format(time.time()-start)
    #if verbose ==1:print(dt,' SECONDS FOR REDUCE (loop phi) CALC')    
    return np.reshape(datasel,(ned*ngx*nbphi*nsearch,8)),outindex
    #return sdb,outindex
###########################################################################       
###########################################################################


@jit(parallel=True,forceobj=True,looplift=True)    
def db_find(y,dbdir,sx,sy,verbose):
## Main script to find and preload all necesary databases. 
## Returns databases as a list along with an encoding coresponding to each voxel of the observation
## this is made compatible to numba object mode with looplifting to read all the necesary database in a parallel fashion
    if verbose >=1: 
        print('------------DB READ START-----------')
        if verbose >= 3:start0=time.time()        
    ######################################################################
    ## Preprocess the string information
    names,numbers,subdirs=file_preprocess(dbdir)
    if verbose >=1:print("DB covers a span of ",numbers.shape[0]," heights between",1+np.min(numbers)/1000,"-",1+np.max(numbers)/1000,"Solar radius")
    
    ######################################################################
    ## Create a file encoding showing which database to read for which voxel
    ## the encoding at this stage does not have na order.
    fenc=np.empty((sx,sy),dtype=np.int_)
    
    for xx in range(sx):
        for yy in prange(sy):                     
            fenc[xx,yy]=find_elongation_num_par(y[xx,yy],numbers)   
    
    sz=np.unique(fenc) ## makes a list of unique databases to read
    if verbose >= 1: print("Load ",sz.shape[0], " DB heights in memory\n----------")
    
    ######################################################################
    ## read the required databases and their header
    
    
    ## numpy large array implementation does not parallelize properly leading to a 5x increase in runtime per 1024 calculations
    ## reverted to use a list to feed the database set to the calculation
    
    if subdirs=="twolineDB":
    
        ## preprocess and return the database header information
        ## when multiple databases we assume all are of the same size.
        ##read the header and dimensions and just pass them as parameters to function calls
        g=readhdr_par(dbdir+'db.hdr')
        dbhdr=[parsehdr_par(g)][0] ## dont understand why the 0 indexing is needed...    
    
        database0=[None]*sz.shape[0]
        for ii in prange(sz.shape[0]):
            database0[ii]=sdbread_par(names[0][sz[ii]],dbhdr,verbose)
            for ij in range(3,-1,-1):
                database0[ii][:,:,:,:,ij]=database0[ii][:,:,:,:,ij]/database0[ii][:,:,:,:,0]
    
    else:
        g=readhdr_par(dbdir+subdirs[0]+'db.hdr') ## assuming the same header info; reading the first DB header
        dbhdr=[parsehdr_par(g)][0]
        
        database0=[None]*sz.shape[0]
        for ii in prange(sz.shape[0]):
            database0[ii]=np.append(sdbread_par(names[0][sz[ii]],dbhdr,verbose),sdbread_par(names[1][sz[ii]],dbhdr,verbose),axis=4)   
            for ij in range(7,-1,-1):
                database0[ii][:,:,:,:,ij]=database0[ii][:,:,:,:,ij]/database0[ii][:,:,:,:,0]

    ## standard reflacted lists will be deprecated in numba 0.54; currently running 0.53. this is a fix!
    database = List()    
    [database.append(x) for x in database0]   

    ## Put the encodings in a readable order# no longer needed
    cnt=0
    for i in sz:
        fenc[np.where(fenc == i)] = cnt
        cnt+=1    
    if verbose >=1:
        if verbose >= 3:
            print("{:4.6f}".format(time.time()-start0),' SECONDS FOR TOTAL DB SEARCH AND FIND')
        print('--------------DB READ FINALIZED----------------')
    return fenc,database,dbhdr
###########################################################################
###########################################################################



###########################################################################
### UTILS for file manipulation and database loading ######################
###########################################################################

def file_preprocess(dbdir):  
#returns filename and index of elongation y in database  
## This prepares the database directory files outside of the numba jitted functions as string operations do not play nice with numba.
    line_str=["fe-xiii_1074/","fe-xiii_1079/","mg-viii_3028/","si-ix_3934/","si-x_1430/"]
    line_bin=[]
    dbfolders=[]
    for i in line_str:
        line_bin.append(os.path.isdir(dbdir+i))
        
    if sum(line_bin) >=2:
        for i in range(len(np.where(line_bin)[0])):
            dbfolders.append(line_str[np.where(line_bin)[0][i]])

        namesA=glob.glob(dbdir+dbfolders[0]+"DB*.DAT")
        namesB=glob.glob(dbdir+dbfolders[1]+"DB*.DAT")
        nn=str.find(namesA[0],'DB')
        numbers=np.empty(len(namesA),dtype=np.int_)
        for i in range(0,len(namesA)):
            numbers[i]=np.int_(namesA[i][nn+2:nn+6])
        return [namesA,namesB],numbers,dbfolders
        
    elif os.path.isfile(dbdir+"db.hdr"):
        names=glob.glob(dbdir+"DB*.DAT")
        nn=str.find(names[0],'DB')
        numbers=np.empty(len(names),dtype=np.int_)
        for i in range(0,len(names)):
            numbers[i]=np.int_(names[i][nn+2:nn+6])
        return [names,names],numbers,'twolineDB' ##double return of names is superfluous; reason is to keep returns consistent regardless of the if case
    
    elif sum(line_bin) ==1:
        print("Requires two individual databases in directory")    
    else: print("No database to read in directory")
###########################################################################
###########################################################################
        
def sdbread_par(fildb,dbhdr,verbose):
## here reading data is done and stored in variable database of size (ncalc,line,4)
## due to calling dcompress nad np.fromfile, this is incompatible with numba    
    if verbose >=3: start=time.time()    
    
    ######################################################################
    ## unpack dbcgrid parameters from the database accompanying db.hdr file 
    ##  
    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    bthetamin, bthetamax, nline, wavel  = dbhdr 

    ## here reading data is done and stored in variable database of size (n,2*nline)

    if verbose >=1:
        print('READING INDIVIDUAL DATA...')
        print("DB file location:", fildb)      
   
    db=dcompress(np.fromfile(fildb, dtype=np.int16),verbose)

    if verbose >= 3: print("{:4.6f}".format(time.time()-start),' SECONDS FOR INDIVIDUAL DB READ\n----------')

    #return np.reshape(db,(ned*ngx*nbphi*nbtheta,nline*4))
    return np.reshape(db,(ned,ngx,nbphi,nbtheta,nline*4))
###########################################################################       
###########################################################################

def dcompress(i,verbose):
## CLEDB writes compressed databases to save storage space.
## This is numba incompatible due to the sys. and ne.evaluate functions

    cnst=-2.302585092994046*15./32767.
    ## Constants here must correspond to those in dbe.f in the CLE main directory
    ## 32767 is the limit of 4-byte ints
    ## 15 is the number of orders of magnitude for the range of intensities
    
    if verbose >=1: print(np.int_(sys.getsizeof(i)/1.e6)," MB in DB file")
    if verbose >=2: 
        if np.int_(sys.getsizeof(i)/1.e6) > 250 : print("Warning: Very large DB; Processing will be slow!")
    negv=np.flatnonzero(i  < 0 )


    f=np.abs(i)*cnst
    f=ne.evaluate("exp(f)")

    if size(negv) > 0:f[negv]=-f[negv]
    del negv

    return f
###########################################################################       
###########################################################################
      

@jit(parallel=True,forceobj=True)    #
def find_elongation_num_par(y,numbers):
# returns index of elongation y in database
# The script likes it more fed as a function rather than an inline calculation! 
           
    return np.argmin(np.abs(1.+ numbers/1000-y))
###########################################################################
###########################################################################


@jit(forceobj=True)
def readhdr_par(file):
## just read the header of the DB*DAT data cubes
## forced as a numba jit object due to np.fromfile
    g=np.fromfile(file,dtype=np.float32,sep=' ')
           
    return g
###########################################################################
###########################################################################


@njit
def parsehdr_par(g):
## parses the header and returns the parameters contained in a database directory           
## db.hdr text file needs to have the specific version format decribed in the CLEDB 2.0.2 readme

## two kinds of data are returned
## 1. the linear coefficents of the form min, max, nx,theta,phi
## 2. logarithmically spaced parameters xed  (electron density array)
## for these also the min and max and number are also returned.
    ned=np.int_(g[0]) 
    ngx=np.int_(g[1])
    nbphi=np.int_(g[2])
    nbtheta=np.int_(g[3])
    emin=np.float64(g[4])
    emax=np.float64(g[5])
    xed=lcgrid_par(emin,emax,ned)    # Ne is log
    gxmin=np.float64(g[6])
    gxmax=np.float64(g[7])
    bphimin=np.float64(g[8])
    bphimax=np.float64(g[9])
    bthetamin=np.float64(g[10])
    bthetamax=np.float64(g[11])
    nline=np.int_(g[12])
    wavel=np.empty(nline,dtype=np.float64)
    
    for k in range(0,nline): wavel[k]=g[13+k]    

    return g, ned, ngx, nbphi, nbtheta, xed, gxmin,gxmax, bphimin, bphimax,\
        bthetamin, bthetamax, nline, wavel 
###########################################################################       
###########################################################################
             

###########################################################################
### Secondary fast processing UTILS #######################################
###########################################################################

@njit
def subst_sort_par(arr,nsearch):
## Simplest possible substitution partial sort. It just returns the first nsearch sorted indexes similar to np.argpartition. 
## It is numba non-python compatible!     
    
    asrt=np.zeros((nsearch),dtype=np.int_)
    for i in range(nsearch):
        a=np.int_(np.argmin(arr[i:])+i)
        asrt[i]=a 
        arr_temp=arr[i]
        arr[i]=arr[a]
        arr[a]=arr_temp
           
    return asrt         ##works like a charm!
###########################################################################       
###########################################################################


@njit
def lcgrid_par(mn,mx,n):
## grid for ned densities from mn to mx logarithms of density.
    
    ff=10.** ((mx-mn)/(n-1))
    xf=np.empty(n,np.float64)
    xf[0]=10.**mn
    for k in range(1,n): xf[k]=xf[k-1]*ff
           
    return xf
###########################################################################      
###########################################################################

           
@njit 
def params_par(index,dbcgrid):
## for index, get i,j,k,l  indices
    
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

    return i,j,k,l
###########################################################################      
###########################################################################   

           
@njit 
def invparams_par(i,j,k,l,dbcgrid):
# for i,j,k,l get index
# Reverse function of params_par
           
    ned=np.int_(dbcgrid[0]) 
    ngx=np.int_(dbcgrid[1])
    nbphi=np.int_(dbcgrid[2])
    nbtheta=np.int_(dbcgrid[3])
           
    return np.int_(i*ngx*nbphi*nbtheta + j*nbphi*nbtheta + k*nbtheta + l)
###########################################################################      
########################################################################### 


@njit
def electron_density_par(r):
## computes the electron radial density using the Baumbach formulation
    baumbach = 1.e8*(0.036/r**1.5 + 1.55/r**6.)
    hscale=   7.18401074e-02  # 50 Mm
           
    return np.float64(3.e8*np.exp(- (r-1.)/hscale) + baumbach)
###########################################################################      
###########################################################################


@njit 
def phys_par(index,gy,dbhdr,b):
## Returns the lvs and LOS geometry and magnetic field physics
    phs=physa_par(index,gy,dbhdr,b)
    bphi=phs[3]
    btheta=phs[4]
    bx=b*np.sin(btheta)*np.cos(bphi)
    by=b*np.sin(btheta)*np.sin(bphi)
    bz=b*np.cos(btheta)
           
    return np.array((phs[0],phs[1],phs[2],b,bphi,btheta,bx,by,bz),dtype=np.float64)
###########################################################################      
###########################################################################


@njit 
def physa_par(index,gy,dbhdr,b):
## Helper for phys_par           
## returns the LVS physics and geometry.

    dbcgrid, ned, ngx, nbphi, nbtheta,  xed, gxmin,gxmax, bphimin, bphimax, \
    bthetamin, bthetamax, nline, wavel  = dbhdr 

    i,j,k,l = params_par(index,dbcgrid)

    gx = gxmin + j*(gxmax-gxmin)/(ngx-1)
    bphi  = bphimin + k*(bphimax-bphimin)/(nbphi-1)
    btheta  = bthetamin + l*(bthetamax-bthetamin)/(nbtheta-1)
    ed = np.float64(xed[i]* electron_density_par(np.sqrt(gy*gy+gx*gx)))    # log of Ne

    return np.array((np.log10(ed),gy,gx,bphi,btheta,b),dtype=np.float64)
###########################################################################      
###########################################################################
