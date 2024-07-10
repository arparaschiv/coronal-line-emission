#######################################################################
# Document name: atmrd.py
# Created by:    Alin Paraschiv arparaschiv@nso.edu
# Original IDL version: Philip Judge
#
# Change log
# 20190827 ARP ported the IDL version to python
#######################################################################
#
#
# Read in atmospheric model for python
# CALLING EXAMPLE:python3 atmrd.py 'path/to/ATMOS_file'
#


import struct
import numpy as np
import sys

## check if input arggument string containing path exists
if len(sys.argv) > 1:
    file = sys.argv[1]
else:
    file = './ATMOS' # if no input argument assumes atmos is in current dir

## read the atmoshpere header to get the data dimensions.
with open(file, 'br') as f3:
    buffer = f3.read()
    print ("Length of buffer is %d" % len(buffer))
    data_h = struct.unpack_from("<6s1d2f1l1d2f1l1d2f1l",buffer, offset=4)
## The format of the file contains 1 empty space, 6 string characters followed by a series of three sets (for x, y, and z directions) 
## of 1 double (empty), two floats (g_min, g_max), and long (ng_) 
## Format example: "<1f6s1d2f1l1d2f1l1d2f1l1d3l1d7f + repeat 3l1d7f as many times as dimensions allow"

cname=data_h[0]
gxmin,gxmax,ngx = float('%.8f'%data_h[2]),float('%.8f'%data_h[3]),data_h[4]  ##truncated so there are no insignificant trailing digits
gymin,gymax,ngy = float('%.7f'%data_h[6]),float('%.8f'%data_h[7]),data_h[8]
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

## put the data in a nice class so it mimics IDL structures
class atmos_params:
	file = file       # file:     file containing atmospheric parameters
	model = cname       # model:    name of fortran routine used to determine coronal parameters
	x     = None       # x,y,z:    positions in solar radii',$
	y     = None        #           (x is along LOS, y is E-W, z N-S in observer''s frame)
	z     = None
	bx    = None        # bx,by,bz: x,y,z-components of magnetic field in G
	by    = None
	bz    = None
	nne   = None        # nne:      electron density in /cm3
	te    = None        # te:       electron temperature in K
	v     = None        # v:        LOS velocity in km/s (+ve is a red-shift)
	vt    = None        # vt:       "turbulent" velocity in km/s

atm=atmos_params()
atm.x=x
atm.y=y
atm.z=z
atm.bx=bx
atm.by=by
atm.bz=bz
atm.nne=nne
atm.te=te
atm.v=v
atm.vt=vt

###other shit
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.ion()

bx1=np.sum(atm.bx,axis=0)
by1=np.sum(atm.by,axis=0)
bz1=np.sum(atm.bz,axis=0)
nne1=np.sum(atm.nne,axis=0)

btot=np.sqrt(bx1*bx1+by1*by1+bz1*bz1)

#######################################################################
## End of 'atmrd.py'
#######################################################################
