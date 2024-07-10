#######################################################################
## Document name: outrd.py
## Created by:    Alin Paraschiv arparaschiv@nso.edu
## Original IDL version: Philip Judge
##
## Change log
## 20190903 ARP ported the IDL version to python
#######################################################################

## Read in the CLE output for python 
## For output from cle (OUT file) queries contact: judge@ucar.edu

## CALLING EXAMPLE: python3 -i outrd.py s kr 'path/to/OUT_file' no_rot

## VARIABLE DESCRIPTION:
##     s is a map structure, where:
##     s.data1 is the intensity I integrated over the a line
##     s.data2 is Q or Px integrated over the line
##     s.data3 is U or Py integrated over the line
##     s.data4 is V integrated over the line (using weights of -1,+1 for wavelengths < and > line center at rest)
##     s.data5 is Vmag=V but computed using the magnetograph formula
##     s.data-q are delta-wavelengths across line profile [Angstrom]
##     s.data-full are stokes vectors full(iquv,q,x,y) units are erg/cm^2/s/Angstrom])
##
##     kr is the transition index (the atom contains several radiative transitions, kr is the
##
##     file takes the path to the CLE OUT file.
##
##     no_rot is a binary argument used to not rotate q and u vectors from CLE.
##            The CLE routine corona.f contains the definition of the reference direction.
##            The Stokes computed using a ref direction parallel to the z-axis [N-S] to a 
##            polarization vector in the frame of (y,z) [E-W, N-S] in the plane of the sky. 
##.           The input numbers q,u are replaced with these components.
##            


import struct
import numpy as np
import sys
import os
import math
import re

#for ulterior plotting
import matplotlib.pyplot as plt
import matplotlib as mpl
plt.ion()

# insert at 1, 0 is the script path (or '' in REPL)
#sys.path.insert(1, '/home/alin/Documents/physics_prog/cle/idl')
#import getwrd      ## local routine 

## check if input arggument string containing path exists
if 4<= len(sys.argv) <= 5:
	s = sys.argv[1]
	kr = int(sys.argv[2])
	if isinstance(sys.argv[3], str):
		file = sys.argv[3]
	else:
		if os.path.isfile('./OUT'):
			file = './OUT' # if no input argument assumes atmos is in current dir
		else:
			ValueError('CLE output file not found.')
	if len(sys.argv) == 5:
		no_rot=int(sys.argv[4])
else:
	raise ValueError('input arguments not correct type or not in right order.')

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
#	cornotes = ''
#	with open('OUT') as f:
#	for line in f:
#		if '1 CORONAL ROUTINE' in re.findall(r'\w+', line):


with open(file) as f:                   # Reads to "1 STOKES DATA OUTPUT"
	for line in f:
		if 'OUTPUTTED' in re.findall(r'\w+', line): 
			nline=int(re.findall(r'\w+\.?\w*', line)[re.findall(r'\w+\.?\w*', line).index('OUTPUTTED')+1] )
			break

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
							if len(sys.argv) < 4 or no_rot == 0: 
								if emerg[k,1,i,j] != 0 or emerg[k,2,i,j] != 0:
									alpha = 0.5*math.atan2(emerg[k,1,i,j],emerg[k,2,i,j])+math.pi/2
									p = math.sqrt(emerg[k,1,i,j]**2+emerg[k,2,i,j]**2)
									emerg[k,1,i,j]=p*math.cos(alpha)
									emerg[k,2,i,j]=p*math.sin(alpha)

## prepare map putputs
print(type(kr))
si  =   emerg[kr,0,:,:].reshape(emerg.shape[2],emerg.shape[3])
px  =   emerg[kr,1,:,:].reshape(emerg.shape[2],emerg.shape[3])
py  =   emerg[kr,2,:,:].reshape(emerg.shape[2],emerg.shape[3])
v   =   emerg[kr,3,:,:].reshape(emerg.shape[2],emerg.shape[3])
vmag=   emerg[kr,4,:,:].reshape(emerg.shape[2],emerg.shape[3])
p = np.sqrt(px**2+py**2)

#b_los from plowman eq. 14 units are cgs
cc=2.99792458e10 #cm s-1 
me= 9.1095e-28   #g
ee=4.8032e-10    #esu  g1/2 cm3/2 s−1
kb=1.38064852e-16  #erg/K
te= 10**6.22
mion=55.845*1.672621E-24  #  g    for FeXIII 
omegat=math.sqrt(2*kb*te/mion)     #thermal velocity without microturbulence
#dlmd_d=omegat*alamb[kr]*1e-8/cc               #Doppler width in frequency units!   1e-8 -->conversion A to cm
dlmd_d=omegat/(1e-8* alamb[kr])               #Doppler width in wavelength units!  1e-8 -->conversion A to cm
#blos= -8*math.pi*me*cc*v*dlmd_d/(3*(si+p)*ee)   #eq 14 LOS magnetic field |B|cos(theta) ##not finished here; will do separately

## output to maps

## units 
##can be either arcsec or solar radii; default is solar radii.
punit = 955.                     #arcsec
unit = 'arcseconds'
punit = 1.                       #solar radii
unit = 'solar radii'

xc = (gymax+gymin)/2*punit 
dx =  (gymax-gymin)/ngy*punit 
yc = (gzmax+gzmin)/2*punit 
dy =  (gzmax-gzmin)/ngz*punit 
x = (gymin + np.arange(ngy)*dx)*punit 
y = (gzmin + np.arange(ngz)*dy)*punit 

if no_rot !=0:
	sx='Q (ref=z-axis)'
	sy='U (ref=z-axis)'
else:
	sx='Px'
	sy='Py'

## Write the output map using a class
class out_params:
	file = file       # file:     file containing atmospheric parameters
	model = crtn       # model:    name of fortran routine used to determine coronal parameters
	x     = None       # x,y    positions in solar radii',
	y     = None       
	xc    = None
	yc    = None
	dx    = None
	dy    = None
	punit    = None
	unit     = None
	wave     = None
	dunits   = None
	data_q   = None
	data_full= None
	data1    = None
	data2    = None
	data3    = None
	data4    = None
	data5    = None
	data1_id  = None
	data2_id  = None
	data3_id  = None
	data4_id  = None
	data5_id  = None


s = out_params()
s.x=x
s.xc=xc
s.dx=dx
s.y=y
s.yc=yc
s.dy=dy
s.punit=punit
s.unit=unit
s.wave=alamb[kr]
s.dunits='erg/cm^2/s'
if iwline > 0: 
	s.data_q =   q[kr,0:nq[kr]]
	s.data_full= full[kr,:,:,:,:].reshape(full.shape[1],full.shape[2],full.shape[3],full.shape[4])
	print('Full stokes profiles were also outputted, stored e.g., in s.q, s.data_full')
s.data1   = si
s.data1_id = 'Stokes I'
s.data2   = px
s.data2_id = 'Stokes '+sx
s.data3   = py
s.data3_id = 'Stokes '+sy
s.data4   = v
s.data4_id = 'Stokes V'
s.data5   = vmag
s.data5_id = 'Stokes V (magnetograph formula)'



#######################################################################
## End of 'outrd.py'
#######################################################################