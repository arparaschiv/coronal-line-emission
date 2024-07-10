July 13, 2006

Notes on fortran-77 program "cle" (coronal line emission) by 
P. Judge.

Directories:
------------

current/       source code (current version)
data/          data needed by the program (atomic models, ionization
               balance, grid, input options)    
idl/           some procedures to read and plot output


Defining the magnetic and thermal coronal structure:
----------------------------------------------------

You must build a subroutine called "cuser.f" which has the 
same argument list and inputs/ outputs as the examples provided,
ffl.f, zl.f or dipole.f.  Consult the headers of these files.  
For solar coordinates used, see AAAREADME in the current/ directory. 

A routine "cuser.f" is included, it simply contains a copy of the 
dipole.f code. An input file DIPOLE.DAT is needed to run this case-
one can be found in the data/ directory. 



Compilation:
------------

use "make cle" in the current directory. executable cle is in the same
directory.

Running:
--------

Several input files are needed, examples are in the data/ directory

The following need not be changed:

IONEQ   equilibrium ionization fractions from CHIANTI

the following should be edited according to preference:

GRID.DAT  specify the uniform, cartesian grid of points on which to
          compute the emission line profiles. See the file
	  AAAREADME in the current/ directory for a description
	  of coordinates used in cle.

INPUT     specify various options within the program. See below

There may be other files which are read in order to specify the 
coronal parameters, for example, the file DIPOLE.DAT specifies the
nature of the dipolar field for this particular case. 

ATOM      atomic data.  Use one of the files atom.fe10, atom.fe13
          atom.fe1, atom.si9, atom.si10 in the data/ directory.
          (The files ending in "o" there are MUCH larger and MUCH
          more time consuming to compute).  


File INPUT:
-----------

The input file looks like the following:

TP2TE=1.0,WLMIN=5000.,WLMAX=40000.,SMALLN=1.E-5,QNORM=10.00,
CECOEF=0.,ISUM=1,ICOLL=1,ISPLIN=3,
IWATOM=1,IWLINE=0,IWEQI=1,IDEBUG=0,IWATMO=0,
CRTN=CUSER

The parameters have the following meaning:

TP2TE    ratio of proton to electron temperature
WLMIN    wavelength [angstroms] must be > this for output
WLMAX    wavelength [angstroms] must be < this for output
SMALLN   skip points with electron densities are < this
QNORM    normalization constant for Doppler widths/shifts
CECOEF*  a "fudge factor" to allow for missing elastic collisions
ISUM     which atomic level to use for particle conservation
ICOLL    include collisions (1), don't (0)
ISPLIN   choice of spline interpolation
IWATOM   write atomic model if 1
IWLINE   write line data (needed for output)
IWEQI    write ionization balance used
IDEBUG   verbose output if 1
IWATMO   write atmospheric data (magnetic & thermal data)
CRTN*    choice of coronal routine supplying data 

A USER WILL TYPICALLY WANT TO CHANGE ONLY THOSE MARKED WITH "*". A
recent paper (Judge, Low, Casini 2006, ApJ in press) suggests that,
for Fe XIII, one should use CECOEF=9.0e-09 with atom.fe13 to mimick
the depolarization seen in the much bigger calculation (atom.fe13o).
For the other models use CECOEF=0.

OUTPUT:
-------

For portability, output is to ASCII files 

OUT    
ATMOS (when IWATMO=1)

While these can be read, the IDL procedures in the idl/ directory 
enable one to look at the data very easily- you'll need the 
mapping software in SSW to make life very easy. 

After running the program check that the string "NORMAL END"
appears in the file 

JOBLOG

- because then the program completed successfully.

The files DUMI and DUMS are used during run time and can be deleted.


IDL procedures:
---------------

see the file AAAREADME in the idl/ subdirectory.

