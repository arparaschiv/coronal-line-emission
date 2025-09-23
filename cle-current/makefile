#
#  MAKEFILE FOR CLE
#
#******************************** PGI Fortran ************************
#FC      =       pgf77
#LINK      =       pgf77
#FFLAGS  =      -fast
#LDFLAGS =       -fast 
#DEBUG=       
#OPT= -fast
#******************************** g77  ************************
#FC=g77
#LINK=g77 
#DEBUG=       
#OPT= -O
#******************************** intel fortran  ************************
#FC=ifort
#LINK=ifort 
#DEBUG=       
#OPT= -fast
#LFLAGS=$(DEBUG)
#FFLAGS=$(OPT) $(DEBUG)
#********************************  gfortran  ************************
FC=gfortran
LINK=gfortran 
DEBUG=       
OPT= -ffixed-line-length-none -Wextra -o3 -static -static-libgcc -static-libstdc++ -static-libgfortran
LFLAGS=$(DEBUG)
FFLAGS=$(OPT) $(DEBUG)

LFLAGS=$(DEBUG)
FFLAGS=$(OPT) $(DEBUG)
.f.o:
	$(FC) $(FFLAGS) $*.f -c -o $*.o

CLE = cle.o m1synth.o atom.o atom1.o c2.o c3.o c5.o c6.o \
ce.o ci.o cs.o close.o colrd.o colcal.o convl.o corona.o sheet.o sheetx.o dvary.o \
cstrip.o cptime.o cuser.o dlamb.o eqsyst.o ffl.o pnl.o freq.o \
freqc.o freql.o getwrd.o grid.o gammae.o gammai.o gammas.o iatomn.o \
eqion.o ljust.o ltepop.o open.o planck.o profil.o \
rewind.o rinput.o rot.o sconst.o second_unix.o start.o stop.o watom.o \
t0tens.o emission.o fctsg.o fieldb.o field_int.o fn1.o \
fn2.o fun3j.o fun6j.o fun9j.o wdebug.o outs.o \
r1.o r2.o r3.o r4.o r5.o r6.o se0_build.o se0_copy.o \
solve.o trap.o uv.o splin.o dipole.o dipolv.o watmos.o zl.o

cle: $(CLE)
	$(LINK) $(CLE) $(LFLAGS) -o cle



DB = db.o dbe.o atom.o atom1.o c2.o c3.o c5.o c6.o \
ce.o ci.o cs.o close.o colrd.o colcal.o convl.o corona.o sheet.o sheetx.o dvary.o \
cstrip.o cptime.o cuser.o dlamb.o eqsyst.o ffl.o pnl.o freq.o \
freqc.o freql.o getwrd.o grid.o gammae.o gammai.o gammas.o iatomn.o \
eqion.o ljust.o ltepop.o open.o planck.o profil.o \
rewind.o rinput.o rot.o sconst.o second_unix.o start.o stop.o watom.o \
t0tens.o emission.o fctsg.o fieldb.o field_int.o fn1.o \
fn2.o fun3j.o fun6j.o fun9j.o wdebug.o outs.o \
r1.o r2.o r3.o r4.o r5.o r6.o se0_build.o se0_copy.o \
solve.o trap.o uv.o splin.o dipole.o dipolv.o watmos.o zl.o

db: $(DB)
	$(LINK) $(DB) $(FFLAGS) -o db

atom.o:PREC PARAM CATOM CATMOS CSLINE CCONST CINTS CLU
atom1.o:PREC PARAM CATOM CATMOS CSLINE CCONST CLU CSE CINTS
c2.o:PREC
c3.o:PREC PARAM CSE SGNM CINTS CINPUT CLU
c5.o:PREC PARAM CSE SGNM CINTS CINPUT CLU
c6.o:PREC
ce.o:PREC PARAM CSE SGNM CINTS
ci.o:PREC PARAM CSE SGNM CINTS
cle.o:PREC PARAM CGRID CCONST CBCNST 
pca.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CSLINE CORON CCONST CBCNST  CINPUT CINTS CLU
pcae.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CSLINE CORON CCONST CBCNST  CINPUT CINTS CLU
db.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CSLINE CORON CCONST CBCNST  CINPUT CINTS CLU
dbe.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CSLINE CORON CCONST CBCNST  CINPUT CINTS CLU
close.o:PREC COPCL
colcal.o:PREC PARAM COLL CATOM CATMOS CATMO2 CCONST CINPUT CLU
colrd.o:PREC PARAM COLL CATOM CATMOS CATMO2 CCONST CINPUT CLU
convl.o:PREC
corona.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CCONST CORON CLU
dvary.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CCONST CORON CLU
sheet.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CCONST CORON CLU
sheetx.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CCONST CORON CLU
cptime.o:PREC CLU
cs.o:PREC PARAM CSE CLU SGNM CINTS
cstrip.o:PREC CLU
cuser.o:CCONST CGRID CINPUT CLU
dipole.o:PREC CCONST CLU CGRID
dipolv.o:PREC CCONST CLU CGRID
dlamb.o:PREC PARAM CATOM CCONST
emission.o:PREC PARAM CSE CSLINE CCONST CATOM CORON SGNM CINTS CLU
eqsyst.o:PREC
fctsg.o:PREC PARAM SGNM CINTS FACT
field_int.o:PREC CATOM CORON CCONST CINTS
fieldb.o:PREC PARAM CATOM CCONST CORON CSE
ffl.o:PREC PARAM CCONST CLU CGRID CINPUT CINTS 
pnl.o:PREC CINTS
fn1.o:PREC PARAM FACT CINTS
fn2.o:PREC PARAM FACT
freq.o:PREC PARAM CATOM CCONST CSLINE 
freqc.o:PREC PARAM CATOM CCONST CSLINE 
freql.o:PREC PARAM CSLINE 
fun3j.o:PREC PARAM SGNM CINTS FACT
fun6j.o:PREC PARAM SGNM CINTS FACT
fun9j.o:PREC PARAM SGNM CINTS
gammae.o:PREC PARAM CSE CLU CATMOS CATOM CINPUT CINTS 
gammai.o:PREC PARAM CSE CLU CATMOS CATOM CINPUT CINTS 
gammas.o:PREC PARAM CSE CLU CATMOS CATOM CINPUT CINTS 
getwrd.o:PREC 
grid.o:PREC GRDPAR CGRID CINPUT CLU
iatomn.o:
eqion.o:PREC PARAM CATOM CLU
ljust.o:PREC
ltepop.o:PREC PARAM CATOM CATMOS CCONST CLU
m1synth.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CSLINE CORON CCONST CINPUT CINTS CLU
open.o:PREC CLU COPCL CBCNST
wdebug.o:PREC PARAM CSE CATOM CSLINE CLU CINPUT CORON CGRID
outs.o:PREC PARAM CATOM CSLINE CINPUT CGRID CCONST CLU
planck.o:PREC CCONST
profil.o:PREC PARAM CATOM CATMOS CSLINE CLU CCONST
r1.o:PREC PARAM CSE CINTS
r2.o:PREC PARAM CSE SGNM CATOM CINTS
r3.o:PREC PARAM CSE SGNM CINTS
r4.o:PREC PARAM CSE SGNM CINTS
r5.o:PREC PARAM CSE SGNM CINTS
r6.o:PREC PARAM CSE SGNM CATOM CINTS
rewind.o:PREC
rinput.o:PREC PARAM CATOM CSLINE CINPUT CLU 
rot.o:PREC
sconst.o:PREC PARAM CCONST CBCNST CINTS
se0_build.o:PREC PARAM CATOM CSE CLU CINTS
se0_copy.o:PREC PARAM CSE CLU CINTS
second_unix.o:PREC
solve.o:PREC PARAM CSE CLU CINPUT CINTS
splin.o:PREC
start.o:PREC PARAM CATOM CATMOS CSLINE CCONST CINPUT CLU COPCL
stop.o:PREC CLU COPCL
t0tens.o:PREC PARAM CORON CCONST SGNM CINTS CLU
trap.o:PREC PARAM CINPUT CSLINE CATOM
uv.o:PREC PARAM 
watmos.o:PREC PARAM CINPUT CGRID CATMOS CATMO2 CATOM CLU 
watom.o:PREC PARAM CATOM CCONST CSLINE CLU CINPUT
zl.o:PREC PARAM CCONST CGRID CINPUT CLU CINTS



PGJTEST = pgjtest.o m1synth.o atom.o atom1.o c2.o c3.o c5.o c6.o \
ce.o ci.o cs.o close.o colrd.o colcal.o convl.o corona.o sheet.o sheetx.o\
cstrip.o cptime.o cuser.o dlamb.o eqsyst.o ffl.o pnl.o freq.o \
freqc.o freql.o getwrd.o grid.o gammae.o gammai.o gammas.o iatomn.o \
eqion.o ljust.o ltepop.o open.o planck.o profil.o \
rewind.o rinput.o rot.o sconst.o second_unix.o start.o stop.o watom.o \
t0tens.o emission.o fctsg.o fieldb.o field_int.o fn1.o \
fn2.o fun3j.o fun6j.o fun9j.o wdebug.o outs.o \
r1.o r2.o r3.o r4.o r5.o r6.o se0_build.o se0_copy.o \
solve.o trap.o uv.o splin.o dipole.o watmos.o zl.o


pgjtest: $(PGJTEST)
	$(LINK) $(PGJTEST) $(LFLAGS) -o pgjtest

atom.o:PREC PARAM CATOM CATMOS CSLINE CCONST CINTS CLU
atom1.o:PREC PARAM CATOM CATMOS CSLINE CCONST CLU CSE CINTS
c2.o:PREC
c3.o:PREC PARAM CSE SGNM CINTS CINPUT CLU
c5.o:PREC PARAM CSE SGNM CINTS CINPUT CLU
c6.o:PREC
ce.o:PREC PARAM CSE SGNM CINTS
ci.o:PREC PARAM CSE SGNM CINTS
pgjtest.o:PREC PARAM CGRID CCONST CBCNST 
close.o:PREC COPCL
colcal.o:PREC PARAM COLL CATOM CATMOS CATMO2 CCONST CINPUT CLU
colrd.o:PREC PARAM COLL CATOM CATMOS CATMO2 CCONST CINPUT CLU
convl.o:PREC
corona.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CCONST CORON CLU
sheet.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CCONST CORON CLU
sheetx.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CCONST CORON CLU
cptime.o:PREC CLU
cs.o:PREC PARAM CSE CLU SGNM CINTS
cstrip.o:PREC CLU
cuser.o:CCONST CGRID CINPUT CLU
dipole.o:PREC CCONST CLU CGRID
dlamb.o:PREC PARAM CATOM CCONST
emission.o:PREC PARAM CSE CSLINE CCONST CATOM CORON SGNM CINTS CLU
eqsyst.o:PREC
fctsg.o:PREC PARAM SGNM CINTS FACT
field_int.o:PREC CATOM CORON CCONST CINTS
fieldb.o:PREC PARAM CATOM CCONST CORON CSE
ffl.o:PREC PARAM CCONST CLU CGRID CINPUT CINTS 
pnl.o:PREC CINTS
fn1.o:PREC PARAM FACT CINTS
fn2.o:PREC PARAM FACT
freq.o:PREC PARAM CATOM CCONST CSLINE 
freqc.o:PREC PARAM CATOM CCONST CSLINE 
freql.o:PREC PARAM CSLINE 
fun3j.o:PREC PARAM SGNM CINTS FACT
fun6j.o:PREC PARAM SGNM CINTS FACT
fun9j.o:PREC PARAM SGNM CINTS
gammae.o:PREC PARAM CSE CLU CATMOS CATOM CINPUT CINTS 
gammai.o:PREC PARAM CSE CLU CATMOS CATOM CINPUT CINTS 
gammas.o:PREC PARAM CSE CLU CATMOS CATOM CINPUT CINTS 
getwrd.o:PREC 
grid.o:PREC GRDPAR CGRID CINPUT CLU
iatomn.o:
eqion.o:PREC PARAM CATOM CLU
ljust.o:PREC
ltepop.o:PREC PARAM CATOM CATMOS CCONST CLU
m1synth.o:PREC PARAM CGRID CATOM CATMOS CATMO2 CSLINE CORON CCONST CINPUT CINTS CLU
open.o:PREC CLU COPCL CBCNST
wdebug.o:PREC PARAM CSE CATOM CSLINE CLU CINPUT CORON CGRID
outs.o:PREC PARAM CATOM CSLINE CINPUT CGRID CCONST CLU
planck.o:PREC CCONST
profil.o:PREC PARAM CATOM CATMOS CSLINE CLU CCONST
r1.o:PREC PARAM CSE CINTS
r2.o:PREC PARAM CSE SGNM CATOM CINTS
r3.o:PREC PARAM CSE SGNM CINTS
r4.o:PREC PARAM CSE SGNM CINTS
r5.o:PREC PARAM CSE SGNM CINTS
r6.o:PREC PARAM CSE SGNM CATOM CINTS
rewind.o:PREC
rinput.o:PREC PARAM CATOM CSLINE CINPUT CLU 
rot.o:PREC
sconst.o:PREC PARAM CCONST CBCNST CINTS
se0_build.o:PREC PARAM CATOM CSE CLU CINTS
se0_copy.o:PREC PARAM CSE CLU CINTS
second_unix.o:PREC
solve.o:PREC PARAM CSE CLU CINPUT CINTS
splin.o:PREC
start.o:PREC PARAM CATOM CATMOS CSLINE CCONST CINPUT CLU COPCL
stop.o:PREC CLU COPCL
t0tens.o:PREC PARAM CORON CCONST SGNM CINTS CLU
trap.o:PREC PARAM CINPUT CSLINE CATOM
uv.o:PREC PARAM 
watmos.o:PREC PARAM CINPUT CGRID CATMOS CATMO2 CATOM CLU 
watom.o:PREC PARAM CATOM CCONST CSLINE CLU CINPUT
zl.o:PREC PARAM CCONST CGRID CINPUT CLU CINTS