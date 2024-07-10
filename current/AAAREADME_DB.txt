DEC 29 2020

TO BUILD AND RUN THE DB (DATABASE) SINGLE POINT "INVERSIONS"

FIRST GET TO MAIN CLE DIRECTORY

cd cle

BUILD EXEC:

      cd current
      make db

(note that this version of 29 dec uses only one value of B=1G
but it reads DB.INPUT that still contains an older array version
for B that is not needed.  the file "db.hdr' does not contain 
the array, it is read by the python codes. 

BUILD DATABASE:

      cd ../seek/dbcle/test
      ../../../current/db          

RUN EXAMPLE FILE TO PRODUCE RESULTS AND PLOTS
    cd ../../
    python example.py


