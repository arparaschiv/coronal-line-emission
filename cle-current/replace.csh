#!/bin/csh
# 
#
foreach f (*.f)
  echo 'doing ' $f
  if(-e a.jou) \rm a.jou
  cat $f | sed -e '/LOUT1/s/LOUT1/LOUT/' \
| sed -e '/LOUT1/s/LOUT1/LOUT/' \
| sed -e '/LOUT1/s/LOUT1/LOUT/'  > a.jou
  mv a.jou $f
end
#

#
