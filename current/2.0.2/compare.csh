#!/bin/csh
# 
#
rm differences
touch differences
foreach f (*.f)
  echo 'doing ' $f
  diff $f ../ORIG/current/$f >> differences
end
#

#
