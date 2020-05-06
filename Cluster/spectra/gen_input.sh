#=============================================================================80
#                    Generate Input Files For Spectra Calculations
#==============================================================================!
#    Discussion:
#Make a new directory (suggested name=Number Basis Functions)
#Generate various subdirectories with varying alpha0 parameter
#==============================================================================!
#   Modified:
# 5 April 2020
#   Author:
# Shane Flynn
#==============================================================================!
mkdir 500
cd 500
d2=3;
potential=tip4p;
NG=500;
GH_order=4;
i_alpha=0.0;
inc=0.1;
mon_num=1;
counter=1;
Ndir=10;
while [ "$counter" -le "$Ndir" ]; do
  mkdir "$counter"
  cp ../*.xyz "$counter"
  cp ../grid.dat "$counter"
  cp ../run.sh "$counter"
  cd "$counter"
#==============================================================================!
#   single arrow makes a file, double arrow appends to an existing file
#==============================================================================!
  echo "$d2" > input
#==============================================================================!
#           can't use floating point with bash, use bc instead
#==============================================================================!
  echo "$potential" >> input
  echo "$NG" >> input
  echo "$GH_order" >> input
  echo "$i_alpha+ $counter*$inc" | bc  >> input
  echo "$mon_num" >> input
  cd ..
  counter=`expr "$counter" + 1`;
done
