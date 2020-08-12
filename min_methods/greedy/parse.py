#=============================================================================80
#                Parse Trajectory File Into Individual Frames
#=============================================================================80
#       Discussion:
#Assumes input file contains snapshots of a trajectory seperated by a blank line
#Reads every line until a new blank line is hit and dumps to a new cycle file
#==============================================================================#
#       Modified:
# 11 August 2020
#       Author:
# Shane Flynn
#==============================================================================#
from sys import exit, argv
script,fname=argv
#for example: python3 parse.py data.dat
#==============================================================================#
try:
    if(len(fname)<2):
        f=open('test.dat')
    else:
        f=open(fname)
except:
    print(f'Error opening file: {fname}: ==> exit')
    exit()
#==============================================================================#
file_counter=0
for line in f:
    if len(line.strip()) < 1:
        file_counter+=1
    with open('cycle' + str(file_counter) + '.dat', 'a+') as output:
        output.write(line)
