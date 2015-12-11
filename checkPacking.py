#!/usr/bin/python

import sys
import numpy as np

file=sys.argv[1]
data=np.loadtxt(file, skiprows=5)

smallSep=[]
tol = 1.0e-5
for i in range(1,len(data)):
    for j in range(i):
        dist=np.sum((data[i,1:]-data[j,1:])**2)**0.5
        sep = dist - (data[i,0] + data[j,0])
        if sep < tol:
            smallSep.append([sep, i, j])
smallSep.sort()
print smallSep
