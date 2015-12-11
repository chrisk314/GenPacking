#!/usr/bin/python

import sys
import numpy as np
import numpy.random as npr
import math
import os
import time

InFile = 'packing.txt'
data = np.loadtxt(InFile, skiprows=5)
partRad = data[:,0]
partPos = data[:,1:]

numPlaced = len(data)

# Output data in format of MF-Unstructured AddHole command
packingFile = 'packing.mf.addhole'
outFile = open(packingFile, 'w')

domainBigger = 1.0e-09
tetMaxVol = domainBigger / 1.0e+05
tetMinAng = 15.0
tetRadEdge = 2.0
ScaleGeo = 1000.0

tol = 0.01 * partRad[-1]

outerBounds = np.zeros((3,2))
for i in range(numPlaced):
    if partPos[i,0] + partRad[i] > outerBounds[0,1]:
        outerBounds[0,1] = partPos[i,0] + partRad[i]
    if partPos[i,0] - partRad[i] < outerBounds[0,0]:
        outerBounds[0,0] = partPos[i,0] - partRad[i]

    if partPos[i,1] + partRad[i] > outerBounds[1,1]:
        outerBounds[1,1] = partPos[i,1] + partRad[i]
    if partPos[i,1] - partRad[i] < outerBounds[1,0]:
        outerBounds[1,0] = partPos[i,1] - partRad[i]
        
    if partPos[i,2] + partRad[i] > outerBounds[2,1]:
        outerBounds[2,1] = partPos[i,2] + partRad[i]
    if partPos[i,2] - partRad[i] < outerBounds[2,0]:
        outerBounds[2,0] = partPos[i,2] - partRad[i]

# Give some space at edges of domain
outerBounds[:,0] = outerBounds[:,0] - partRad[0]
outerBounds[:,1] = outerBounds[:,1] + partRad[0]

outFile.write("\tAutoMesh TET %+1.7e %f %f %+1.7e %+1.7e %+1.7e\n\n"\
        %(tetMaxVol * ScaleGeo**3, tetMinAng, tetRadEdge,\
        (outerBounds[0,1]-outerBounds[0,0] + 2*tol) * ScaleGeo,\
        (outerBounds[1,1]-outerBounds[1,0] + 2*tol) * ScaleGeo,\
        (outerBounds[2,1]-outerBounds[2,0] + 2*tol) * ScaleGeo))

# Determine refinement classes
refineMax = 3
Amax = math.pi * partRad[0]**2
Amin = math.pi * partRad[-1]**2
AmaxElem = Amax / (20 * 4**(refineMax - 1))
minDelta = 1.0E30

refineMin = refineMax
for k in range(refineMax,0,-1):
    AminElem = Amin / (20 * 4**(k - 1))
    if abs(AmaxElem - AminElem) < minDelta:
        minDelta = abs(AmaxElem - AminElem)
        refineMin = k

radClasses = [partRad[-1]]
if refineMin != refineMax:
    refineDelta = refineMax - refineMin
    bounds = 0
    for i in range(numPlaced-1, 0, -1):
        AElem0 = (math.pi * partRad[i]**2) / (20 * 4**(refineMin + bounds - 1))
        AElem1 = (math.pi * partRad[i]**2) / (20 * 4**(refineMin + bounds))
        if abs(AmaxElem - AElem0) > abs(AmaxElem - AElem1):
            radClasses.append(partRad[i])
            bounds += 1
        if bounds == refineDelta:
            radClasses.append(partRad[0])
            break

# Translate positions to place corner of cell at origin
AddHolePos = np.empty(partPos.shape)
AddHolePos[:,0] = ScaleGeo * (partPos[:,0] + abs(outerBounds[0,0]) + tol)
AddHolePos[:,1] = ScaleGeo * (partPos[:,1] + abs(outerBounds[1,0]) + tol)
AddHolePos[:,2] = ScaleGeo * (partPos[:,2] + abs(outerBounds[2,0]) + tol)
AddHoleRad = ScaleGeo * partRad

refine = refineMax
for i in range(numPlaced):
    for j in range(len(radClasses)-1):
        if partRad[i] > radClasses[j] and partRad[i] < radClasses[j+1]:
            refine = refineMin + j
    outFile.write("\tAddHole %+1.7e %+1.7e %+1.7e %+1.7e 1.0 1.0 1.0 0.0 0.0 0.0 %d\n"\
                 %(AddHolePos[i,0], AddHolePos[i,1], AddHolePos[i,2], AddHoleRad[i], refine))
outFile.write("\n\tScaleGeometryFactor %+1.7e\n" %(ScaleGeo))
outFile.close() 
