#!/usr/bin/python

import sys
import numpy as np
import numpy.random as npr
import math
import os
import time

# Get current script name.
thisScript = os.path.basename(__file__)

# ------------------------------------------------------------------------------
# Function definitions
# ------------------------------------------------------------------------------
def volFrac(r,a,p):
    if p == 0:
        return phiFCC(r,a)
    if p == 1:
        return phiBCC(r,a)
    if p == 2:
        return phiSC(r,a)

def phiFCC(r,a):
    return 4.0 * (4.0 * math.pi / 3.0) * r**2

def phiBCC(r,a):
    return 2.0 * (4.0 * math.pi / 3.0) * r**2

def phiSC(r,a):
    return (4.0 * math.pi / 3.0) * r**2

# Generate primitive cell positions
def primCell(r,a,p):
    posData = np.empty((14,3))

    # Allocate SC positions
    posData[0] = [0.0, 0.0, 0.0]
    posData[1] = [a, 0.0, 0.0]
    posData[2] = [0.0, a, 0.0]
    posData[3] = [a, a, 0.0]
    posData[4] = [0.0, 0.0, a]
    posData[5] = [a, 0.0, a]
    posData[6] = [0.0, a, a]
    posData[7] = [a, a, a]

    if p == 0:
        posData[8] = [0.5*a, 0.5*a, 0.0]
        posData[9] = [0.5*a, 0.0, 0.5*a]
        posData[10] = [0.0, 0.5*a, 0.5*a]
        posData[11] = [0.5*a, a, 0.5*a]
        posData[12] = [a, 0.5*a, 0.5*a]
        posData[13] = [0.5*a, 0.5*a, a]
        return posData

    elif p == 1:
        posData[8] = [0.5*a, 0.5*a, 0.5*a]
        return posData[:9,:]

    else:
        return posData[:8,:]

def uniqueRows(a):
  a = np.ascontiguousarray(a)
  unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
  return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

# ------------------------------------------------------------------------------
# Program inputs
# ------------------------------------------------------------------------------
# default values.
r = 0.01            # Particle radius.
#a = 2 * 2**0.5 * r  # Interlayer spacing.
a = 4.5 * r
p = 0               # Packing type: 0 FCC, 1 BCC, 2 SC.
n = 14              # Target number of particles. 

# Read values from command line.
if len(sys.argv) != 5:
  print("%s: Usage: GenArray.py <r> <a> <p> <n>"%thisScript)
  print("%s: Using default options."%thisScript)
else:
  r = float(sys.argv[1])
  a = float(sys.argv[2])
  p = int(sys.argv[3])
  n = int(sys.argv[4])
  density = 1.0 # Required for IB

# Check inputs -----------------------------------------------------------------
if p == 0:
    if a < 2.0 * 2.0**0.5 * r:
        print "Invalid value of a. a E [2*2**0.5*r, inf] for FCC"
        print "Program will exit."
        exit(1)
    if n < 14:
        print "Invalid value of n. n E [14, inf] for FCC"
        print "Program will exit."
        exit(1)

elif p == 1:
    if a < (4.0/3.0**0.5) * r:
        print "Invalid value of a. a E [(4/3**0.5)*r, inf] for BCC"
        print "Program will exit."
        exit(1)
    if n < 9:
        print "Invalid value of n. n E [9, inf] for BCC"
        print "Program will exit."
        exit(1)
        
elif p == 2:
    if a < 2.0 * r:
        print "Invalid value of a. a E [2**r, inf] for SC"
        print "Program will exit."
        exit(1)
    if n < 8:
        print "Invalid value of n. n E [8, inf] for SCC"
        print "Program will exit."
        exit(1)

else:
    print "Invalid value of p: 0 FCC, 1 BCC, 2 SC."
    print "Program will exit."
    exit(1)

# Create primitive cell --------------------------------------------------------
primCellPos = primCell(r,a,p)
partPos = primCellPos
#config = [8,8,8]
config = [2,1,1]
#space = 1.0e-07
space = 0.0
Lx = a * (config[0] + 1) + r
Ly = a * (config[1] + 1) + space
Lz = a * (config[2] + 1) + space

for i in range(config[0]):
    for j in range(config[1]):
        for k in range(config[2]):
            transPos = primCellPos + [i*a, j*a, k*a]
            partPos = np.vstack((partPos, transPos))

# Eliminate duplicate positions
partPos = uniqueRows(partPos)
# Translate particles off origin
partPos = partPos + [0.5 * (a + r), 0.5 * (a + space), 0.5 * (a + space)]
#partPos = partPos + [0.5 * (a + r), 0.0 * (a + space), 0.0 * (a + space)]
numParts = partPos.shape[0]

# Calculate porosity
porosity = 1.0 - (numParts * 4.0 * math.pi * r**3.0 / 3.0 ) / (Lx * Ly * Lz)

# Write out data ---------------------------------------------------------------
packingFile = "packings/packing_p%d_%d.txt"%(p,numParts)
outFile = open(packingFile,'w')

outFile.write("NUMPARTS\n%d\n"%numParts)
outFile.write("POROSITY\n%1.4f\n"%porosity)
outFile.write("BOXDIMS\n%+1.15e %+1.15e %+1.15e\n"%(Lx, Ly, Lz))
outFile.write("PARTDATA\n")
for i in range(numParts):
  outFile.write("%+1.15e %+1.15e %+1.15e %+1.15e\n"\
                %(r, partPos[i,0], partPos[i,1], partPos[i,2]))        
outFile.close()


# Output data in .vtp format readable by ParaView. -----------------------------
packingFile = "packings/packing_p%d_%d.vtp"%(p,numParts)
outFile = open(packingFile, 'w')
outFile.write("<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" format=\"ascii\">\n")
outFile.write("<PolyData>\n\t<Piece NumberOfPoints=\"%d\">\n\t\t<Points>\n\t\t\t"%(numParts))
outFile.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
for i in range(numParts):
  outFile.write("\t\t\t\t%f %f %f\n"%(partPos[i,0], partPos[i,1], partPos[i,2]))
outFile.write("\t\t\t</DataArray>\n\t\t</Points>\n")
outFile.write("\t\t<PointData Scalars=\"Diameter\">\n\t\t\t<DataArray type=\"Float32\" Name=\"Diameter\" format=\"ascii\">\n")
for i in range(numParts):
  outFile.write("\t\t\t\t%f\n"%(2*r))
outFile.write("\t\t\t\t</DataArray>\n")
outFile.write("\t\t\t<DataArray type=\"Int32\" Name=\"ID\" format=\"ascii\" NumberOfComponents=\"1\">\n")
for i in range(numParts):
  outFile.write("\t\t\t\t%d\n"%(i))
outFile.write("\t\t\t\t</DataArray>\n")
outFile.write("\t\t\t</PointData>\n\t\t</Piece>\n\t</PolyData>\n</VTKFile>")
outFile.close()

# Output data in format of MF-Unstructured AddHole command. --------------------
packingFile = 'packings/packing_p%d_%d.mf.addhole'%(p,numParts)
outFile = open(packingFile, 'w')

tetMaxVol = (Lx * Ly * Lz) / 1.0e+05
tetMinAng = 15.0
tetRadEdge = 2.0
ScaleGeo = 1.0e+02

outFile.write("\tAutoMesh TET %+1.7e %f %f %+1.7e %+1.7e %+1.7e\n\n"\
        %(tetMaxVol * ScaleGeo**3, tetMinAng, tetRadEdge, Lx * ScaleGeo, Ly * ScaleGeo,\
        Lz * ScaleGeo))

outFile.write("\tAddHoleCount %d\n\n" %(numParts))
refine = 3
for i in range(numParts):
    outFile.write("\tAddHole %+1.7e %+1.7e %+1.7e %+1.7e 1.0 1.0 1.0 0.0 0.0 0.0 %d\n"\
        %(partPos[i,0] * ScaleGeo, partPos[i,1] * ScaleGeo, partPos[i,2] * ScaleGeo, r * ScaleGeo, refine))
outFile.write("\n\tScaleGeometryFactor %+1.7e\n" %(1.0/ScaleGeo))
outFile.close() 

# Output data in format of MF-Unstrucutred IB command
packingFile = 'packings/packing_p%d_%d.mf.ib'%(p,numParts)
outFile = open(packingFile, 'w')

outFile.write("\tAutoMesh HEX 100 100 100 %1.7e %1.7e %1.7e\n\n"%(Lx, Ly, Lz))
outFile.write("IMMERSEDBOUNDARY\n\tCOUNT %d\n"%(numParts))
refine = 3
for i in range(numParts):
    outFile.write("\tOBJECT %d\n\t\tDENSITY %+1.7e\n\t\tTYPE MFTL\n\t\tMOVE NO\n\t\tSHAPE SPHERE\n\t\tRADIUS %+1.7e\n"\
        "\t\tREFINELEVEL %d\n\t\tINITCENTER %+1.7e %+1.7e %+1.7e\n\tENDOBJECT\n"\
        %(i, density, r, refine, partPos[i,0], partPos[i,1], partPos[i,2]))
outFile.write("END\n")
outFile.close()
