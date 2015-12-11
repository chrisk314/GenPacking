#! /usr/bin/python

import sys
import numpy as np
import os

# Get current script name.
thisScript = os.path.basename(__file__)

# Take as argument the number of particles in the x, y and z directions and radius
if len(sys.argv) != 6:
    print "%s: Usage: %s NumX NumY NumZ Radius RowSep"%(thisScript, thisScript)
    exit(1)
else:
    NumX = int(sys.argv[1])   # Number of particles in x dir.
    NumY = int(sys.argv[2])   # Number of particles in y dir.
    NumZ = int(sys.argv[3])   # Number of particles in z dir.
    Radius = float(sys.argv[4]) # Radius of particles
    RowSep = float(sys.argv[5])	
    
    if NumZ % 2 != 0:
        print "%s: NumZ must be even."%(thisScript)
        exit(1)
    
    if RowSep < 2*Radius:
        print "%s: Require RowSep >= 2*Radius."%(thisScript)
            
	
# Set origin of coordinate system
Origin = (0.0, 0.0, 0.0)

# Set material parameters for particle type
Density = 10000
YoungsMod = 7.001e10
Poisson = 0.2
Friction = 0.0
Alpha = 0.0

# Allocate array for particle coordinates
TotalParticles = NumX * NumY * NumZ
partPos = np.zeros((TotalParticles,3))

# Create first plane of particles
for i in range(NumX):
    for j in range(NumY):
        partPos[i*NumY + j,0] = Origin[0] + (i + 0.5) * RowSep
        partPos[i*NumY + j,1] = Origin[1] + (j + 0.5) * RowSep
        partPos[i*NumY + j,2] = Origin[2] + 0.5 * 2.0**-0.5 * RowSep
		
# Now repeat for subsequent planes
for k in range(1,NumZ):
    for i in range(NumX*NumY):
	# Shift odd planes into depression of plane beneath
	if k % 2 != 0:
            partPos[k*NumX*NumY + i][0] = partPos[i][0] + 0.5 * RowSep
            partPos[k*NumX*NumY + i][1] = partPos[i][1] + 0.5 * RowSep
	# Leave even planes unaltered
	else:
            partPos[k*NumX*NumY + i][0] = partPos[i][0]
            partPos[k*NumX*NumY + i][1] = partPos[i][1]
		
	# Set z coord	
	partPos[k*NumX*NumY + i][2] = Origin[2] + k * 2.0**-0.5 * RowSep

# Calculate coordinates of the boundaries
Lx = NumX * RowSep
Ly = NumY * RowSep
Lz = 2**-0.5 * (NumZ + 1) * RowSep

# Write data to file			
packingFile = "packings/packing_FCC_%d.txt"%(TotalParticles)
outFile = open(packingFile,'w')

outFile.write("NUMPARTS\n%d\n"%TotalParticles)
outFile.write("BOXDIMS\n%+1.15e %+1.15e %+1.15e\n"%(Lx, Ly, Lz))
outFile.write("PARTDATA\n")
for i in range(TotalParticles):
  outFile.write("%+1.15e %+1.15e %+1.15e %+1.15e\n"\
                %(Radius, partPos[i,0], partPos[i,1], partPos[i,2]))        
outFile.close()


# Output data in format of MF-Unstructured AddHole command. --------------------
packingFile = 'packings/packing_FCC_%d.mf.addhole'%(TotalParticles)
outFile = open(packingFile, 'w')

tetMaxVol = (Lx * Ly * Lz) / 1.0e+05
tetMinAng = 15.0
tetRadEdge = 2.0
ScaleGeo = 1.0

outFile.write("\tAutoMesh TET %+1.7e %f %f %+1.7e %+1.7e %+1.7e\n\n"\
        %(tetMaxVol * ScaleGeo**3, tetMinAng, tetRadEdge, Lx * ScaleGeo,\
        Ly * ScaleGeo, Lz * ScaleGeo))

outFile.write("\tAddHoleCount %d\n\n" %(TotalParticles))
refine = 3
for i in range(TotalParticles):
    outFile.write("\tAddHole %+1.7e %+1.7e %+1.7e %+1.7e 1.0 1.0 1.0 0.0 0.0 0.0 %d\n"\
        %(partPos[i,0], partPos[i,1], partPos[i,2], Radius, refine))
outFile.write("\n\tScaleGeometryFactor %+1.7e\n" %(ScaleGeo))
outFile.close() 

# Output data in format of MF-Unstrucutred IB command
packingFile = 'packings/packing_FCC_%d.mf.ib'%(TotalParticles)
outFile = open(packingFile, 'w')

outFile.write("\tAutoMesh HEX 100 100 100 %1.7e %1.7e %1.7e\n\n"%(Lx, Ly, Lz))
outFile.write("IMMERSEDBOUNDARY\n\tCOUNT %d\n"%(TotalParticles))
refine = 3
for i in range(TotalParticles):
    outFile.write("\tOBJECT %d\n\t\tTYPE MFTL\n\t\tSHAPE SPHERE\n\t\tRADIUS %1.7e\n"\
        "\t\tREFINELEVEL %d\n\t\tINITCENTER %1.7e %1.7e %1.7e\n\tENDOBJECT\n"\
        %(i, Radius, refine, partPos[i,0], partPos[i,1], partPos[i,2]))
outFile.write("END\n")
outFile.close()
