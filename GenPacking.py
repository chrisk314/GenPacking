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
# Takes the position of a particle. Returns indices of containing subdomains.
def getCubeIndex(pos):
  cubeIndex = []
  subCubes = []
  for i in range(lCube):
    if pos[2] > cubeData[i*lCube2].bounds[4] and pos[2] < cubeData[i*lCube2].bounds[5]:
      subCubes = [i*lCube2 + x for x in range(lCube2)]
      if not i == lCube - 1:
        if pos[2] > cubeData[(i+1)*lCube2].bounds[4]:
          subCubes.extend([(i+1)*lCube2 + x for x in range(lCube2)])
      break
  for i in subCubes:
    if (pos[0] > cubeData[i].bounds[0] and pos[0] < cubeData[i].bounds[1])\
          and (pos[1] > cubeData[i].bounds[2] and pos[1] < cubeData[i].bounds[3]):
      cubeIndex.append(i)
  return cubeIndex

# Generate random position.
def genPos():
  return [npr.rand() * Lx, npr.rand() * Ly, npr.rand() * Lz] 

# Check if particle crosses periodic boundaries.
def checkPbc(pos, pbc):
  check = False
  if pos[0] < partRad[0]:
    pbc[0] = 1;
    check = True
  elif pos[0] > Lx - partRad[0]:
    pbc[0] = -1
    check = True
  if pos[1] < partRad[0]:
    pbc[1] = 1;
    check = True
  elif pos[1] > Ly - partRad[0]:
    pbc[1] = -1
    check = True
  if pos[2] < partRad[0]:
    pbc[2] = 1;
    check = True
  elif pos[2] > Lz - partRad[0]:
    pbc[2] = -1
    check = True
  return check
      
def uniqueRows(a):
  a = np.ascontiguousarray(a)
  unique_a = np.unique(a.view([('', a.dtype)]*a.shape[1]))
  return unique_a.view(a.dtype).reshape((unique_a.shape[0], a.shape[1]))

# ------------------------------------------------------------------------------
# Program inputs
# ------------------------------------------------------------------------------
# default values.
gradFile = "Grading_ToyouraSand_q0.txt"
poros = 0.4
Lx = 1.0e-3
Ly = 1.0e-3
Lz = 1.0e-3
minDia = 1.2e-4
volExcess = 0.1
#overlapTol = 1.0e-6
overlapTol = 0.0

# Read values from command line.
if len(sys.argv) != 9:
  print("%s: Usage: GenPacking.py <grading file> <porosity> <Lx> <Ly> <Lz>"\
    " <min. dia.> <vol. excess> <overlap tol.>"%thisScript)
  print("%s: Using default options."%thisScript)
else:
  gradFile = sys.argv[1]
  poros = float(sys.argv[2])
  Lx = float(sys.argv[3])
  Ly = float(sys.argv[4])
  Lz = float(sys.argv[5])
  minDia = float(sys.argv[6])
  volExcess = float(sys.argv[7])
  overlapTol = float(sys.argv[8])  

domainVol = Lx * Ly * Lz
domainPartVol = domainVol * (1 - poros)
# Give some extra space for the particles.
domainBigger = domainVol * (1.0 + volExcess)
lengthExcess = (1.0 + volExcess)**(1.0/3.0)
Lx *= lengthExcess
Ly *= lengthExcess
Lz *= lengthExcess

# ------------------------------------------------------------------------------
# Handle PSD grading file and generate particle diameters
# ------------------------------------------------------------------------------
gradData = np.genfromtxt(gradFile, usecols=(0,1,2))
gradData[:,:2] = np.divide(gradData[:,:2], 1000) # Convert mm to m.
gradData[:,2] = np.divide(gradData[:,2], 100) # Convert % to prob.

# Eliminate small particle classes.
for gradClass in gradData:
  if gradClass[1] < minDia:
    gradClass[2] = 0.0
  elif gradClass[1] > minDia and gradClass[0] < minDia:
    gradClass[2] *= (gradClass[1] - minDia) / (gradClass[1] - gradClass[0])
    gradClass[0] = minDia

meanDia = 0.5 * (gradData[:,0] + gradData[:,1])             # Class mean particle diameters.
classPartVol = (math.pi/6) * meanDia[:]**3 * gradData[:,2]  # Class mean particle volumes.
meanPartVol = sum(classPartVol)                             # Global mean particle volume.
numPartsEst = domainPartVol / meanPartVol                   # Estimated number of particles.
classVol = numPartsEst * classPartVol                       # Class volumes.

# Generate particle diameters.
partDia = []
overflowVol = 0.0
for i in reversed(range(len(gradData))):
  accumVol = 0.0
  if classVol[i] > 0.0:
    while accumVol < classVol[i] + overflowVol:
      partDia.append(gradData[i,0] + npr.rand() * (gradData[i,1] - gradData[i,0]))
      addVol = (math.pi/6) * partDia[-1]**3
      accumVol += addVol
      if accumVol - classVol[i] > classVol[i] - (accumVol - addVol):
        partDia.pop()
        overflowVol += classVol[i] - (accumVol - addVol)

# Sort diameters from largest to smallest.
numPartsAttempt = len(partDia)
partDia.sort()
partDia.reverse()
partDia = np.array(partDia)
partRad = 0.5*partDia
partData = np.zeros((numPartsAttempt,4))

volActual = (math.pi/6) * sum(partDia**3)

print("\n%s: %d particles generated with volume %+1.4e within %+1.4e of target volume\n" \
      %(thisScript, numPartsAttempt, volActual, abs(domainPartVol - volActual)/domainPartVol))

# ------------------------------------------------------------------------------
# Split domain into subdomains to reduce neighbour search overhead
# ------------------------------------------------------------------------------
lCube = int(math.ceil((float(numPartsAttempt)/100)**(1.0/3.0)))
while lCube > 1:
  Sx = Lx/lCube
  if (Sx+partDia[0])/(2*partDia[0]) < 3:
    lCube -= 1
      
# Subdomain lengths.
Sx = Lx/lCube
Sy = Ly/lCube
Sz = Lz/lCube
  
lCube2 = lCube**2
numCubes = lCube**3

tol = 1.0e-8*partDia[0]

numPartsCubeEst = math.floor((1.4 * numPartsAttempt) / numCubes)
if numPartsCubeEst > numPartsAttempt:
  numPartsCubeEst = numPartsAttempt

# Cube data structure storing boundaries of subdomains and particle lists.
class Cube:
  # Initialise subdomain bounds and allocate memory.
  def __init__(self, index):
    self.partData = np.empty((numPartsCubeEst,4))
    self.bounds = np.empty(6)
    self.numParts = 0
    # Allow extra space or some overlap between particles by specifying
    # a positive or negative value of overlapTol respectively.
    self.overlapTol = overlapTol
    # Specify cube boundaries.
    self.bounds[0] = (index % lCube) * Sx - (partRad[0] + tol);
    self.bounds[2] = (math.floor(index/lCube) % lCube) * Sy - (partRad[0] + tol);
    self.bounds[4] = (math.floor(index/lCube2)) * Sz - (partRad[0] + tol);
    self.bounds[1] = self.bounds[0] + Sx + partDia[0] + tol;
    self.bounds[3] = self.bounds[2] + Sy + partDia[0] + tol;
    self.bounds[5] = self.bounds[4] + Sz + partDia[0] + tol;

  # Add particle to subdomain.
  def addPart(self, data):
    if self.numParts == self.partData.shape[0]:
      self.partData = np.vstack((self.partData, data))
    else:
      self.partData[self.numParts] = data
    self.numParts += 1

  # Check for overlaps between last added cube and other members of subdomain.
  def checkOverlaps(self):
    if self.numParts < 2:
      return False
    else:
      dist = np.sum((self.partData[:self.numParts-1,1:] - \
                     self.partData[self.numParts-1,1:])**2, axis=1)
      if np.min(dist - (self.partData[:self.numParts-1,0] + \
                        self.partData[self.numParts-1,0] + self.overlapTol)**2) < 0.0:
        return True
      return False
    
  def checkOverlapsAlt(self):
    if self.numParts < 2:
      return False
    else:
      dist = np.sum((self.partData[:self.numParts-1,1:] - \
                     self.partData[self.numParts-1,1:])**2, axis=1)
      if np.min(dist - (self.partData[:self.numParts-1,0] + \
                        self.partData[self.numParts-1,0] + self.overlapTol)**2) < 0.0:
        overlaps = np.abs(dist - (self.partData[:self.numParts-1,0] + \
                        self.partData[self.numParts-1,0] + self.overlapTol)**2)
        testRads = np.empty((self.numParts-1,2))
        testRads[:,0] = self.partData[:self.numParts-1,0]**2
        testRads[:,1] = self.partData[self.numParts-1,0]**2
        minRads = np.min(testRads, axis=1)
        overlapRatio = np.divide(minRads, overlaps)
        if np.any(np.greater(overlapRatio, 2.0)):
          return True
      return False

# Initialise subdomains.
cubeData = [Cube(i) for i in range(numCubes)]

# DEBUG
#file = open('cube_bounds.txt','w')
#for i, cube in enumerate(cubeData):
#  file.write("%d %f %f %f %f %f %f\n"%(i, cube.bounds[0], cube.bounds[1],\
#    cube.bounds[2], cube.bounds[3], cube.bounds[4], cube.bounds[5]))
#file.close()

# All possible ghost particle translation unit vectors.
pbcPerms = []
for i in range(2):
  for j in range(2):
    for k in range(2):
      pbcPerms.append([i, j, k])

posPerms = np.empty((8,3))
maxTries = 1.0e6
numPlaced = 0

# ------------------------------------------------------------------------------
# Main placement loop
# ------------------------------------------------------------------------------
for i in range(numPartsAttempt):
  tries = 0
  retry = True
  
  # DEBUG
  tStartClock = time.clock()
  tStartWall = time.time()
  
  # Attempt to place particle. Loop until success.
  while retry and tries < maxTries:
    tries += 1
    pos = genPos()            # Generate a random position.
    cubeIndex = getCubeIndex(pos)   # Determine subdomain membership.
    tryIndex = []       # List of all cubes in which placement has been attempted.
    
    # Check for overlaps with particles in containing subdomains.
    for j in cubeIndex:
      tryIndex.append(j)
      cubeData[j].addPart([partRad[i], pos[0], pos[1], pos[2]])
      retry = cubeData[j].checkOverlaps()
      if retry:
        # Remove particle from all cubes in which it has been placed.
        for k in tryIndex:
          cubeData[k].numParts -= 1
        break
      
    # Determine if particle crosses PBCs.
    if not retry:
      pbc = np.zeros(3)
      if checkPbc(pos, pbc):
        # Determine ghost translations and generate all possible ghost particles.
        posPerms = pos + np.multiply(np.multiply(pbcPerms, pbc), [Lx, Ly, Lz])
        # Keep only unique ghost positions and remove real particle from ghost list.
        posPerms = uniqueRows(posPerms)
        posPerms = np.delete(posPerms, np.where((posPerms == pos).all(axis=1)), axis=0) 
          
        # Check for overlaps in ghost positions.
        for ghost in posPerms:
          cubeIndex = getCubeIndex(ghost)
          for j in cubeIndex:
            tryIndex.append(j) 
            cubeData[j].addPart([partRad[i], ghost[0], ghost[1], ghost[2]])
            retry = cubeData[j].checkOverlaps()
            if retry:
              # Remove particle from all cubes in which it has been placed.
              for k in tryIndex:
                cubeData[k].numParts -= 1
              break
          if retry:
            break
    
    # Succesful placement.
    if not retry:
      partData[numPlaced,0] = partRad[i]
      partData[numPlaced,1:] = pos
      numPlaced += 1
      # DEBUG
      tEndClock = time.clock()
      tEndWall = time.time()
      timeTotalClock = tEndClock - tStartClock
      timeTotalWall = tEndWall - tStartWall
      timePerIterClock = (timeTotalClock)/tries 
      timePerIterWall = (timeTotalWall)/tries 
      print("%s: Placed particle %d at %f %f %f in %d attempts"\
            %(thisScript, i, pos[0], pos[1], pos[2], tries))
      #print "clock time: %f %f"%(timeTotalClock,timePerIterClock)
      #print "wall time: %f %f"%(timeTotalWall,timePerIterWall)
      break
  
  if tries >= maxTries:
    print("%s: Bailed out on particle %d after %d tries"%(thisScript, i, tries))
  
  
# ------------------------------------------------------------------------------
# Check for errors and output data to file
# ------------------------------------------------------------------------------

for i in range(numPlaced):
  if ((partData[i,1] < -tol or partData[i,1] > Lx+tol) or\
      (partData[i,2] < -tol or partData[i,2] > Ly+tol) or\
      (partData[i,3] < -tol or partData[i,3] > Lz+tol)):
    print("%s: Error: Particle %d out of bounds. Program will exit."%(thisScript, i))
    exit(1)
  for j in range(i):
    dist = sum((partData[i,1:] - partData[j,1:])**2)
    if dist - (partData[i,0] + partData[j,0] + tol)**2 < 0.0:
      print("%s: Error: Overlap detected between particles %d and %d."\
            " Program will exit."%(thisScript, i, j))
      exit(1)

partData = partData[:numPlaced,:]

# Calculate mesh domain size.
outerBounds = np.zeros((3,2))
for i in range(numPlaced):
    if partData[i,1] + partData[i,0] > outerBounds[0,1]:
        outerBounds[0,1] = partData[i,1] + partData[i,0]
    if partData[i,1] - partData[i,0] < outerBounds[0,0]:
        outerBounds[0,0] = partData[i,1] - partData[i,0]

    if partData[i,2] + partData[i,0] > outerBounds[1,1]:
        outerBounds[1,1] = partData[i,2] + partData[i,0]
    if partData[i,2] - partData[i,0] < outerBounds[1,0]:
        outerBounds[1,0] = partData[i,2] - partData[i,0]
        
    if partData[i,3] + partData[i,0] > outerBounds[2,1]:
        outerBounds[2,1] = partData[i,3] + partData[i,0]
    if partData[i,3] - partData[i,0] < outerBounds[2,0]:
        outerBounds[2,0] = partData[i,3] - partData[i,0]

# Give some space at edges of domain.
extraSpace = 0.1 * partData[0,0]
outerBounds[:,0] = outerBounds[:,0] - extraSpace
outerBounds[:,1] = outerBounds[:,1] + extraSpace

meshLx = outerBounds[0,1] - outerBounds[0,0]
meshLy = outerBounds[1,1] - outerBounds[1,0]
meshLz = outerBounds[2,1] - outerBounds[2,0]

# Calculate porosity
partVol = 0.0
for rad in partData[:,0]:
    partVol += 4.0 * math.pi * rad**3.0 / 3.0
porosity = 1.0 - partVol / (meshLx * meshLy * meshLz) 

packingFile = 'packings/packing_%d.txt'%(numPlaced)
print("\n%s: Placed %d of %d particles. Writing data to file %s."\
      %(thisScript, numPlaced, numPartsAttempt, packingFile))
outFile = open(packingFile, 'w')
outFile.write("NUMPARTS\n%d\n"%numPlaced)
outFile.write("POROSITY\n%1.4f\n"%porosity)
outFile.write("BOXDIMS\n%+1.15e %+1.15e %+1.15e\n"%(Lx, Ly, Lz))
outFile.write("PARTDATA\n")
for i in range(numPlaced):
  outFile.write("%+1.15e %+1.15e %+1.15e %+1.15e\n"\
                %(partData[i,0], partData[i,1], partData[i,2], partData[i,3]))        
outFile.close()

# Output data in .vtp format readable by ParaView. -----------------------------
packingFile = 'packings/packing_%d.vtp'%(numPlaced)
outFile = open(packingFile, 'w')
outFile.write("<?xml version=\"1.0\"?>\n<VTKFile type=\"PolyData\" version=\"0.1\" format=\"ascii\">\n")
outFile.write("<PolyData>\n\t<Piece NumberOfPoints=\"%d\">\n\t\t<Points>\n\t\t\t"%(numPlaced))
outFile.write("<DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">\n")
for i in range(numPlaced):
  outFile.write("\t\t\t\t%f %f %f\n"%(partData[i,1], partData[i,2], partData[i,3]))
outFile.write("\t\t\t</DataArray>\n\t\t</Points>\n")
outFile.write("\t\t<PointData Scalars=\"Diameter\">\n\t\t\t<DataArray type=\"Float32\" Name=\"Diameter\" format=\"ascii\">\n")
for i in range(numPlaced):
  outFile.write("\t\t\t\t%f\n"%(2*partData[i,0]))
outFile.write("\t\t\t\t</DataArray>\n")
outFile.write("\t\t\t<DataArray type=\"Int32\" Name=\"ID\" format=\"ascii\" NumberOfComponents=\"1\">\n")
for i in range(numPlaced):
  outFile.write("\t\t\t\t%d\n"%(i))
outFile.write("\t\t\t\t</DataArray>\n")
outFile.write("\t\t\t</PointData>\n\t\t</Piece>\n\t</PolyData>\n</VTKFile>")
outFile.close()

# Output data in format of MF-Unstructured AddHole command. --------------------
packingFile = 'packings/packing_%d.mf.addhole'%(numPlaced)
outFile = open(packingFile, 'w')


tetMaxVol = (meshLx * meshLy * meshLz) / 1.0e+05
tetMinAng = 15.0
tetRadEdge = 2.0
ScaleGeo = 1000.0

outFile.write("\tAutoMesh TET %+1.7e %f %f %+1.7e %+1.7e %+1.7e\n\n"\
        %(tetMaxVol * ScaleGeo**3, tetMinAng, tetRadEdge,\
        meshLx * ScaleGeo, meshLy * ScaleGeo, meshLz * ScaleGeo))

outFile.write("\tAddHoleCount %d\n\n" %(numPlaced))

# Determine refinement classes.
refineMax = 3
Amax = math.pi * partData[0,0]**2
Amin = math.pi * partData[-1,0]**2
AmaxElem = Amax / (20 * 4**(refineMax - 1))
minDelta = 1.0E30

refineMin = refineMax
for k in range(refineMax,0,-1):
    AminElem = Amin / (20 * 4**(k - 1))
    if abs(AmaxElem - AminElem) < minDelta:
        minDelta = abs(AmaxElem - AminElem)
        refineMin = k

radClasses = [partData[-1,0]]
if refineMin != refineMax:
    refineDelta = refineMax - refineMin
    bounds = 0
    for i in range(numPlaced-1, 0, -1):
        AElem0 = (math.pi * partData[i,0]**2) / (20 * 4**(refineMin + bounds - 1))
        AElem1 = (math.pi * partData[i,0]**2) / (20 * 4**(refineMin + bounds))
        if abs(AmaxElem - AElem0) > abs(AmaxElem - AElem1):
            radClasses.append(partData[i,0])
            bounds += 1
        if bounds == refineDelta:
            radClasses.append(partData[0,0])
            break

# Translate positions to place corner of cell at origin.
AddHoleData = np.empty(partData.shape)
AddHoleData[:,1] = ScaleGeo * (partData[:,1] + abs(outerBounds[0,0]))
AddHoleData[:,2] = ScaleGeo * (partData[:,2] + abs(outerBounds[1,0]))
AddHoleData[:,3] = ScaleGeo * (partData[:,3] + abs(outerBounds[2,0]))
AddHoleData[:,0] = ScaleGeo * partData[:,0]

refine = refineMax
for i in range(numPlaced):
    for j in range(len(radClasses)-1):
        if partData[i,0] > radClasses[j] and partData[i,0] < radClasses[j+1]:
            refine = refineMin + j
    outFile.write("\tAddHole %+1.7e %+1.7e %+1.7e %+1.7e 1.0 1.0 1.0 0.0 0.0 0.0 %d\n"\
                 %(AddHoleData[i,1], AddHoleData[i,2], AddHoleData[i,3], AddHoleData[i,0], refine))
outFile.write("\n\tScaleGeometryFactor %+1.7e\n" %(ScaleGeo))
outFile.close() 

# Output data in format of MF-Unstrucutred IB command
packingFile = 'packings/packing_%d.mf.ib'%(numPlaced)
outFile = open(packingFile, 'w')
outFile.write("\tAutoMesh HEX 100 100 100 %1.7e %1.7e %1.7e\n\n"%(meshLx, meshLy, meshLz))
outFile.write("IMMERSEDBOUNDARY\n\tCOUNT %d\n"%(numPlaced))
refine = refineMax
for i in range(numPlaced):
  for j in range(len(radClasses)-1):
    if partData[i,0] > radClasses[j] and partData[i,0] < radClasses[j+1]:
      refine = refineMin + j
  outFile.write("\tOBJECT %d\n\t\tTYPE MFTL\n\t\tSHAPE SPHERE\n\t\tRADIUS %1.7e\n"\
          "\t\tREFINELEVEL %d\n\t\tINITCENTER %1.7e %1.7e %1.7e\n\tENDOBJECT\n"\
          %(i, partData[i,0], refine, partData[i,1] + abs(outerBounds[0,0]),
            partData[i,2] + abs(outerBounds[1,0]), partData[i,3] + abs(outerBounds[2,0])))
outFile.write("END\n")
outFile.close()

