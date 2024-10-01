# Import basic Python modules
import sys                      # IO Basics
import os                       # access file system

# Import from my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
from qmmm import printf
from qmmm import fprintf

# Import mathematical tools
import math as ma               # access math functions
import numpy as np              # numerical nathematics

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob     # basic OpenBabel
ob.obErrorLog.SetOutputLevel(0)           # suppress warnings

# define control variables
xyz_folder = "/dicos_ui_home/tlankau/Structures/TingYi_Data/ABASEE_1/"
xyz_file = "ABASEE-oxo_1_1_tpss_d4_lanl2dz_631gpol_Co3.opt.xyz"

################################################################################

def get_bond_vec(mol, at1, at2):
  # bond vector from at1 to at2
  vec = np.array([at2.GetX() - at1.GetX(),
                  at2.GetY() - at1.GetY(),
                  at2.GetZ() - at1.GetZ()])
  return(vec)

def contribution(p1, p2):
  contrib = 0.0
  d  = (p1[0]-p2[0])*(p1[0]-p2[0])
  d += (p1[1]-p2[1])*(p1[1]-p2[1])
  d += (p1[2]-p2[2])*(p1[2]-p2[2])
  d  = ma.sqrt(d)
  contrib = 1.0 / d
  return(contrib)

################################################################################

def vec_angle(v1, v2):
  u1 = v1 / np.linalg.norm(v1)
  u2 = v2 / np.linalg.norm(v2)
  angle = np.arccos(np.dot(u1, u2))
  return(angle)

def vec_axis(v1, v2):
  rv = np.cross(v1, v2)
  rv = rv / np.linalg.norm(rv)
  return(rv)

# a: normalized rotational axis, t: rotational angle
# https://en.wikipedia.org/wiki/Rotation_matrix
def build_rot_mat(a, t):
  st = ma.sin(t)
  ct = ma.cos(t)
  mt = (1.0 - ct)
  mat = np.array([[a[0]*a[0] * mt +      ct,
                   a[0]*a[1] * mt - a[2]*st,
                   a[0]*a[2] * mt + a[1]*st],
                  [a[0]*a[1] * mt + a[2]*st,
                   a[1]*a[1] * mt +      ct,
                   a[1]*a[2] * mt - a[0]*st],
                  [a[0]*a[2] * mt - a[1]*st,
                   a[1]*a[2] * mt + a[0]*st,
                   a[2]*a[2] * mt +      ct]])
  return(mat)

################################################################################

# https://en.wikipedia.org/wiki/Spherical_coordinate_system
# r: radius, t: theta angle in xy plane, p: phi angle to the zaxis
# r: 0 to inf, t: 0 to 360, p: 0 to 180
# r in Angs, t and p in radian
def sphere_to_xyz(r, t, p):
  tx = r * ma.sin(p) * ma.cos(t)
  ty = r * ma.sin(p) * ma.sin(t)
  tz = r * ma.cos(p)
  return([tx, ty, tz])

################################################################################

# read xyz file
inpf = xyz_folder + xyz_file
obConversion = ob.OBConversion()
obConversion.SetInFormat("xyz")
obmol = ob.OBMol()
obConversion.ReadFile(obmol, inpf)
atom_number = obmol.NumAtoms()
print(xyz_file)
printf("%i atoms read\n", atom_number)
print()

# locate metal atom of metal-oxygen bond
for at in ob.OBMolAtomIter(obmol):
  if at.GetAtomicNum() >= 19:
    metal = at
    break
metal_idx = metal.GetIdx()
metal_symbol = ob.GetSymbol(metal.GetAtomicNum())
printf("%s found at index: %i\n", metal_symbol, metal_idx)

# locate ozygen atom of metal-oxygen bond
for at in ob.OBAtomAtomIter(metal):
  if at.GetAtomicNum() == 8:
    oxy = at
    break
oxy_idx = oxy.GetIdx()
printf("O atom at index %i\n", oxy_idx)
print()

# get bond properties
bond = metal.GetBond(oxy)
bond_order = bond.GetBondOrder()
bond_length = bond.GetLength()
bond_vec = get_bond_vec(obmol, oxy, metal)
printf("Properties of the %s%i-O%i bond\n",
       metal_symbol, metal_idx, oxy_idx)
print ("Bond vec   :", bond_vec)
printf("Bond order : %i\n", bond_order)
printf("Bond length: %5.3f A\n", bond_length)
print()

# move O atom to the center
ox = oxy.GetX()
oy = oxy.GetY()
oz = oxy.GetZ()
for at in ob.OBMolAtomIter(obmol):
  at.SetVector(at.GetX()-ox, at.GetY()-oy, at.GetZ()-oz)

# get bond properties
bond = metal.GetBond(oxy)
bond_order = bond.GetBondOrder()
bond_length = bond.GetLength()
bond_vec = get_bond_vec(obmol, oxy, metal)
printf("Properties of the %s%i-O%i bond\n",
       metal_symbol, metal_idx, oxy_idx)
print ("Bond vec   :", bond_vec)
printf("Bond order : %i\n", bond_order)
printf("Bond length: %5.3f A\n", bond_length)
print()

# prepare for the rotation
target = np.array([0.0, 0.0, -1.0])

# Get the surroundings (rot angle & axis)
rot_angle = vec_angle(target, bond_vec)
rot_axis = vec_axis(bond_vec, target)
print("Rotation parameter")
printf("angle: %.4f rad (%.1f deg)\n", rot_angle, ma.degrees(rot_angle))
print ("axis :", rot_axis)
print()

# Build rotation matix
mat = build_rot_mat(rot_axis, rot_angle)
print("Rotation matrix")
print(mat)
print()

# Apply rotation matix
for at in ob.OBMolAtomIter(obmol):
  old_pos=np.array([at.GetX(), at.GetY(), at.GetZ()])
  new_pos = np.dot(mat, old_pos)
  at.SetVector(new_pos[0], new_pos[1], new_pos[2])

# get bond properties
bond = metal.GetBond(oxy)
bond_order = bond.GetBondOrder()
bond_length = bond.GetLength()
bond_vec = get_bond_vec(obmol, oxy, metal)
printf("Properties of the %s%i-O%i bond\n",
       metal_symbol, metal_idx, oxy_idx)
print ("Bond vec   :", bond_vec)
printf("Bond order : %i\n", bond_order)
printf("Bond length: %5.3f A\n", bond_length)
print()

# write xyz file with the result from the rotation
print("Write rotated molecule to disk")
obmol.SetTitle("molecule after translation and rotation")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_circle_01.xyz"
obConversion.WriteFile(obmol, xyz_file)
print()

# draw a circle of dummy atoms
print("Add a circle of dummy atoms to the molecule")
for theta in range(0, 360, 10):
  printf("theta: %5.1f\n", theta)
  coor = sphere_to_xyz(1.5, ma.radians(theta), ma.radians(180.0-109.47))
  # add a dummy atom to the molecule
  new_atom = obmol.NewAtom()
  new_atom.SetAtomicNum(0)
  new_atom.SetVector(coor[0], coor[1], coor[2])
print()

# write xyz file with the dummy atom
print("Write augmented molecule to disk")
obmol.SetTitle("molecule with a ring of dummy atoms")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_circle_02.xyz"
obConversion.WriteFile(obmol, xyz_file)
print()

# delete the dummy atoms for a fresh start
del_at=bool(True)
while del_at == bool(True):
  del_at=bool(False)
  for at in ob.OBMolAtomIter(obmol):
    if at.GetAtomicNum() == 0:
      obmol.DeleteAtom(at)
      del_at=bool(True)
      break

# do a circle scan
print("Search for the best point on the circle")
rs = 10000.0 # final rsults sum
rl = []      # final results vector
for theta in range(0, 360, 2):
  # pick a point on the circle
  pt = sphere_to_xyz(1.5, ma.radians(theta), ma.radians(180.0-109.47))
  # calculate sum
  sum = 0.0
  for at in ob.OBMolAtomIter(obmol):
    if at == oxy: continue
    sum += contribution([at.GetX(), at.GetY(), at.GetZ()], pt)
  if sum < rs:
    print(sum, rs)
    rs = sum
    rl = pt
if rl == []:
  print("Locating the H position failed")
print("Best point found")
printf("Sum     : %6.4f\n", rs)
print ("Position:", rl)
print()

# add the optimized point as dummy to the molecule
new_atom = obmol.NewAtom()
new_atom.SetAtomicNum(0)
new_atom.SetVector(rl[0], rl[1], rl[2])

# write xyz file with the optimized dummy atom
print("Write optimization test to disk")
obmol.SetTitle("molecule with best H atom position")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_circle_03.xyz"
obConversion.WriteFile(obmol, xyz_file)

exit()