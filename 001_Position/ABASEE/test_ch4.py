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

def get_bond_vec(at1, at2):
  # bond vector from at1 to at2
  vec = np.array([at2.GetX() - at1.GetX(),
                  at2.GetY() - at1.GetY(),
                  at2.GetZ() - at1.GetZ()])
  return(vec)

# function to estimate the interacrion between two atoms for finding the best
# position of the H atom
def contribution(p1, p2):
  contrib = 0.0
  d  = (p1[0]-p2[0])*(p1[0]-p2[0])
  d += (p1[1]-p2[1])*(p1[1]-p2[1])
  d += (p1[2]-p2[2])*(p1[2]-p2[2])
  d  = ma.sqrt(d)
  contrib = 1.0 / d
  return(contrib)

# delete all dummy atoms for a fresh start
def del_dummies(obmol):
  del_at=bool(True)
  while del_at == bool(True):
    del_at=bool(False)
    for at in ob.OBMolAtomIter(obmol):
      if at.GetAtomicNum() == 0:
        obmol_clu.DeleteAtom(at)
        del_at=bool(True)
        break
  return

################################################################################

def vec_angle(v1, v2):
  u1 = v1 / np.linalg.norm(v1)
  u2 = v2 / np.linalg.norm(v2)
  dot_prod = np.dot(u1, u2)
  if dot_prod >= 1.0:
    return(0.0)
  if dot_prod <= -1.0:
    return(ma.pi)
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
obmol_clu = ob.OBMol()
obConversion.ReadFile(obmol_clu, inpf)
atom_number = obmol_clu.NumAtoms()
print(xyz_file)
printf("%i atoms read\n", atom_number)
print()

# locate metal atom of metal-oxygen bond
for at in ob.OBMolAtomIter(obmol_clu):
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
bond_vec = get_bond_vec(oxy, metal)
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
for at in ob.OBMolAtomIter(obmol_clu):
  at.SetVector(at.GetX()-ox, at.GetY()-oy, at.GetZ()-oz)

# get bond properties
bond = metal.GetBond(oxy)
bond_order = bond.GetBondOrder()
bond_length = bond.GetLength()
bond_vec = get_bond_vec(oxy, metal)
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
for at in ob.OBMolAtomIter(obmol_clu):
  old_pos=np.array([at.GetX(), at.GetY(), at.GetZ()])
  new_pos = np.dot(mat, old_pos)
  at.SetVector(new_pos[0], new_pos[1], new_pos[2])

# get bond properties
bond = metal.GetBond(oxy)
bond_order = bond.GetBondOrder()
bond_length = bond.GetLength()
bond_vec = get_bond_vec(oxy, metal)
printf("Properties of the %s%i-O%i bond\n",
       metal_symbol, metal_idx, oxy_idx)
print ("Bond vec   :", bond_vec)
printf("Bond order : %i\n", bond_order)
printf("Bond length: %5.3f A\n", bond_length)
print()

# write xyz file with the result from the rotation
print("Write rotated molecule to disk (test_ch4_01.xyz)")
obmol_clu.SetTitle("molecule after translation and rotation")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_ch4_01.xyz"
obConversion.WriteFile(obmol_clu, xyz_file)
print()

# do a circle scan
print("Search for the best point on the circle")
rs = 10000.0 # final rsults sum
rl = []      # final results vector
for theta in range(0, 360, 2):
  # pick a point on the circle
  pt = sphere_to_xyz(1.5, ma.radians(theta), ma.radians(180.0-109.47))
  # calculate sum
  sum = 0.0
  for at in ob.OBMolAtomIter(obmol_clu):
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
oh_vec = np.array([rl[0], rl[1], rl[2]])
print ("OH_vec  :", oh_vec)
print()

# add the optimized point as dummy to the molecule
new_atom = obmol_clu.NewAtom()
new_atom.SetAtomicNum(0)
new_atom.SetVector(rl[0], rl[1], rl[2])

# write xyz file with the optimized dummy atom
print("Write optimization test to disk (test_ch4_02.xyz)")
obmol_clu.SetTitle("molecule with best H atom position (test_ch4_02.xyz)")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_ch4_02.xyz"
obConversion.WriteFile(obmol_clu, xyz_file)
print()

# delete all dummy atoms prior to the addition of methane
del_dummies(obmol_clu)

################################################################################

# read methane.xyz in a new molecule
obConversion = ob.OBConversion()
obConversion.SetInFormat("xyz")
obmol_met = ob.OBMol()
obConversion.ReadFile(obmol_met, "./methane.xyz")
atom_number = obmol_met.NumAtoms()
printf("%i atoms read from file\n", atom_number)
print("Formula:", obmol_met.GetFormula())
print("Molecular Weight:", obmol_met.GetMolWt())
print()

# move the H atom into the center
hmet = obmol_met.GetAtom(2)
hx = hmet.GetX()
hy = hmet.GetY()
hz = hmet.GetZ()
for at in ob.OBMolAtomIter(obmol_met):
  at.SetVector(at.GetX()-hx, at.GetY()-hy, at.GetZ()-hz)

# prepare for the rotation of the methane molecule
cmet = obmol_met.GetAtom(1)
hmet = obmol_met.GetAtom(2)
print("oh_vec", oh_vec, "(target)")
hc_vec = get_bond_vec(hmet, cmet)
print("hc_vec", hc_vec, "(bond)")
print()

# Get the surroundings (rot angle & axis)
rot_angle = vec_angle(oh_vec, hc_vec)
rot_axis = vec_axis(hc_vec, oh_vec)
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
for at in ob.OBMolAtomIter(obmol_met):
  old_pos=np.array([at.GetX(), at.GetY(), at.GetZ()])
  new_pos = np.dot(mat, old_pos)
  at.SetVector(new_pos[0], new_pos[1], new_pos[2])

# move the H atom of the methane molecule to the position on the ring
for at in ob.OBMolAtomIter(obmol_met):
  at.SetVector(at.GetX()+oh_vec[0], at.GetY()+oh_vec[1], at.GetZ()+oh_vec[2])

# write xyz file of the methane mokecule to disk
print("Write reorientated methane molecule to disk (test_ch4_03.xyz)")
obmol_met.SetTitle("CH4 generated from smiles 'C'")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_ch4_03.xyz"
obConversion.WriteFile(obmol_met, xyz_file)
print()

# do a sanity test before adding the CH4 coordinates to the cluster
# metal metal atom
# oxy   oxo ligand atom
# hmet  CH4 H atom
# cmet  CH4 C atom
oh_v = get_bond_vec(oxy, hmet)
oc_v = get_bond_vec(oxy, cmet)
ang  = ma.degrees(vec_angle(oh_v, oc_v))
print("Quick sanity test")
print("OH vec", oh_v, f"{np.linalg.norm(oh_v):5.3f}", "A")
print("OC vec", oc_v, f"{np.linalg.norm(oc_v):5.3f}", "A")
print("angle ", ang, "deg")
if abs(ang) > 1.0:
  print("The OH and OC vectots are misaligned")
  exit()
if np.linalg.norm(oh_v) >= np.linalg.norm(oc_v):
  print("The OH bond would be longer than the OC bond")
  exit()
print("Sanity check passed -> merge molecules")
print()

# add the methane molecule to the cluster
for at in ob.OBMolAtomIter(obmol_met):
  new_atom = obmol_clu.NewAtom()
  new_atom.SetAtomicNum(at.GetAtomicNum())
  new_atom.SetVector(at.GetX(), at.GetY(), at.GetZ())

# write xyz file of the augmented cluster to disk
print("Write extended cluster to disk (test_ch4_04.xyz)")
obmol_clu.SetTitle("oxo-cluster with CH4 molecule")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_ch4_04.xyz"
obConversion.WriteFile(obmol_clu, xyz_file)
print()

exit()