
# Import basic Python modules
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions
import operator                 # used to sort dictionaries
import pickle                   # save Python objects to file

import math

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob           # basic OpenBabel
ob.obErrorLog.SetOutputLevel(0)           # suppress warnings

# Import my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
from qmmm import printf
from qmmm import fprintf

inpf="/dicos_ui_home/tlankau/Structures/TingYi_Data/ABASEE_1/ABASEE-oxo_1_1_tpss_d4_lanl2dz_631gpol_Co3.opt.xyz"

def scan_inverse(at, rad):
  # td 0 to 360 deg
  # pd 0 to 180 deg
  rs = 10000.0 # final rsults sum
  rl = []      # final results vector
  for td in range(0, 360, 5):
    for pd in range(0, 180, 5):
      # pick a point on the sphere
      tx = at.GetX() + rad * math.sin(math.radians(pd)) * math.cos(math.radians(td))
      ty = at.GetY() + rad * math.sin(math.radians(pd)) * math.sin(math.radians(td))
      tz = at.GetZ() + rad * math.cos(math.radians(pd))
      sum = 0.0
      # calculate distances
      for tat in ob.OBMolAtomIter(obmol):
        if at == tat: continue
        d  = (tx - tat.GetX()) * (tx - tat.GetX())
        d += (ty - tat.GetY()) * (ty - tat.GetY())
        d += (tz - tat.GetZ()) * (tz - tat.GetZ())
        d  = math.sqrt(d)
        sum += 1.0/d
      # print(sum, rs)
      if sum < rs:
        # print(sum, rs)
        rs = sum
        rl = [tx, ty, tz]
      if rl == []:
        print("Locating the H position failed")
        exit()
  return(rl, rs)

def scan_simple_sum(at, rad):
  # td 0 to 360 deg
  # pd 0 to 180 deg
  rs = 0.0     # final rsults sum
  rl = []      # final results vector
  for td in range(0, 360, 5):
    for pd in range(0, 180, 5):
      # pick a point on the sphere
      tx = at.GetX() + rad * math.sin(math.radians(pd)) * math.cos(math.radians(td))
      ty = at.GetY() + rad * math.sin(math.radians(pd)) * math.sin(math.radians(td))
      tz = at.GetZ() + rad * math.cos(math.radians(pd))
      sum = 0.0
      # calculate distances
      for tat in ob.OBMolAtomIter(obmol):
        if at == tat: continue
        d  = (tx - tat.GetX()) * (tx - tat.GetX())
        d += (ty - tat.GetY()) * (ty - tat.GetY())
        d += (tz - tat.GetZ()) * (tz - tat.GetZ())
        d  = math.sqrt(d)
        sum += d
      # print(sum, rs)
      if sum > rs:
        # print(sum, rs)
        rs = sum
        rl = [tx, ty, tz]
      if rl == []:
        print("Locating the H position failed")
        exit()
  return(rl, rs)

printf("Read input file\n")
obConversion = ob.OBConversion()
obConversion.SetInAndOutFormats("xyz", "xyz")
obmol = ob.OBMol()
obConversion.ReadFile(obmol, inpf)
atom_number = obmol.NumAtoms()
printf("%i atoms read\n", atom_number)
printf("Search for MeO pairs\n")

for at in ob.OBMolAtomIter(obmol):
  if at.GetAtomicNum() >= 19:
    metal_idx = at.GetIdx()
    metal_symbol = ob.GetSymbol(at.GetAtomicNum())
    printf("%s found at index: %i\n", metal_symbol, metal_idx)
    break
for nat in ob.OBAtomAtomIter(at):
  ox_idx = nat.GetIdx()
  if nat.GetAtomicNum() == 8:
    printf("O atom at index %i\n", ox_idx)
    bond = at.GetBond(nat)
    bond_order = bond.GetBondOrder()
    bond_length = bond.GetLength()
    printf("The bond order between both is: %i\n", bond_order)
    printf("%sO bond length: %5.3f A\n", metal_symbol, bond_length)
    break

# scan using the inverse distance
rad = 1.5
printf("Scan O-%i with a radius of %.1f A (inverse distance)\n", ox_idx, rad)
hpinv, tval = scan_inverse(nat, rad)
print(hpinv, tval)

# scan using the inverse distance
rad = 1.5
printf("Scan O-%i with a radius of %.1f A (sum of distances)\n", ox_idx, rad)
hpsum, tval = scan_simple_sum(nat, rad)
print(hpsum, tval)

# add a new dummy atoms to mark the H atom position results
dummy_atom = obmol.NewAtom()  # this creates a new atom and adds it to the molecule
dummy_atom.SetAtomicNum(0)    # atomic number 0 indicates a dummy atom (X)
dummy_atom.SetVector(hpinv[0], hpinv[1], hpinv[2])  # add coordinates to the new atom

dummy_atom = obmol.NewAtom()  # this creates a new atom and adds it to the molecule
dummy_atom.SetAtomicNum(0)    # atomic number 0 indicates a dummy atom (X)
dummy_atom.SetVector(hpsum[0], hpsum[1], hpsum[2])  # add coordinates to the new atom

# write new xyz file
print("Write the extended molecule to disk")
obConversion.SetOutFormat("xyz")
xyz_file = "./test_sphere.xyz"
obConversion.WriteFile(obmol, xyz_file)

exit()