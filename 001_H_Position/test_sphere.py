
# Import basic Python modules
import sys                      # IO Basics
import glob                     # wild cards in file names
import os                       # access file system
import math                     # mathematical  functions
import operator                 # used to sort dictionaries
import pickle                   # save Python objects to file

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob           # basic OpenBabel
ob.obErrorLog.SetOutputLevel(0)           # suppress warnings

# Import my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
from qmmm import printf
from qmmm import fprintf

inpf="/dicos_ui_home/tlankau/Structures/TingYi_Data/ABASEE_1/ABASEE-oxo_1_1_tpss_d4_lanl2dz_631gpol_Co3.opt.xyz"

def scan_sphere(at, rad):
  return()

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

rad = 1.5
printf("Scan O-%i with a radius of %f A\n", ox_idx, rad)
exit()