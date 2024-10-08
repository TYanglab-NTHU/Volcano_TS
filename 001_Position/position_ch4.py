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
# xyz_folder = "/dicos_ui_home/tlankau/Structures/TingYi_Data/ABASEE_1/"
# xyz_file   = "ABASEE-oxo_1_1_tpss_d4_lanl2dz_631gpol_Co3.opt.xyz"
xyz_folder = "/dicos_ui_home/tlankau/Structures/TingYi_Data/RUHCON_2/"
xyz_file   = "RUHCON-oxo_2_1_tpss_d4_lanl2dz_631gpol_Co3.opt.xyz"

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
def read_xyz_file(inpf):
  obConversion = ob.OBConversion()
  obConversion.SetInFormat("xyz")
  obmol = ob.OBMol()
  obConversion.ReadFile(obmol, inpf)
  return(obmol)