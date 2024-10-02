# Import basic Python modules
import sys                      # IO Basics
import os                       # access file system

# Import mathematical tools
import math as ma               # access math functions
import numpy as np              # numerical nathematics

# Import my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
from qmmm import printf
from qmmm import fprintf

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

# Define both vectors
target = np.array([0.0, 0.0, 1.0]) # z axis
bond   = np.array([1.0, 0.0, 0.0]) # x axis
print("Start geometries")
print("target:", target)
print("bond  :", bond)
print()

# Normalize both vectors
target = target / np.linalg.norm(target)
bond = bond / np.linalg.norm(bond)
print("Normalized vectors")
print("target:", target)
print("bond  :", bond)
print()

# Get the surroundings (rot angle & axis)
rot_angle = vec_angle(target, bond)
rot_axis = vec_axis(bond, target)
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
new_bond = np.dot(mat, bond)
print("After rotation")
print("target    :", target)
print("bond      :", bond)
print("new bond  :", new_bond)
exit()