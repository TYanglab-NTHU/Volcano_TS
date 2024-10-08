# Import basic Python modules
import sys                      # IO Basics
import os                       # access file system
import glob                     # wild cards in file names

# Import mathematical tools
import math as ma               # access math functions
import numpy as np              # numerical nathematics

# Import OpenBabel to assist with atom classification
from openbabel import openbabel as ob     # basic OpenBabel
ob.obErrorLog.SetOutputLevel(0)           # suppress warnings

# Import from my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
from qmmm import printf
from qmmm import fprintf

# Import from my code to add the CH4 molecule
sys.path.append('/dicos_ui_home/tlankau/Volcano_TS/001_Position')
from position_ch4 import get_bond_vec
from position_ch4 import contribution
from position_ch4 import del_dummies
from position_ch4 import vec_angle
from position_ch4 import vec_axis
from position_ch4 import build_rot_mat
from position_ch4 import sphere_to_xyz

from position_ch4 import read_xyz_file

# set main control variabla
xyz_dir = "/dicos_ui_home/tlankau/Structures/TingYi_Data"
cvs_dat = xyz_dir + "/volcano_plot_muls.csv"

# numbers in an xyz file name
# 1) total charge of the cluster
# 2) multiplicity for the cluster
# 3) oxidation state of the metal ion as copied from the paper

# importatnt lists to handle the files
ccdc   = []   # CCDC name tags
mult   = []   # multiplicity of the cluster
charge = []   # total charge of the cluster
xyz    = []   # absolute path to the xyz files 

# open csv file
print("Sanity check of the csv fil:")
with open(cvs_dat) as file:
  for line in file:
    line = line.rstrip()
    line = line.split(',')
    col1 = line.pop(0)
    if col1 == "":
      ccdc = line
    if col1 == "oxo_mul":
      mult = line
printf("ccdc columns: %i\n", len(ccdc))
printf("mult columns: %i\n", len(mult))
if (len(ccdc) != len(mult)):
  print("Column numbers don't match")
  exit()
print()

# correlate the csv file with the xyz files
for i in range(0, len(ccdc), 1):
  mult[i] = int(mult[i]) # type cast multiplicity
  cdir = glob.glob(xyz_dir+'/'+ccdc[i]+"*")
  cdir = cdir[0]
  cdir = cdir.split('/')
  cdir = cdir[-1]
  cha = cdir.split("_")
  cha = int(cha[-1])
  charge.append(cha)
  file = "%s-oxo_%i_%i" % (ccdc[i], charge[i], int(mult[i]))
  file = glob.glob(xyz_dir+'/'+cdir+"/"+file+"*")
  if len(file) != 1:
    print("More then 1 file found")
    exit()
  file = file[0]
  xyz.append(file)

print("CCDC charge mult file")
for i in range(0, len(ccdc), 1):
  printf("%-8s  ", ccdc[i])
  printf("%i  ", charge[i])
  printf("%i  ", mult[i])
  printf("%s\n", xyz[i])
print()

# work through the file list
for i in range(0, len(ccdc), 1):
  printf("File   %s\n", xyz[i])
  printf("Charge %i\n", charge[i])
  printf("Mult   %i\n", mult[i])
  
  # read cluster molecule
  obmol_clu = read_xyz_file(xyz[i])
  atom_number = obmol_clu.NumAtoms()
  printf("%i atoms read for %s\n", atom_number, ccdc[i])

  # locate metal atom of metal-oxygen bond
  for at in ob.OBMolAtomIter(obmol_clu):
    if at.GetAtomicNum() >= 19 and at.IsMetal():
      metal = at
      break
  metal_idx = metal.GetIdx()
  metal_symbol = ob.GetSymbol(metal.GetAtomicNum())
  printf("%s found at index: %i\n", metal_symbol, metal_idx)

  # locate ozygen atom of metal-oxygen bond
  for at in ob.OBAtomAtomIter(metal):
    if at.GetAtomicNum() == 8:
      valence = at.GetTotalValence()
      if valence == 1:
        # print(valence)
        oxy = at
        break
  oxy_idx = oxy.GetIdx()
  printf("O atom (valence %i) at index %i\n", valence, oxy_idx)
  print()

  # move O atom to the center
  ox = oxy.GetX()
  oy = oxy.GetY()
  oz = oxy.GetZ()
  for at in ob.OBMolAtomIter(obmol_clu):
    at.SetVector(at.GetX()-ox, at.GetY()-oy, at.GetZ()-oz)

  # build the rotation environment
  bond_vec = get_bond_vec(oxy, metal)
  target = np.array([0.0, 0.0, -1.0])
  print("Vectors")
  print ("Bond   vec   :", bond_vec)
  print ("Target vec   :", target)
  print()

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

  # read cluster molecule
  xyz_met = "/dicos_ui_home/tlankau/Volcano_TS/003_Start_Geometry_Opt/methane.xyz"
  obmol_met = read_xyz_file(xyz_met)
  atom_number = obmol_met.NumAtoms()
  printf("%i atoms read for methane\n", atom_number)

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
  # print("oh_vec", oh_vec, "(target)")
  hc_vec = get_bond_vec(hmet, cmet)
  # print("hc_vec", hc_vec, "(bond)")
  # print()

  # Get the surroundings (rot angle & axis)
  rot_angle = vec_angle(oh_vec, hc_vec)
  rot_axis = vec_axis(hc_vec, oh_vec)
  # print("Rotation parameter")
  # printf("angle: %.4f rad (%.1f deg)\n", rot_angle, ma.degrees(rot_angle))
  # print ("axis :", rot_axis)
  # print()

  # Build rotation matix
  mat = build_rot_mat(rot_axis, rot_angle)
  # print("Rotation matrix")
  # print(mat)
  # print()

  # Apply rotation matix
  for at in ob.OBMolAtomIter(obmol_met):
    old_pos=np.array([at.GetX(), at.GetY(), at.GetZ()])
    new_pos = np.dot(mat, old_pos)
    at.SetVector(new_pos[0], new_pos[1], new_pos[2])

  # move the H atom of the methane molecule to the position on the ring
  for at in ob.OBMolAtomIter(obmol_met):
    at.SetVector(at.GetX()+oh_vec[0], at.GetY()+oh_vec[1], at.GetZ()+oh_vec[2])

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
  titel_str = (f'%s %s oxo CH4, charge %i mult %i bond %i %i') % (ccdc[i], ob.GetSymbol(metal.GetAtomicNum()), charge[i], mult[i], oxy.GetIdx()-1, obmol_clu.NumAtoms()-4)
  obmol_clu.SetTitle(titel_str)
  obConversion = ob.OBConversion()
  obConversion.SetOutFormat("xyz")
  xyz_file = "./Result_XYZ_Files/"+ccdc[i].lower()
  xyz_file += (f'-oxo-met_%i_%i.xyz') % (charge[i], mult[i])
  printf("Write %s\n", xyz_file)
  obConversion.WriteFile(obmol_clu, xyz_file)
  print()
  print()

  # clean up
  del obmol_clu, obmol_met

exit()