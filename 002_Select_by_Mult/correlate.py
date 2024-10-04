# Import basic Python modules
import sys                      # IO Basics
import os                       # access file system
import glob                     # wild cards in file names

# Import from my modules
sys.path.append('/dicos_ui_home/tlankau/QMMM/Q15_CCDC')
from qmmm import printf
from qmmm import fprintf

# set main control variabla
cvs_dat = "/dicos_ui_home/tlankau/Volcano_TS/002_Select_by_Mult/volcano_plot_muls.csv"
xyz_dir = "/dicos_ui_home/tlankau/Structures/TingYi_Data"

# numbers in an xyz file name
# 1) total charge of the cluster
# 2) multiplicity for the cluster
# 3) oxidation state of the metal ion as copied from the paper

# open csv file
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

print("Ground state mulitplicities of the oxo comples")
for i in range(0, len(ccdc), 1):
  printf("%-8s", ccdc[i])
  printf("  ")
  printf("%i", int(mult[i]))
  print(" | ", end="")
  cdir = glob.glob(xyz_dir+'/'+ccdc[i]+"*")
  cdir = cdir[0]
  cdir = cdir.split('/')
  cdir = cdir[-1]
  printf("%-10s", cdir)
  charge = cdir.split("_")
  charge = int(charge[-1])
  printf("  ")
  printf("%i", charge)
  printf("  ")
  file = "%s-oxo_%i_%i" % (ccdc[i], charge, int(mult[i]))
  file = glob.glob(xyz_dir+'/'+cdir+"/"+file+"*")
  if len(file) != 1:
    print("More then 1 file found")
    exit()
  file = file[0]
  printf("%s", file)
  print()
