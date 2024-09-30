
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

printf("Experiment to scan a shere for the position\n")
exit()