# Volcano_TS
- The original volcano data are not part of the repository to make tracking easier.
  /dicos_ui_home/tlankau/Structures/TingYi_Data
- 001_Position
  A collection of scripts to test methods to find the best position for the H atom of the nethane CH bond.
  + chat_gpt_example.py: The solution suggested by ChatGpt, which didn't work. Rubbish!
  + test_sphere.py: Find the point with the least amount of interactions with cluster on a sphere around the oxygen atom. Two methods were used: a) Minimizing the sum of inverted distances. b) Maximizing the sum of distances.
  + test_vec.py: The ChatGpt solution for the rotation matrix produced rubbish results. The target vector and the bond vector were misaligend by several degrees. I there replaced the bond matrix with an old fahioned one from wikipedia and changed the sign of the rotation axis by interchanging the vectoes in the cross product. The test code using a single worked well. The observed deviation from a zero result was 1D-17.
  + test_circle.py: Tesyting positions on ring around the oxygen atom. The final position is determined by the minimum of of the sum of ynverse distances.
  + test_ch4.py: Replace the H atom on the ring with an orientated CH4 molecule. Do some sanitychecks prior to merging the 2 molecules.
- 002_Select_by_Mult
  This directory contains code to get the multiplicity of the cluster lowest in energy.
  + volcano_plot_muls.csv: This file provided by Ting Yi contains the ground state multiplicities of the volcano cluster.
  + correlate.py: This code correlate the entries in volcano_plot_muls.csv with the xyz files in /dicos_ui_home/tlankau/Structures/TingYi_Data/
