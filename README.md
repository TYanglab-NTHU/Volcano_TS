# Volcano_TS
- The original volcano data are not part of the repository to make tracking easier.
  /dicos_ui_home/tlankau/Structures/TingYi_Data
- 001_H_Position
  A collection of scripts to test methods to find the best position for the H atom of the CH bond.
  test_sphere.py
    Find the point with the least amount of interactions with cluster on a sphere around the oxygen atom. Two methods were used: a) Minimizing the sum of inverted distances. b) Maximizing the sum of distances.
  chat_gpt_example.py
    The solution suggested by ChatGpt, which didn't work. Rubbish!
  test_vec.py
    The ChatGpt solution for the rotation matrix produced rubbish results. The target vector and the bond vector were misaligend by several degrees. I there replaced the bond matrix with an old fahioned one from wikipedia and changed the sign of the rotation axis by interchanging the vectoes in the cross product. The test code using a single worked well. The observed deviation from a zero result was 1D-17.