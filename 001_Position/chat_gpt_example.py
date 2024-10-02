import numpy as np
import openbabel

# Helper function to calculate the cross product
def cross_product(vec1, vec2):
    return np.cross(vec1, vec2)

# Helper function to calculate the dot product
def dot_product(vec1, vec2):
    return np.dot(vec1, vec2)

# Function to compute a rotation matrix for aligning a vector with the Z axis
def rotation_matrix_to_align_with_z(bond_vector):
    z_axis = np.array([0.0, 0.0, 1.0])
    
    bond_norm = bond_vector / np.linalg.norm(bond_vector)  # Normalize the bond vector
    z_norm = z_axis / np.linalg.norm(z_axis)  # Normalize Z-axis vector
    
    # Find the angle between the bond and the Z-axis
    cos_theta = dot_product(bond_norm, z_norm)
    theta = np.arccos(cos_theta)
    
    # Find the axis of rotation (cross product of bond vector and Z-axis)
    axis_of_rotation = cross_product(bond_norm, z_norm)
    
    if np.linalg.norm(axis_of_rotation) != 0:
        axis_of_rotation = axis_of_rotation / np.linalg.norm(axis_of_rotation)  # Normalize axis of rotation
    
    # Create a rotation matrix around the axis of rotation by theta angle
    a = np.cos(theta / 2)
    b, c, d = -axis_of_rotation * np.sin(theta / 2)
    
    rotation_matrix = np.array([[a*a + b*b - c*c - d*d, 2*(b*c - a*d), 2*(b*d + a*c)],
                                [2*(b*c + a*d), a*a + c*c - b*b - d*d, 2*(c*d - a*b)],
                                [2*(b*d - a*c), 2*(c*d + a*b), a*a + d*d - b*b - c*c]])
    
    return rotation_matrix

# Function to apply a rotation matrix to all atoms in the molecule
def apply_rotation(mol, rotation_matrix):
    for atom in openbabel.OBMolAtomIter(mol):
        coord = np.array([atom.GetX(), atom.GetY(), atom.GetZ()])
        new_coord = np.dot(rotation_matrix, coord)
        atom.SetVector(new_coord[0], new_coord[1], new_coord[2])

# Open a molecule using OpenBabel
obConversion = openbabel.OBConversion()
obConversion.SetInAndOutFormats("xyz", "xyz")

mol = openbabel.OBMol()
obConversion.ReadFile(mol, "input.xyz")  # Read your molecule file here

# Select two atoms forming the bond to align with the Z-axis
atom1 = mol.GetAtom(1)  # Replace with correct atom index
atom2 = mol.GetAtom(2)  # Replace with correct atom index

# Compute the bond vector
bond_vector = np.array([atom2.GetX() - atom1.GetX(),
                        atom2.GetY() - atom1.GetY(),
                        atom2.GetZ() - atom1.GetZ()])

# Compute the rotation matrix to align the bond with the Z-axis
rot_matrix = rotation_matrix_to_align_with_z(bond_vector)

# Apply the rotation to the molecule
apply_rotation(mol, rot_matrix)

# Write the rotated molecule to a new file
obConversion.WriteFile(mol, "rotated_output_z_axis.xyz")
