import numpy
import os
import Atomic_Number_Dictionary as AnD
import networkx as nx
import networkx.algorithms.isomorphism as iso
import Is_a_bond as Iab
import Find_Protons 
try:
    import matplotlib.pyplot as plt
except:
    print('matplotlib not found. Will not be able to print graphs')

csv_location = "/A_PATH/"

def get_equivalent_atoms(mol_file1, mol_file2):
    """Get a dictionary of equivalent atom pairs based on bonds
    
    mol_file1 -- path to a .mol file
    mol_file2 -- path to a second .mol file
    
    Will return an empty dictionary if no mapping is found """
    
    # Get coordinates and number atoms from a mol file
    coords_and_bonds1 = Iab.coords_to_lists(mol_file1)

    # Construct graph where nodes are atom numbers and edges are bonds
    G1 = Iab.construct_graph(coords_and_bonds1)
    
    # Get coordinates and number atoms for molecule 2
    coords_and_bonds2 = Iab.coords_to_lists(mol_file2)
    
    # Construct graph where nodes are atom numbers and edges are bonds 
    # for molecule 2
    G2 = Iab.construct_graph(coords_and_bonds2)
    
    # Find the combination of atoms that exclude the hydroxyl protons
    # and return the mapping dictionary
    
    return [Iab.get_redox_graph_mapping(G1, G2), coords_and_bonds1, \
            coords_and_bonds2, G1, G2]
        
    
def get_max_delta_and_bonds(mapping_and_coords, name):
    """Return top three bond length changes and their atom number bonds
    
    **Updated to also return highest non carbonyl to hydroxyl oxygen bond 
    length change**"""
    #mapping_dictionary
    md = mapping_and_coords[0][0]
    #First molecule coordinates
    coords0 = mapping_and_coords[1][0]
    #Second molecule coordinates
    coords1 = mapping_and_coords[2][0]
    #First molecule bonds
    #This is the only bond set needed because they should be matched up
    bonds = mapping_and_coords[1][1]
    coord_dict0, coord_dict1 = {}, {}
    G1, G2 = mapping_and_coords[3], mapping_and_coords[4]
    
    for atom, x, y, z in coords0:
        coord_dict0[atom[1]] = (x, y, z)
    
    for atom, x, y, z in coords1:
        coord_dict1[atom[1]] = (x, y, z)
            
    atom_dict = Find_Protons.build_atom_dict(coords0)

    atoms_used_list = []
    bond_lengths_list = []
    
    for bond in nx.edges(G1):
        atom_1 = coord_dict0.get(bond[0])
        atom_2 = coord_dict0.get(bond[1])
        atom_3 = coord_dict1.get(md.get(bond[0]))
        atom_4 = coord_dict1.get(md.get(bond[1]))

        print(bond)
        print(atom_1)
        print(atom_2)
        print('.... .... .... ....')
        print(atom_3)
        print(atom_4)
        
        # Length change calculation outlined
        length0 = Iab.atom_distance(*(atom_1 + atom_2))
        length1 = Iab.atom_distance(*(atom_3 + atom_4))
        print((abs(length0-length1)))
        #Produce a list of bond lengths and bond identites
        bond_lengths_list.append((abs(length0-length1), bond))
        bond_lengths_list = sorted(bond_lengths_list, key = lambda tup: tup[0], reverse=True)

    #Variables assigned for readability
    third_to_max  = []
    bond3         = []
    second_to_max = []
    bond2         = []
    max_delta_bond_length = []
    bond1         = []

    
    #Try//except statements here for creating a list of failing files
    #Exclude this for normal use
    fourth_to_max         = bond_lengths_list[3][0]
    bond4                 = bond_lengths_list[3][1]
    third_to_max          = bond_lengths_list[2][0]
    bond3                 = bond_lengths_list[2][1]
    second_to_max         = bond_lengths_list[1][0]
    bond2                 = bond_lengths_list[1][1]
    max_delta_bond_length = bond_lengths_list[0][0]
    bond1                 = bond_lengths_list[0][1]

    #This is where I find atom pairs my program is unable to match
    #Get the highest bond not directly involved with the reduction mechanism
    highest_non = None
    Redox_involved_oxygen = mapping_and_coords[0][1]
    print(Redox_involved_oxygen)
    for bond in [bond1, bond2, bond3]:
        if Redox_involved_oxygen[-1] == True:
            a = bond[0] not in Redox_involved_oxygen
            b = bond[1] not in Redox_involved_oxygen
        else:
            a = md.get(bond[0]) not in Redox_involved_oxygen
            b = md.get(bond[1]) not in Redox_involved_oxygen
            
        if a and b:
            if bond == bond1:
                non_bond = bond1
                highest_non = max_delta_bond_length
            elif bond == bond2:
                non_bond = bond2
                highest_non = second_to_max
            elif bond == bond3:
                non_bond = bond3
                highest_non = third_to_max
            break
        else:
            continue
        
    print(bond_lengths_list)
    print(highest_non)

    template_string = "{}, {}-{},"

    first = template_string.format(max_delta_bond_length, atom_dict.get(bond1[0]), atom_dict.get(bond1[1]))
    second = template_string.format(second_to_max, atom_dict.get(bond2[0]), atom_dict.get(bond2[1]))
    third = template_string.format(third_to_max, atom_dict.get(bond3[0]), atom_dict.get(bond3[1]))
    fourth = template_string.format(fourth_to_max, atom_dict.get(bond4[0]), atom_dict.get(bond4[1]))
    fifth = template_string.format(highest_non, atom_dict.get(non_bond[0]), atom_dict.get(non_bond[1]))
    return first + second + third + fourth + fifth
    
