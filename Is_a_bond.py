import numpy
import Atomic_Number_Dictionary as AnD
import networkx as nx
from string import ascii_letters, digits
from networkx.linalg.graphmatrix import adjacency_matrix, incidence_matrix
import networkx.algorithms.isomorphism as iso
import matplotlib.pyplot as plt
from networkx.classes.function import number_of_edges, number_of_nodes, edges
from itertools import combinations

def atom_distance(x1, y1, z1, x2, y2, z2):
    return numpy.sqrt((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)

def coords_to_lists(mol_file):
    """Takes a .mol file and converts it to a python list of lists.
    The first list in the list of lists will be of all the atoms accompanied.
    By their elements and coordinates.
    The second will be bonds based on atom numbers"""
    coords = list(open(mol_file, 'r'))[4:-1]
    coord_list = []
    element_count = 0
    bond_pairs = []
    for line in coords:
        coord = []
        element = ""
        number = ""
        for char in line:
            if char in ascii_letters:
                element += char
            elif char == ' ':
                if number == "" or len(number) < 3:
                    number = ""
                    continue
                else: 
                    coord.append(numpy.float32(number))
                    number = ""
            elif char in digits or char == '-' or char =='.':
                number += char
            else:
                continue
        
        element_count += 1
        if element == '':
            atom1 = ''
            atom2 = ''
            Found = False
            for char in line:
                if not Found:
                    if char in digits:
                        atom1 +=char
                    elif char == " " and len(atom1) >=1:
                        Found = True
                elif Found:
                    if char in digits:
                        atom2 += char
                    elif char == " " and len(atom2) >=1:
                        break
                else:
                    continue
            bond_pairs.append((atom1, atom2))
        elif len(element) <= 2:
            coord.insert(0, (element, str(element_count)))
            coord_list.append(coord)
        else:
            continue

    return [coord_list, bond_pairs]


def construct_graph(coords_and_bonds):
    """Takes a 2 element list of coords and bonds and constructs
    a graph where the nodes are the atoms in the first list and 
    the edges are the bonds in the second list"""    
    coords = coords_and_bonds[0]
    bonds = coords_and_bonds[1]
    element_list = [] 

    for atom,x,y,z in coords:
        element_list.append(atom)
    
    G = nx.Graph()
    # Add atoms as nodes with attribute 'Element'
    # 'Element' will be the element of the atom
    for element in element_list:
        G.add_node(element[1], {'Element': element[0]})
    # Bonds will be edges in the graph
    # Pretty straightforward
    # Can show graph with Networkx and matplotlib
    for atom1, atom2 in bonds:    
        G.add_edge(atom1, atom2)
    return G 

def add_definite_protons(G1, Protons, Bonds):
    """ Pass a graph, a list of protons and bonds they're involved in
    to add protons to a graph missing them.
    Will likely not be an issue with .mol files but left for
    illustrative purposes"""
    for proton in Protons:
        G1.add_node(proton, {'Element': 'H'})
    for atom1, atom2 in Bonds:    
        G1.add_edge(atom1, atom2)


def get_isomorphic_subgraph_mapping(G1, G2):
    """Returns a dictionary mapping one graph onto another if one graph is a 
    subgraph of the other"""
    Gm = iso.GraphMatcher(G1,G2,node_match=iso.categorical_node_match(['Element'],['C']))
    Gm.is_isomorphic()
    return Gm.mapping
    
def get_redox_graph_mapping(G1, G2):
    """Given two graphs, returns subgraph mapping (after removing redox protons)
    ***Now also returns oxygens binded to redox protons in a tuple with the mapping***
    This is a function written for the ease of my research process"""
    if nx.number_of_nodes(G1) > nx.number_of_nodes(G2):
        G1 = G1
        x = remove_redox_protons(G1)
        G1, Hydroxy_oxygen = x[0], x[1]
        Hydroxy_oxygen.append(True)
    else:
        G2 = G2
        x = remove_redox_protons(G2)
        G2, Hydroxy_oxygen = x[0], x[1]
        Hydroxy_oxygen.append(False)
    Gm = iso.GraphMatcher(G1,G2,node_match=iso.categorical_node_match(['Element'],['C']))
    Gm.is_isomorphic()
    return (Gm.mapping,Hydroxy_oxygen)
