def get_connections(atom_number, bonds):
    # Most readable way of finding all the bond connections for a specific atom
    connections = []
    for bond in bonds:
        if atom_number in bond:
            connections.append(bond)
        else:
            continue
    return connections

def build_atom_dict(coords_and_bonds):
    # Creates useful dictionary of atom number keys and element, xyz coordinate
    # values.
    atom_dict = {}
    for atom, x, y, z in coords_and_bonds:
        atom_dict[atom[1]] = atom[0]     
    return atom_dict


def remove_added_protons(coords_and_bonds): 
    coords = coords_and_bonds[0]
    bonds = coords_and_bonds[1]   
    protons = []
    atomdict = build_atom_dict(coords_and_bonds)
   
    for atom in atomdict:
        element = atomdict.get(atom)
        number = atom
        if element == 'H':
            for bond in (get_connections(atom, bonds)):
                
                if bond[0] == atom:
                    if atomdict.get(bond[1]) == 'O':
                        protons.append(bond[0])

                elif bond[1] == atom:
                    if atomdict.get(bond[0]) == 'O':
                        protons.append(bond[1])
                else:
                    pass
        else:
            pass        
        for proton in protons:
            coords = [x for x in coords if proton not in x[0]]
            bonds  = [y for y in bonds if proton not in y]
    
    return [coords, bonds]
