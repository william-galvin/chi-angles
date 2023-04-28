import json
import numpy as np
import math

from reconstruction_utils import *

class Reconstructor():
    """
    Usage:
    >>> r = Reconstructor("/gscratch/scrubbed/wgalvin/python/reconstruction.json")
    >>> angles = [90, 75, -40, 30]
    >>> nb = ...
    >>> r.reconstruct(nb, angles)
    """
    def __init__(self, path: str):
        """
        Takes path to reconstruction.json
        """
        
        with open(path, 'r') as f:
            reconstruction = json.loads(f.read())
        
        # =====Beta carbon virtualization==========
        self.CA_CB_dict = reconstruction['CA_CB_dict'] # averge length of CA-CB bond length, by AA type
        self.N_C_CA_dict = reconstruction['N_C_CA_dict'] # average N-C-CA bond angle
        self.N_C_CA_CB_dict = reconstruction['N_C_CA_CB_dict'] # avereg N-C-CA-CB dihedral angle
        
        # =====Chi Angle reconstruction===========
        self.ideal_bond_lengths = reconstruction['ideal_bond_lengths'] # average bonds lengths by AA + chi
        self.ideal_bond_angles = reconstruction['ideal_bond_angles'] # averge bond angle (NOT dihedral) by AA + chi
        
        #======misc lookup tables=================
        self.aa_symbols = reconstruction['aa_symbols']
        self.chi_atoms = reconstruction['chi_atoms']
        
        
    def reconstruct(self, neighborhood: np.ndarray, chi_angles: list[float]) -> dict:
        """
        Takes a neighborhood with standard dt as used in get_neighborhood_pipeline, 
        and a list of chi angles
        
        returns a dict of {name -> coords} for sidechain atoms that were placed, 
        including the beta carbon
        """
        
        res_id = decode_id(neighborhood['res_id'])
        AA = self.aa_symbols[res_id[0]]
        
        if AA in ['ALA', 'GLY']: return None # no Chi angles
    
        atoms = {'CA': np.zeros(3, dtype=float)}
        for atom, _res_id, coords in zip(neighborhood['atom_names'],neighborhood['res_ids'], neighborhood['coords']):
            if (decode_id(_res_id) != res_id): continue
            atom = atom.decode('utf-8').strip()
            atoms[atom] = coords
        
        #=====virtualize beta carbon============
        if 'N' not in atoms or'C' not in atoms or 'CA' not in atoms: # check that all needed atoms are found in neighborhood
            return None
    
        CB_norm = get_normal_vector(atoms['N'], atoms['C'], atoms['CA'])
        
        atoms['CB'] = get_atom_place(
             CB_norm, self.N_C_CA_CB_dict[AA], 
             atoms['C'], atoms['CA'], 
             self.CA_CB_dict[AA],
             self.N_C_CA_dict[AA]
        )
        
        placed = {'CB': atoms['CB']}
        
        
        #====place side chain atoms===========
        for chi_num in range(1, 5):
            
            if AA not in self.chi_atoms[f'chi{chi_num}']: break
            
            a1, a2, a3, a4 = self.chi_atoms[f'chi{chi_num}'][AA]
            if a2 not in atoms or a3 not in atoms or a4 not in atoms: # check that all needed atoms are found in neighborhood
                continue
                
            p1_norm = get_normal_vector(atoms[a1], atoms[a2], atoms[a3])
            chi_angle = chi_angles[chi_num - 1]
            bond_length, bond_angle = self.ideal_bond_lengths[f'{AA}{chi_num - 1}'], self.ideal_bond_angles[f'{AA}{chi_num - 1}']

            predicted_place = get_atom_place(p1_norm, chi_angle, atoms[a2], atoms[a3], bond_length, bond_angle)

            # Use predicted place downstream
            atoms[a4] = predicted_place
            
            placed[a4] = atoms[a4]
        

        return placed
    