"""
    Module containing all the constants of the pdb-tools package.

    Created on Thu 12-10-2023 14:27:15
    Version: 0.0.1
    Author: Arend Frederik Koens
    Title: Constants for atomic masses and acid masses
"""

from typing import Iterable

# https://www.lenntech.com/periodic-chart-elements/atomic-mass.htm
# Atomic mass in g/mol | Atomic mass in Dalton
PERIODIC_ELEMENT_MASSES = {"H": 1.0079, "He": 4.0026, "Li": 6.941, "Be": 9.0122, "B": 10.811,
                           "C": 12.0107, "N": 14.0067, "O": 15.9994, "F": 18.9984, "Ne": 20.1797,
                           "Na": 22.9897, "Mg": 24.305, "Al": 26.9815, "Si": 28.0855, "P": 30.9738,
                           "S": 32.065, "Cl": 35.453, "K": 39.0983, "Ar": 39.948, "Ca": 40.078,
                           "Sc": 44.9559, "Ti": 47.867, "V": 50.9415, "Cr": 51.9961, "Mn": 54.938,
                           "Fe": 55.845, "Ni": 58.6934, "Co": 58.9332, "Cu": 63.546, "Zn": 65.39,
                           "Ga": 69.723, "Ge": 72.64, "As": 74.9216, "Se": 78.96, "Br": 79.904,
                           "Kr": 83.8, "Rb": 85.4678, "Sr": 87.62, "Y": 88.9059, "Zr": 91.224,
                           "Nb": 92.9064, "Mo": 95.94, "Tc": 98, "Ru": 101.07, "Rh": 102.9055,
                           "Pd": 106.42, "Ag": 107.8682, "Cd": 112.411, "In": 114.818, "Sn": 118.71,
                           "Sb": 121.76, "I": 126.9045, "Te": 127.6, "Xe": 131.293, "Cs": 132.9055,
                           "Ba": 137.327, "La": 138.9055, "Ce": 140.116, "Pr": 140.9077,
                           "Nd": 144.24, "Pm": 145, "Sm": 150.36, "Eu": 151.964, "Gd": 157.25,
                           "Tb": 158.9253, "Dy": 162.5, "Ho": 164.9303, "Er": 167.259,
                           "Tm": 168.9342, "Yb": 173.04, "Lu": 174.967, "Hf": 178.49,
                           "Ta": 180.9479, "W": 183.84, "Re": 186.207, "Os": 190.23, "Ir": 192.217,
                           "Pt": 195.078, "Au": 196.9665, "Hg": 200.59, "Tl": 204.3833, "Pb": 207.2,
                           "Bi": 208.9804, "Po": 209, "At": 210, "Rn": 222, "Fr": 223, "Ra": 226,
                           "Ac": 227, "Pa": 231.0359, "Th": 232.0381, "Np": 237, "U": 238.0289,
                           "Am": 243, "Pu": 244, "Cm": 247, "Bk": 247, "Cf": 251, "Es": 252,
                           "Fm": 257, "Md": 258, "No": 259, "Rf": 261, "Lr": 262, "Db": 262,
                           "Bh": 264, "Sg": 266, "Mt": 268, "Rg": 272, "Hs": 277, "?": 0}

# https://aminoacidsguide.com/
# Molar mass in Da or g/mol
AMINO_MOL_MASS_SLC = {'A': 89.09318, 'R': 174.20096, 'N': 132.11792, 'D': 133.10268, 'C': 121.15818,
                      'E': 147.12926, 'Q': 146.1445, 'G': 75.0666, 'H': 155.15456, 'I': 131.17292,
                      'L': 131.17292, 'K': 146.18756, 'M': 149.21134, 'F': 165.18914,
                      'P': 115.13046,
                      'S': 105.09258, 'T': 119.11916, 'W': 204.22518, 'Y': 181.18854,
                      'V': 117.14634,
                      "?": 0}

AMINO_TLC_SLC = {"ALA": 'A', "ARG": 'R', "ASN": 'A', "ASP": 'D', "CYS": 'C', "GLU": 'E', "GLN": 'Q',
                 "GLY": 'G', "HIS": 'H', "ILE": 'I', "LEU": 'L', "LYS": 'K', "MET": 'M', "PHE": 'F',
                 "PRO": 'P', "SER": 'S', "THR": 'T', "TRP": 'W', "TYR": 'Y', "VAL": 'V', "UNK": '?'}

MASS_OF_WATER = PERIODIC_ELEMENT_MASSES["H"] * 2 + PERIODIC_ELEMENT_MASSES["O"]


def tlc_to_slc(seq: Iterable) -> list:
    """
    Convert a Triple Letter Code to a Single Letter Code for amino acids
    :param seq:
    :return:
    """
    return list(map(lambda x: AMINO_TLC_SLC[x], seq))


def protein_to_acid_masses(protein: Iterable[str]) -> list:
    """
    Create a generator that finds the molar mass for all the acids of a protein.
    :param protein: A string containing single lette codes for amino acids.
    :return: A list containing all the molar masses of the acids in a protein.
    """

    protein_mass_generator = [AMINO_MOL_MASS_SLC[acid] for acid in protein if
                              acid in AMINO_MOL_MASS_SLC]

    return protein_mass_generator
