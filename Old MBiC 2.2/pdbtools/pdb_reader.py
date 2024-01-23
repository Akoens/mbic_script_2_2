"""
    Functions for reading and manipulating PDB files and Data.

    Created on Wed 1-11-2023 09:21:00
    Author: A F Koens
    Version: 0.0.1
    Title: PDB Reader
"""
import itertools

from pdbtools.constants import PERIODIC_ELEMENT_MASSES, MASS_OF_WATER, tlc_to_slc, \
    protein_to_acid_masses


def _get_seqres_from_line(line: str, chains: tuple[str],
                          previous_sequences: dict[str, list[str]]) -> tuple[
    tuple[str], dict[str, list[str]]]:
    # Check if line is a SEQ-RES
    if not line.startswith("SEQRES"):
        return chains, previous_sequences

    if chains is None:
        chains = ()

    if previous_sequences is None:
        previous_sequences = {}

    # Index the chain id
    chain_id = line[11:12].strip()

    # Append new chains
    chains += chain_id

    # Check if chain id has an instance in sequences
    if chain_id not in previous_sequences:
        previous_sequences[chain_id] = []

    # Get the proteins form sequence line
    amino_acids = line[19:70].rstrip().split(" ")

    # extend the chain ids protein sequence
    previous_sequences[chain_id].extend(amino_acids)

    return chains, previous_sequences


def _get_atom_mass_from_line(line: str, previous_mass: float) -> any:
    # Check if file line codes for an ATOM
    if not line.startswith("ATOM"):
        return previous_mass

    if previous_mass is None:
        previous_mass = 0.0

    # Get the Periodic Element of the ATOM
    element = line[76:78].rstrip()

    # Add the mass of the element to the total mass
    previous_mass += PERIODIC_ELEMENT_MASSES[element]

    return previous_mass


def _get_seqres_mass_from_line(line: str, previous_mass: float) -> float:
    # Check for SEQRES Lines
    if not line.startswith("SEQRES"):
        return previous_mass

    if previous_mass is None:
        previous_mass = 0.0

    # Extract Proteins
    acids = line[19:70].rstrip().split(" ")
    # Translate to Single Letter Code
    acids = tlc_to_slc(acids)
    # Count Proteins
    protein_count = len(acids)
    # Sum the mass of Proteins
    protein_mass = sum(protein_to_acid_masses(acids))

    # Calc water weight
    peptide_bounds = protein_count - 1
    water_weight = peptide_bounds * MASS_OF_WATER

    # Remove Excess water weight
    protein_mass -= water_weight

    # add previous mass
    protein_mass += previous_mass

    return protein_mass


def _get_secondary_structures_from_line(line: str, sequence: str,
                                        previous_secondary_structures: dict[
                                            str, dict[str, str]]):
    # Find helices in file
    if line.startswith("HELIX"):
        # Get helix values from line
        helix_id = line[11:14].rstrip()

        start_seq = int(line[21:25].rstrip())
        end_seq = int(line[33:37].rstrip())

        # Get the Amino acid from sequence
        previous_secondary_structures['helix'][helix_id] = sequence[start_seq - 1: end_seq]

    # Find sheets in file
    if line.startswith("SHEET"):
        # Get sheet values from line
        sheet_strand = line[7:10].rstrip()
        sheet_id = line[11:14].rstrip()

        start_seq = int(line[22:26].rstrip())
        end_seq = int(line[33:37].rstrip())

        # Get the Amino acids from sequence
        previous_secondary_structures['sheet'][sheet_strand + sheet_id] = sequence[
                                                                          start_seq - 1: end_seq]

    return previous_secondary_structures


def get_sequences_from_pdb_file(filepath: str) -> tuple[tuple[str], dict[str, list[str]]]:
    """
    Extracts the SEQRES into different chains from a pdb file.
    :param filepath: Name of pdb file to be read.
    :return: A dict containing all the TLC of the different sequence
             chains with chainID as key and sequence as value.
    """

    # Open a .pdb file and read the lines
    with open(filepath, 'r', encoding="UTF-8") as file:
        for line in file:
            chains, sequences = _get_seqres_from_line(line, chains, sequences)

    return chains, sequences


def estimate_protein_mass_from_atoms(filepath: str) -> float:
    """
    Estimate the total molar mass of a protein from a pdb file by counting the elements.
    :param filepath: A filepath to a pdb file.
    :return: The total molar mass as an int.
    """

    # Open a .pdb file and read the lines
    with open(filepath, 'r', encoding="UTF-8") as pdb_file:
        for line in pdb_file:
            total_protein_mass = _get_atom_mass_from_line(line, total_protein_mass)

    return total_protein_mass


def estimate_mass_from_seq(filepath: str) -> float:
    """
    Calculates the total molar mass of a protein from a pdb file using the SEQRES.
    :param filepath: A filepath to a pdb file.
    :return: The total molar mass as an int.
    """

    protein_mass = 0

    # Open a .pdb file and read the lines
    with open(filepath, 'r', encoding="UTF-8") as pdb_file:
        for line in pdb_file:
            protein_mass = _get_seqres_mass_from_line(line, protein_mass)
    return protein_mass


def get_secondary_structures(filepath: str) -> dict[str, dict[str, list[str]]]:
    """
    Calculate the preference of amino acids for secondary structures composition.
    :param filepath: A filepath to a pdb file.
    :return: The Preferences of amino acids in Sheets and Helices.
    """

    secondary_structures = {"helix": {}, "sheet": {}, "turn": {}}
    sequence = get_sequences_from_pdb_file(filepath)

    # Open a .pdb file and read the lines
    with open(filepath, 'r', encoding="UTF-8") as pdb_file:
        for line in pdb_file:
            _get_secondary_structures_from_line(line, sequence, secondary_structures)

    return secondary_structures


def calc_preferences(ss_w_p: dict[str: dict[str, float]]) -> dict[str, list[str]]:
    """
    Compares the prominences of all the acids in secondary structures and assigns the acids to
        the structures where they are preferred.
    :type ss_w_p: dict[str: dict[str, float]
    :param ss_w_p: Secondary Structures With Prominences is A dict width the secondary structure name
        as key and a dict of all the amino acids and their prominences as value.
    :return: a dict with the secondary structure and the list of preferred aminoacids as value.
    """

    # Initialise preferences for all the structures, structure_name as key: empty list as value.
    preferences = {structure: [] for structure in ss_w_p}

    sets_of_acids = (set(structure) for structure in ss_w_p)

    # itertools chain combines alle iterables into a single iterable, in this cast it's an
    # iterable of sets
    for acid in itertools.chain(sets_of_acids):
        for structure, prominences in ss_w_p.items():
            # os: other_structure, osp: other_structure_prominences

            # Check if acid is not in any of the other structures
            # If it is, assign the acid to itself
            if not any(acid in osp for os, osp in ss_w_p.items() if os is not structure):
                preferences[structure].append(acid)

            # If the acid is also in other structures
            # Assign the acid to this structure if its prominence is bigger than any other.
            if any(prominences[acid] > osp[acid] for os, osp in ss_w_p.values() if
                   acid in osp and os is not structure):
                preferences[structure].append(acid)

    return preferences
