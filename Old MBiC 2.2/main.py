"""
    Module containing all the main functionalities for the project.

    Created on Wed 1-11-2023 09:21:00
    Author: A F Koens
    Version: 0.0.1
    Title: Fasta Reader
"""
import textwrap

from pdbtools import pdb_reader, sinumber, fasta_reader, pretty_f
from pdbtools.constants import tlc_to_slc
from utilities import unpack_2di_to_list, calc_prominence, find_matches

LINE_LENGTH: int = 70


def __calc_prot_mass(pdb_source: str) -> None:
    # Protein Calculations
    estimated_mass = pdb_reader.estimate_protein_mass_from_atoms(pdb_source)
    estimate_seq_res_mass = pdb_reader.estimate_mass_from_seq(pdb_source)

    # Print Protein Data
    print(f"# Estimated Mass from ATOMS\n{sinumber.abbreviate_float(estimated_mass, 2)} g/Mol")
    print(f"# Mass from SEQ-RES\n{sinumber.abbreviate_float(estimate_seq_res_mass, 2)} g/Mol")
    print()


def __helix_sheet(pdb_source: str) -> None:
    # Find sequences of helices and sheets
    helices, sheets = pdb_reader.get_helix_sheet(pdb_source)

    # TODO Generalize secondary structure prominence and preference
    for secondary_structure in (helices.values(), sheets.values()):
        # Unpack all sequences to list containing all the amino acids in helices and in sheets
        all_amino_acids_in_secondary_structure = unpack_2di_to_list(secondary_structure)

        # Count the occurrences of amino acids in helices and sheets
        secondary_structure_acid_prominence = calc_prominence(all_amino_acids_in_secondary_structure)

        # Calculate the preference for acids in helices and sheets
        helix_pref, sheet_pref = pdb_reader.calc_preferences(helices_prominences, sheets_prominences)

    print("# Prominence of amino acids in helices")
    formatted_string = ", ".join([f"{acid}: {prom}%" for acid, prom in helices_prominences.items()])
    if LINE_LENGTH > 0:
        formatted_string = textwrap.fill(formatted_string, LINE_LENGTH)
    print(formatted_string, end="\n\n")

    print("# Prominence of amino acids in sheets")
    formatted_string = ", ".join([f"{acid}: {prom}%" for acid, prom in sheets_prominences.items()])
    if LINE_LENGTH > 0:
        formatted_string = textwrap.fill(formatted_string, LINE_LENGTH)
    print(formatted_string, end="\n\n")

    # Amino acid preference
    formatted_string = ", ".join([f"{acid}" for acid in helix_pref])
    print(f"# Helix Preference\n{formatted_string}")
    formatted_string = ", ".join([f"{acid}" for acid in sheet_pref])
    print(f"# Sheet Preference\n{formatted_string}\n")


def __seq_compare(pdb_source: str, fasta_source: str) -> None:
    # Get sequence from pdb and fasta for comparing
    pdb_sequences = pdb_reader.get_sequences_from_pdb_file(pdb_source)
    fasta_sequences = fasta_reader.get_sequences_from_fasta_file(fasta_source)

    print(pdb_sequences)

    print("# PDB/FASTA Sequence match")
    # TODO convert to chain_id system
    for chain_id, pdb_sequence in pdb_sequences.items():
        # Convert Pdb sequence from TLC to SLC
        pdb_sequence = tlc_to_slc(pdb_sequence)

        # Make a string from pdb sequence list
        pdb_sequence = "".join(pdb_sequence)

        # Compare and format sequences
        index_matches = find_matches(pdb_sequence, fasta_sequences[chain_id])
        match_string = pretty_f.format_match_sequence(pdb_sequence,
                                                      fasta_sequences[chain_id],
                                                      index_matches, line_length=LINE_LENGTH)
        print(f"Chain {chain_id}")
        print(match_string)


def main() -> None:
    """
    The Main function
    :return: None
    """

    pdb_source = "data/5uak.pdb"
    fasta_source = "data/rcsb_pdb_5UAK.fasta"

    __calc_prot_mass(pdb_source)
    __helix_sheet(pdb_source)
    __seq_compare(pdb_source, fasta_source)


if __name__ == '__main__':
    main()
