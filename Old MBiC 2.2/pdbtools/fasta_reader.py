"""
    Functions for reading and manipulating Fasta files and Data.

    Created on Wed 1-11-2023 09:21:00
    Author: A F Koens
    Version: 0.0.1
    Title: Fasta Reader
"""


def get_sequences_from_fasta_file(filepath: str) -> dict[str, str]:
    """
    Extracts the sequences from a fasta file.
    :param filepath: A filepath to a fasta file.
    :return: A dict with a string of chain_id as key and str sequence as value.
    """

    sequences = {}
    with open(filepath, 'r', encoding="UTF-8") as fasta_file:
        for line in fasta_file:
            # Check for fasta header
            if line.startswith(">"):

                # Find chain ids and remove "Chain " identifier
                chain_ids_segment = line.split("|")[1][6:]

                # Chain id can have alternate ids, here remove alternate id identifiers
                chain_ids_segment = chain_ids_segment.replace("[auth ", ",").replace("]", "")
                # Split into a list of ids
                chain_ids = chain_ids_segment.split(',')
                # strip IDS
                chain_ids = map(lambda i: i.strip(), chain_ids)

                continue

            sequence = line.rstrip()
            sequences.update({chain_id: sequences[chain_id] + sequence for chain_id in chain_ids})

    return sequences
