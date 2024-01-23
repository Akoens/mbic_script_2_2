"""
    This script is used to analyze pdb files and understand more about a given protein.

    Created on Wed 1-11-2023 09:21:00
    Author: A F Koens
    Version: 0.0.1
    Title: PDB Analyzer

    Functionalities:
    1. Measure the weight of the protein by counting the elements.
    2. Measure the weight of the protein by adding all the amino acids, plus water.
    3. Find the amino acids that are most prominent in each of the secondary structures.
    4. Compare the PDB with a FASTA
    5. Calculate the hydrophobicity of the protein and its secondary structures.

    Usage:
    main.py <pdb-file> [--histo]
    main.py <pdb-file> <fasta-file> [--histo]


"""
import logging
# Builtin.
from pathlib import Path

# Local.
from pdb_analyzer import PDBAnalyzer

# External
from docopt import docopt


def main():
    args = docopt(__doc__, version='PDB Analyzer 2.0')
    print(args)

    pdb = Path(args["<pdb-file>"])
    fasta = Path(args["<fasta-file>"])
    if not pdb.is_file() or pdb.suffix != ".pdb":
        logging.exception("File is not a .pdb")
        exit(-1)

    if not fasta.is_file() or fasta.suffix != ".fasta":
        logging.exception("File is not a .fasta")
        exit(-1)

    analyzer = PDBAnalyzer(pdb)
    analyzer.load_pdb()
    print(analyzer.analyze())

    if args["--histo"]:
        print(analyzer.compositions())

    if args["<fasta-file>"]:

        print(analyzer.compare_sequences(fasta))


if __name__ == '__main__':
    main()
