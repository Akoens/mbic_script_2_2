"""
    All the dependencies for the pdb-tools package.

    Created on Wed 1-10-2023 08:40:00
    Author: A F Koens
    Version: 0.0.1

"""

# Imports modules that are part of the package
import pdbtools.pdb_reader
import pdbtools.fasta_reader
import pdbtools.pretty_f
import pdbtools.sinumber
import pdbtools.constants

# Allows for * import of package
__all__ = ["fasta_reader", "pdb_reader", "sinumber", "pretty_f", "constants"]
