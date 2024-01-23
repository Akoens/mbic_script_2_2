This script is used to analyze pdb files and understand more about a given protein.

Created on Wed 1-11-2023 09:21:00  
Author: A F Koens

How to use:
> pip install -r requirements.txt  
> python main.py "pdb_file_path" "fasta_file_path"


Functionalities:
1. Measure the weight of the protein by counting the elements.
2. Measure the weight of the protein by adding all the amino acids, plus water.
3. Find the amino acids that are most prominent in each of the secondary structures.
4. Compare the PDB with a FASTA
5. Calculate the hydrophobicity of the protein and its secondary structures.
