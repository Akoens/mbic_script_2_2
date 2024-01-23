"""
    The module containing the Heart of operation.
"""

# Builtin.
from itertools import repeat
from pathlib import Path

# Local.
from utilities import find_matches, split_str, calc_prominence

# External.
from Bio.Data.IUPACData import protein_letters_1to3, protein_letters_3to1, \
    protein_weights
from Bio.SeqUtils.ProtParamData import kd
from periodictable import elements
from humanize import scientific


class PDBAnalyzer:
    def __init__(self, pdb_file: Path):
        """
        Initializing pdb variables
        :param pdb_file: Path to a pdb file.
        """
        self.pdb_file_path: Path = pdb_file
        self.chain_order = []
        self.chain_sequence = {}
        self.sequence = ""
        self.helix = {}
        self.sheet = {}
        self.turn = {}
        self.protein_mass = 0.0
        self.protein_atom_mass = 0.0

        # String formatting.
        self.min_line_length = 10
        self.max_line_length = 70

    def _scan_helix_line(self, line: str):
        """
        Reads and stores the desired values, id and sequence from a Helix record.
        :param line: A line with record name HELIX
        """
        # Get helix identifier from line.
        helix_id = line[11:14].strip().rstrip()

        # Get helix position in sequence.
        start_seq = int(line[21:25].strip().rstrip())
        end_seq = int(line[33:37].strip().rstrip())

        # Get a coding section form sequence for the helix.
        self.helix[helix_id] = self.sequence[start_seq - 1: end_seq]

    def _scan_sheet_line(self, line: str):
        """
        Reads and stores the desired values, id, strand and sequence from a Sheet record.
        :param line: A line with record name SHEET
        """
        # Get identifier values from line.
        sheet_strand = line[7:10].strip().rstrip()
        sheet_id = line[11:14].strip().rstrip()

        # Get sheet position in sequence.
        start_seq = int(line[22:26].strip().rstrip())
        end_seq = int(line[33:37].strip().rstrip())

        # Get a coding section form sequence for the sheet.
        self.sheet[f"{sheet_id}-{sheet_strand}"] = self.sequence[start_seq - 1: end_seq]

    def _scan_turn_line(self, line: str):
        """
        Reads and stores the desired values, id and sequence from a Turn record.
        :param line: A line with record name TURN
        """
        # Get turn identifier from line.
        turn_id = line[11:14].strip().rstrip()

        # Get turn position in sequence.
        start_seq = int(line[22:26].strip().rstrip())
        end_seq = int(line[33:37].strip().rstrip())

        # Get a coding section form sequence for the turn.
        self.sheet[turn_id] = self.sequence[start_seq - 1: end_seq]

    def _scan_sequence_line(self, line):
        """
        Reads and stores the sequence of a protein. One as chain snippets and one as the complete
        sequence.
        :param line: A line with record name HELIX
        """
        # Get value form line.
        chain_id = line[11:12].strip().rstrip()
        sequence_chunk = line[19:70].strip().rstrip().split(" ")

        # Check if amino is in protein letters,
        # Convert to single letter code amino.
        sequence_chunk = [
            protein_letters_3to1[
                amino.capitalize()] if amino.capitalize() in protein_letters_3to1 else "?" for
            amino in sequence_chunk]

        # Convert to string.
        sequence_chunk = "".join(sequence_chunk)

        # Check if amino is in protein weights,
        # Sum the weight of all the amino acids.
        sequence_amino_mass = sum(
            [protein_weights[amino] for amino in sequence_chunk if amino in protein_weights])

        # Get the mass of water and remove the mass for every peptide bond.
        mass_water = (elements.symbol("H").mass * 2 + elements.symbol("O").mass)
        sequence_amino_mass -= (len(sequence_chunk) - 1) * mass_water

        # Initialize chain sequence if chain is not present.
        if chain_id not in self.chain_sequence:
            self.chain_sequence[chain_id] = ""

        # Update values.
        if chain_id not in self.chain_order:
            self.chain_order += chain_id
        self.chain_sequence[chain_id] += sequence_chunk
        self.sequence += sequence_chunk
        self.protein_mass += sequence_amino_mass

    def _scan_atom_line(self, line: str):
        """
        Reads and stores the desired values, Atom symbol from an Atom record.
        :param line: A line with record name ATOM
        """
        # Get the Periodic Element of the ATOM.
        element = line[76:78]
        element = element.strip().rstrip()

        # Add the mass of the element to the total mass.
        self.protein_atom_mass += elements.symbol(element).mass

    def _calc_preference(self) -> list[list[tuple[str, float]]]:
        """
        Loops over all secondary structures and finds for every amino acid it's preferred structure.
        :return: A list with lists for every structure containing tuples with the amino acids, and
         its preference percentage
        """
        # Glue all the acids of the structures together.
        structures = (self.helix, self.sheet, self.turn)
        structures = ["".join(structure.values()) for structure in structures]

        structures_prom = [[] for _ in structures]

        for acid1, acid3 in protein_letters_1to3.items():
            highest_prom_structure_index = 0
            highest_prom = 0.0
            for index, structure in enumerate(structures):

                # Check if acid exists in structure else ignore.
                if acid1 not in structure:
                    continue

                acid_prom = calc_prominence(structure, acid1)

                # Compare the value with other structures.
                if acid_prom > highest_prom:
                    highest_prom_structure_index = index
                    highest_prom = acid_prom
            structures_prom[highest_prom_structure_index].append((acid3, highest_prom))

        return structures_prom

    @staticmethod
    def _histogram(percentages: dict[str, float]) -> str:
        """
        Formats a dict of percentages to a histogram.
        :param percentages: Dict with name as string and its percentage as float.
        :return: A string format of the histogram.
        """
        formatted_histogram = ""
        for element, percentage in percentages.items():
            bar = '#' * round(percentage * 100)
            formatted_histogram += f"{element}: {bar}\n"
        return formatted_histogram

    @staticmethod
    def _calc_composition(sequence: str):
        """
            Calculates the amino acid composition of a sequence.
            :param sequence: String of single letter amino acids.
            :return: Percentage of composition per amino acid.
        """

        composition = {}
        total = len(sequence)

        if total <= 0:
            return composition

        for acid1, acid3 in protein_letters_1to3.items():
            composition[acid3] = sequence.count(acid1) / total
        return composition

    @staticmethod
    def _calc_hydrophobicity(sequence: str) -> float:
        """
        Calculates the total hydrophobicity of a sequence.
        :param sequence: String of single letter amino acids.
        :return: Score for hydrophobicity
        """
        return sum(kd[acid] for acid in sequence)

    @staticmethod
    def _load_fasta_sequence(fasta_file_path: Path) -> str:
        """
        Load a sequence from a fasta file.
        :param fasta_file_path: A Path to a fasta file.
        :return: A string of one letter amino acids.
        """
        sequence = ""
        with fasta_file_path.open(mode="r", encoding="UTF-8") as fasta_file:
            for line in fasta_file:
                if line.startswith(">"):
                    continue
                sequence += line.strip().rstrip()
        return sequence

    def clear_data(self) -> None:
        """
        Reset PDB variables.
        :return: None.
        """
        self.chain_order = []
        self.chain_sequence = {}
        self.sequence = ""
        self.helix = {}
        self.sheet = {}
        self.turn = {}
        self.protein_mass = 0.0
        self.protein_atom_mass = 0.0

    def load_pdb(self) -> None:
        """
        Load all the desired variables from the PDB records.
        :return: None.
        """
        # Clear previous data before loading new data,
        # Prevent conflicts in data.
        self.clear_data()

        # Read the PDB file line by line and handle records accordingly.
        with self.pdb_file_path.open(mode="r", encoding="UTF-8") as pdb_file:
            for line in pdb_file:
                if line.startswith("ATOM"):
                    self._scan_atom_line(line)
                elif line.startswith("SEQRES"):
                    self._scan_sequence_line(line)
                elif line.startswith("HELIX"):
                    self._scan_helix_line(line)
                elif line.startswith("SHEET"):
                    self._scan_sheet_line(line)
                elif line.startswith("TURN"):
                    self._scan_turn_line(line)

    def get_protein_mass(self) -> str:
        """
        Get the Protein mass from its sequence in grams per mol as a string.
        :return: A formatted string with protein mass.
        """
        return f"Protein mass: {scientific(self.protein_mass)} g/mol\n"

    def get_atom_mass(self) -> str:
        """
            Get the Protein mass from its atoms in grams per mol as a string.
            :return: A formatted string with protein atom mass.
        """
        return f"Atom mass: {scientific(self.protein_atom_mass)} g/mol\n"

    def hydrophobicity(self) -> str:
        """
        Get all the hydrophobicity scores, for the entire protein and the secondary structures.
        :return: A formatted string with all the hydrophobicity scores.
        """
        total_hydro = self._calc_hydrophobicity(self.sequence)
        helix_hydro = self._calc_hydrophobicity("".join(self.helix.values()))
        sheet_hydro = self._calc_hydrophobicity("".join(self.sheet.values()))
        turn_hydro = self._calc_hydrophobicity("".join(self.turn.values()))

        return f"Total Hydrophobicity: {total_hydro:.2f}\n" \
               f"Helix Hydrophobicity: {helix_hydro:.2f}\n" \
               f"Sheet Hydrophobicity: {sheet_hydro:.2f}\n" \
               f"Turn Hydrophobicity: {turn_hydro:.2f}\n"

    def preferences(self) -> str:
        """
        Get all the preferences of the amino acids in the secondary structures.
        :return: A formatted string of preferences.
        """
        preference_string = ""
        structure_names = ["HELIX", "SHEET", "TURN"]
        for index, preference in enumerate(self._calc_preference()):
            name = structure_names[index]
            structure_preference = ", ".join([acid for acid, _ in preference])
            preference_string += f"{name}:\n{structure_preference}\n"

        return f"Preferences: \n{preference_string}"

    def analyze(self):
        return (f"{self.get_protein_mass()}"
                f"{self.get_atom_mass()}"
                f"{self.hydrophobicity()}\n"
                f"{self.preferences()}")

    def compositions(self):
        total_comp = self._calc_composition(self.sequence)
        helix_comp = self._calc_composition("".join(self.helix.values()))
        sheet_comp = self._calc_composition("".join(self.sheet.values()))
        turn_comp = self._calc_composition("".join(self.turn.values()))

        return f"Total Composition:\n{self._histogram(total_comp)}\n" \
               f"Helix Composition:\n{self._histogram(helix_comp)}\n" \
               f"Sheet Composition:\n{self._histogram(sheet_comp)}\n" \
               f"Turn Composition:\n{self._histogram(turn_comp)}\n"

    def compare_sequences(self, fasta_file_path: Path):
        """
        Compares the PDB Sequence with a Fasta sequence.
        :param fasta_file_path: Path to a fasta file.
        :return: A formatted human-readable sequence comparison. Where a bar(|) indicates a match.
        """
        # Load fasta Sequence.
        fasta_sequence = self._load_fasta_sequence(fasta_file_path)

        # Find the longest sequence for formatting.
        match_length = max(len(self.sequence), len(fasta_sequence))
        match_format = list(repeat(" ", match_length))

        # Determan the line length to use.
        line_length = self.max_line_length
        if self.max_line_length <= 0:
            line_length = match_length

        # Add a link for each match.
        matches = find_matches(self.sequence, fasta_sequence)
        for match in matches:
            match_format[match] = "|"
        match_format = "".join(match_format)

        formatted_sequence = (self.sequence, match_format, fasta_sequence)

        # Split into segments when the lines are too long.
        if match_length > line_length:
            # Split the lines into line_length long segments.
            lines = []
            for line in formatted_sequence:
                split_line = split_str(line, line_length, "\n")
                lines.append(split_line)

            # Make zipped sets for the lines and join these together.
            formatted_sequence = []
            for line in zip(*lines):
                formatted_sequence.append("".join(line))

        formatted_sequence = "\n".join(formatted_sequence)

        return f"Compared Sequence: \n{formatted_sequence}"
