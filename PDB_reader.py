# pdb file lecture
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Structure import Structure


class PDB_reader:

    def __init__(self, protein_id: str, file: str):
        self.protein_name = protein_id
        self.file_name = file

    def read(self) -> Structure:
        # try:
        parser = PDBParser(PERMISSIVE=1)
        structure: Structure = parser.get_structure(self.protein_name, self.file_name)
        # except OSError:
        # print("Can´t find file path")
        return structure


# print(proteinStructure[0])
# structure_x1 = PDB_reader("2KPA", "db").read()
