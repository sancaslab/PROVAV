import unittest
from proteinStructure import proteinStructure
from Bio.PDB.Structure import Structure

class proteinStructure_Case(unittest.TestCase):

    def test_ReadingPdbFileShouldReturnTrue(self):
        structure = proteinStructure("2kpa", "PDB_proteins_files/2kpa.pdb")
        self.assertFalse(type(structure) == type(None))

    def test_ReadingCifFileShouldReturnTrue(self):
        structure = proteinStructure("6si7", "PDB_proteins_files/6si7.cif")
        self.assertFalse(type(structure) == type(None))

    def test_NotReadingFileshouldReturnTrue(self):
        structure = proteinStructure("6si7", "PDB_proteins_files/6si7.xml")
        self.assertFalse(structure is None)

    def test_FileNotFoundShouldPrintFileNotFound(self):
        structure = proteinStructure("2kpa", "PDB_proteins_files/2kp.pdb")
        self.assertFalse(type(structure) == type(None))

if __name__ == '__main__':
    unittest.main()
