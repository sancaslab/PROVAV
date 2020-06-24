from Bio.PDB import PDBParser, MMCIFParser
#This package was included in the initial version but had to be commented due to a problem with the installation
#from Bio.PDB.mmtf import MMTFParser


class proteinStructure:

    # default constructor
    def __init__(self, proteinId, proteinFile):

        if proteinFile.endswith(".cif"):
            structure = self.cifReader(proteinId, proteinFile)
        elif proteinFile.endswith(".ent") or proteinFile.endswith(".pdb"):
            structure = self.pbdReader(proteinId, proteinFile)
        #elif proteinFile.endswith(".mmtf"):
        #    structure = self.mmtfReader(proteinFile)


    def cifReader(self, name, file):
        try:
            parser = MMCIFParser()
            structure = parser.get_structure(name, file)
            return structure
        except:
            print("Something went wrong: File not found")


    def pbdReader(self, name, file):
        try:
            parser = PDBParser(PERMISSIVE=1)
            structure = parser.get_structure(name, file)
            return structure
        except:
            print("Something went wrong: File not found")

    #def mmtfReader(self, file):
     #   try:
      #      structure=MMTFParser.get_structure(file)
       #     return structure
        #except:
         #   print("Something went wrong: File not found")