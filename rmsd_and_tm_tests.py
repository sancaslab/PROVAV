import unittest
import warnings
from PDB_reader import PDB_reader
from root_mean_squared_deviation_and_tm import root_mean_squared_deviation, tm_score
from Bio.PDB.PDBExceptions import PDBConstructionWarning


# It's normal to get a lot of variations from these files -> for that we omit deprecation and constructing warning


class RMSD_Test_Case(unittest.TestCase):
    # Testing root mean squared deviation

    # Creatine Kinase from Human Muscle - Creatine Kinase from Human Brain
    # High number of atoms variation -> this may add changes
    def test_CKHB_CKHM(self):
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        structure_x = PDB_reader("3b6r", "PDB_proteins_files/3b6r.pdb").read()
        structure_y = PDB_reader("1i0e", "PDB_proteins_files/1i0e.pdb").read()
        rmsd = root_mean_squared_deviation(structure_x, structure_y).compare_rmsd_models()
        res = False
        # pyMol = 92.95 -> Same type of proteins, same specialized function in different places,
        # however they have very diverse structures
        # 92.72
        if 90 <= rmsd <= 93:
            res = True

        self.assertEqual(res, True)

    # V-ATPase a2-subunit isoforms
    # Type H+ -> Vacuolate proton pumping
    # mammalian proteins that share a highly conserved structure
    def test_2KPA_2KPB(self):
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        structure_x = PDB_reader("2kpa", "PDB_proteins_files/2kpa.pdb").read()
        structure_y = PDB_reader("2kpb", "PDB_proteins_files/2kpb.pdb").read()
        rmsd = root_mean_squared_deviation(structure_x, structure_y).compare_rmsd_models()
        res = False
        # pyMol = 3.46
        # RMSD =  2.04 (software)
        # 3.92

        if 3 <= rmsd <= 5:
            res = True

        self.assertEqual(res, True)

    # Isoforms of Adenylyl cyclase bound to protein G
    def test_ACG1_ACG2(self):
        warnings.filterwarnings("ignore", category=PDBConstructionWarning)
        structure_x = PDB_reader("6r3q", "PDB_proteins_files/6r3q.pdb").read()
        structure_y = PDB_reader("6r4o", "PDB_proteins_files/6r4o.pdb").read()
        rmsd = root_mean_squared_deviation(structure_x, structure_y).compare_rmsd_models()
        res = False
        # pyMol = 1.61
        # RMSD =  1.98 (software)
        # 2.46
        if 1 <= rmsd <= 3:
            res = True

        self.assertEqual(res, True)


class TM_Test_Case(unittest.TestCase):
    # Testing tm scoring

    def test_ACG1_ACG2(self):
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        # TM - score = 0.97
        # 0.97 (now)
        tm = tm_score("PDB_proteins_files/6r3q.pdb", "PDB_proteins_files/6r4o.pdb").get_tm()
        res = False
        # We start with random structural similarity
        if 0.5 < tm < 1.00:
            # Proteins in the same fold
            res = True
        self.assertEqual(res, True)

    def test_CKHB_CKHM(self):
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        # TM - score = 0.9453
        # 0.97 (now)
        tm = tm_score("PDB_proteins_files/3b6r.pdb", "PDB_proteins_files/1i0e.pdb").get_tm()
        res = False
        # We start with random structural similarity
        if 0.5 < tm < 1.00:
            # Proteins in the same fold
            res = True
        self.assertEqual(res, True)

    # Human Leukocyte Antigen -> some chains of major histocompatibility complex (MHC)
    def test_HLA_proteins(self):
        warnings.filterwarnings("ignore", category=DeprecationWarning)
        # TM - score =  0.5078
        # 0.53 (now)
        tm = tm_score("PDB_proteins_files/1hdm.pdb", "PDB_proteins_files/4i0p.pdb").get_tm()
        res = False
        # We start with random structural similarity
        if 0.5 < tm < 1.00:
            # Proteins in the same fold
            res = True
        self.assertEqual(res, True)


if __name__ == '__main__':
    unittest.main()
