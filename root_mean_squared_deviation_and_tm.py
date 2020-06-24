import numpy
import tmscoring
from Bio.PDB import Structure
from Bio.PDB.PDBExceptions import PDBConstructionWarning
from Bio.SVDSuperimposer import SVDSuperimposer
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=PDBConstructionWarning)


# Method that calculate rmsd with SVDSuperimposer BioPython library
# This function has two atoms arrays like parameters from which it will
# calculate the best position to get the value
def calculate_rmsd(atoms_x_coord, atoms_y_coord) -> float:
    super_imposer = SVDSuperimposer()  # Calling the class that superposes atoms arrays
    super_imposer.set(atoms_x_coord, atoms_y_coord)  # Vector y will be rotated and translated on vector x
    super_imposer.run()
    value = super_imposer.get_rms()  # Get the value of RMSD

    return value


# Method to get atoms coordinates (are situated in CA residue, alpha carboxyl)
def get_alpha_carboxyl_atoms_coord(model):
    atoms = []
    for chain in model:
        for res in chain:
            for atom in res:
                if atom.get_name() == "CA":
                    atoms.append(atom.get_coord())
    return numpy.array(atoms)


# All atoms -> to obtain some functionality tests
""""
def get_atoms_coord(model):
    atoms = []
    for chain in model:
        for res in chain:
            for atom in res:
                atoms.append(atom.get_coord())
    return numpy.array(atoms)
"""


# We want to get the same number of coordinates in both vectors
def get_atoms(atoms_x, atoms_y):
    res = []
    if len(atoms_x) > len(atoms_y):
        new_atoms_x = atoms_x[:len(atoms_y)]
        res.append(new_atoms_x)
        res.append(atoms_y)
    else:
        new_atoms_y = atoms_y[:len(atoms_x)]
        res.append(atoms_x)
        res.append(new_atoms_y)

    return res


class root_mean_squared_deviation:

    # Class constructor
    # It receives two proteins structures
    def __init__(self, structure_1: Structure, structure_2: Structure):
        self.structure_x = structure_1
        self.structure_y = structure_2

    # Function that calls other methods to get the best value in all the models we find in structure protein
    def compare_rmsd_models(self) -> float:
        atoms_x = []
        atoms_y = []
        total_rmsd = []
        # For each model we get in the proteins, we combine all possibilities, get atoms coordinates
        # and calculate rmsd
        for model_x in self.structure_x:
            for model_y in self.structure_y:
                atoms_x.append(get_alpha_carboxyl_atoms_coord(model_x))
                atoms_y.append(get_alpha_carboxyl_atoms_coord(model_y))
        # Selection of several atoms when arrays have different sizes
        for i in range(len(atoms_x)):
            if len(atoms_x[i]) != len(atoms_y[i]):
                new_atoms = get_atoms(atoms_x[i], atoms_y[i])
                atoms_x[i] = new_atoms[0]
                atoms_y[i] = new_atoms[1]

            rmsd = calculate_rmsd(atoms_x[i], atoms_y[i])
            total_rmsd.append(rmsd)

        # We select the best score
        # print("Min RMSD obtained: " + str(min(total_rmsd)))
        return min(total_rmsd)


class tm_score:
    # Class constructor that receives like parameters two files names
    def __init__(self, file_name_1: str, file_name_2: str):
        self.file_1 = file_name_1
        self.file_2 = file_name_2

    # We call python tmscoring to get best alignment and
    # It calculates tm value
    def get_tm(self) -> float:
        alignment = tmscoring.TMscoring(self.file_1, self.file_2)
        alignment.optimise()
        tm = alignment.tmscore(**alignment.get_current_values())
        # tm = tmscoring.get_tm(file_pdb_1, file_pdb_2)
        # rmsd = tmscoring.get_rmsd(pdb1[0][0], pdb2[0][0])

        return tm

