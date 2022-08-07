import enum
import numpy as np


class LatticeSystem(enum.Enum):
    TRICLINIC = (1, "Triclinic", 6)
    MONOCLINIC = (2, "Monoclinic", 4)
    ORTHORHOMBIC = (3, "Orthorhombic", 3)
    HEXAGONAL = (4, "Hexagonal", 2)
    RHOMBOHEDRAL = (5, "Rhombohedral", 2)
    TETRAGONAL = (6, "Tetragonal", 2)
    CUBIC = (7, "Cubic", 1)

    def __init__(self, id, title, d):
        self.id = id
        self.title = title
        self.d = d


root = "/home/daniil/PycharmProjects/1prog/indexing/icsd-new/systems/"
files = {lattice_system: "{}/{}.txt".format(root, lattice_system.name.lower()) for lattice_system in LatticeSystem}


def get_data(lattice_system, unpack=True):
    """
    :param lattice_system:
    :param unpack:
    :return:
    tuple of numpy arrays
    angles are in degrees
    TRICLINIC
    a, b, c, alpha, beta, gamma
    MONOCLINIC
    a, b, c, beta
    ORTHORHOMBIC
    a, b, c
    HEXAGONAL
    a, c
    RHOMBOHEDRAL
    a, c
    TETRAGONAL
    a, c
    CUBIC
    a
    """
    data = np.loadtxt(files[lattice_system], unpack=unpack)
    if (lattice_system is LatticeSystem.CUBIC) and unpack:
        # to keep uniform: always a tuple
        return data,
    else:
        return data


lengths = ("a", "b", "c")
angles = ("alpha", "beta", "gamma")

parameter_names = ("a", "b", "c", "alpha", "beta", "gamma")

independent_parameter_names = {
    LatticeSystem.TRICLINIC: ("a", "b", "c", "alpha", "beta", "gamma"),
    LatticeSystem.MONOCLINIC: ("a", "b", "c", "beta"),
    LatticeSystem.ORTHORHOMBIC: ("a", "b", "c"),
    LatticeSystem.HEXAGONAL: ("a", "c"),
    LatticeSystem.RHOMBOHEDRAL: ("a", "c"),
    LatticeSystem.TETRAGONAL: ("a", "c"),
    LatticeSystem.CUBIC: ("a", )
}


def get_independent_parameter_names(lattice_system):
    return independent_parameter_names[lattice_system]


parameter_groups = {
    LatticeSystem.TRICLINIC: (["a", "b", "c"], ["alpha", "beta", "gamma"]),
    LatticeSystem.MONOCLINIC: (["a", "c"], ["b"], ["beta"]),
    LatticeSystem.ORTHORHOMBIC: (["a", "b", "c"],),
    LatticeSystem.HEXAGONAL: (["a"],  ["c"]),
    LatticeSystem.RHOMBOHEDRAL: (["a"], ["c"]),
    LatticeSystem.TETRAGONAL: (["a"], ["c"]),
    LatticeSystem.CUBIC: (["a"],)
}


def get_parameter_groups(lattice_system):
    return parameter_groups[lattice_system]
