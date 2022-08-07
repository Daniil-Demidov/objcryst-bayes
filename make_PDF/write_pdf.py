import numpy as np
import global_data
from global_data import LatticeSystem
import split_dist
from pathlib import Path
import statsmodels.api as sm


thresholds = {
    LatticeSystem.TRICLINIC:
        {"a": 20},
    LatticeSystem.MONOCLINIC:
        {"a": 20, "b": 15},
    LatticeSystem.ORTHORHOMBIC:
        {"a": 20},
    LatticeSystem.HEXAGONAL:
        {"a": 15, "c": 50},
    LatticeSystem.RHOMBOHEDRAL:
        {"a": 25, "c": 90},
    LatticeSystem.TETRAGONAL:
        {"a": 15, "c": 25},
    LatticeSystem.CUBIC:
        {"a": 20}
}

for lattice_system in LatticeSystem:
    print(lattice_system.name)
    Path("PDF/{}".format(lattice_system.name)).mkdir(parents=True, exist_ok=True)
    all_params = global_data.get_data(lattice_system)
    param_names = global_data.get_independent_parameter_names(lattice_system)
    dict_params = {}
    for param_name, data in zip(param_names, all_params):
        dict_params[param_name] = data
    for par_group in global_data.get_parameter_groups(lattice_system):
        first_par_name = par_group[0]
        print("{}: {}".format(lattice_system.name, first_par_name))
        if first_par_name in global_data.angles:
            data = np.concatenate([dict_params[name] for name in par_group])
            data = np.delete(data, data == 90)
            data = np.maximum(180 - data, data)
            data = np.delete(data, data > 120)
            data = np.concatenate((data, 180 - data))
            data = np.cos(np.radians(data))
            dist = sm.nonparametric.KDEUnivariate(data)
            dist.fit(adjust=0.1)
            dx = 0.001
            xs = np.arange(0, 0.55, dx)
            dens = dist.evaluate(xs)
        else:
            data = np.concatenate([dict_params[name] for name in par_group])
            thresh = thresholds[lattice_system][first_par_name]
            dist = split_dist.fit_dist(data, thresh)
            dx = 0.05
            xs = np.arange(0, data.max() + 2*dx, dx)
            dens = dist.pdf(xs)
            min_x = 0.95*data.min()
            opt_x = xs[dens.argmax()]
            doubtful_mask = np.logical_and(min_x <= xs, xs < opt_x)
            dens[doubtful_mask] = np.maximum(dens[doubtful_mask], 0.005*dens.max())
        if first_par_name in global_data.angles:
            first_par_title = "cos_" + first_par_name
        else:
            first_par_title = first_par_name
        file = "PDF/{}/{}-PDF.txt".format(lattice_system.name, first_par_title)
        header = "{}, density. dx={}, size={}".format(first_par_title, dx, xs.size)
        np.savetxt(file, np.vstack((xs, dens)).T, header=header)
