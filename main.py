import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pyvista as pv
from pyvista import examples

mesh = examples.load_airplane()
print('cells:', mesh.n_cells)
print('points:', mesh.n_points)
print('n_arrays:', mesh.n_arrays)
print('bounds:', mesh.bounds)
print('center:', mesh.center)