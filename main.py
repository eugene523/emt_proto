import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pyvista as pv
from pyvista import examples

points = np.array([[0, 0, 0]], dtype=float)
point_cloud = pv.PolyData(points)
point_cloud.plot(eye_dome_lighting=True)