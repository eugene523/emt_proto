
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pyvista as pv
from pyvista import examples
import vtk

from panel import Panel

p = Panel(200, 100)
p.elem_length = 1
p.do_mesh()
p.show()