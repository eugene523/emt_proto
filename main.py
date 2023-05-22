
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pyvista as pv
from pyvista import examples
import vtk

from panel import Panel

import demo

#demo.demo_1x1()

from scipy import sparse

I = np.array([0,0,1,3,1,0,0])
J = np.array([0,2,1,3,1,0,0])
V = np.array([1,1,1,1,1,1,1])
B = sparse.coo_matrix((V,(I,J)),shape=(4,4)).tocsr()
print(B.toarray())