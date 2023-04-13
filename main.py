import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

import pyvista as pv
from pyvista import examples
import vtk

vtk_array = vtk.vtkDoubleArray()
vtk_array.SetNumberOfComponents(3)
vtk_array.SetNumberOfValues(9)
vtk_array.SetValue(0, 0)
vtk_array.SetValue(1, 0)
vtk_array.SetValue(2, 0)
vtk_array.SetValue(3, 1)
vtk_array.SetValue(4, 0)
vtk_array.SetValue(5, 0)
vtk_array.SetValue(6, 0.5)
vtk_array.SetValue(7, 0.667)
vtk_array.SetValue(8, 0)

np_points = np.array([[0, 0, 0],
                      [1, 0, 0],
                      [0.5, 0.667, 0]])

cells = [3, 0, 1, 2]

mesh = pv.PolyData(np_points, cells)

# mesh.plot(cpos='xy', show_edges=True)

grid = pv.UniformGrid(dimensions=(3, 3, 1))
ugrid = grid.cast_to_unstructured_grid()
pl = pv.Plotter()
pl.add_mesh(ugrid, show_edges=True, line_width=5)
label_coords = ugrid.points + [0, 0, 0.02]
point_labels = [f'Point {i}' for i in range(ugrid.n_points)]
pl.add_point_labels(label_coords, point_labels, font_size=25, point_size=20)
cell_labels = [f'Cell {i}' for i in range(ugrid.n_cells)]
pl.add_point_labels(ugrid.cell_centers(), cell_labels, font_size=25)
pl.camera_position = 'xy'

simple_range = range(ugrid.n_cells)
ugrid.cell_data['my-data'] = simple_range

ugrid.plot(cpos='xy', show_edges=True)