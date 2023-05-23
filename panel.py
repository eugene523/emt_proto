from enum import Enum
import math
from math_utils import CooBuilder
from shellmat import ShellMaterial
from meshing import *
from validation import *
from fea import Fe3
import numpy as np
import pyvista


class NodeGroup(Enum):
    N00 = 1
    N01 = 2
    N10 = 3
    N11 = 4
    LFT = 5
    RGT = 6
    TOP = 7
    BOT = 8


class Panel:
    ERR_LENGTH_NOT_SET      = "Panel length is not setted."
    ERR_WIDTH_NOT_SET       = "Panel width is not setted."
    ERR_ELEM_LENGTH_NOT_SET = "Mesh element length is not setted."

    def __init__(self, length: float, width: float):
        self.length:      float = length
        self.width:       float = width
        self.material:    ShellMaterial = None
        self.elem_length: float = 0.0
        self.mesh:        Mesh = None
        self.node_groups: dict[NodeGroup, list[Node]] = None
        self.constraints: dict[NodeGroup, ConstraintVector] = {}
        self.forces:      dict[NodeGroup, ForceVector] = {}
        self.fin_elems:   list[Fe3] = []
        self.cb:          CooBuilder = None

        # --- fea results --- #
        self.disp_vector: np.ndarray = None

    def validate_before_meshing(self) -> ValidationResult:
        v = ValidationResult()

        if self.length == 0:
            v.add_issue(self.ERR_LENGTH_NOT_SET)

        if self.width == 0:
            v.add_issue(self.ERR_WIDTH_NOT_SET)

        if self.elem_length == 0:
            v.add_issue(self.ERR_ELEM_LENGTH_NOT_SET)

        return v
        
    def __get_n_elems_on_edge(self, edge_len: float) -> int:
        n = math.ceil(edge_len / self.elem_length)
        return n

    def __create_node_groups(self):
        self.node_groups = {}
        L = self.length
        W = self.width
        eps = 1e-6

        self.node_groups[NodeGroup.N00] = self.mesh.select_nodes_near_point(0, 0, 0, eps)
        self.node_groups[NodeGroup.N01] = self.mesh.select_nodes_near_point(0, W, 0, eps)
        self.node_groups[NodeGroup.N10] = self.mesh.select_nodes_near_point(L, 0, 0, eps)
        self.node_groups[NodeGroup.N11] = self.mesh.select_nodes_near_point(L, W, 0, eps)

        self.node_groups[NodeGroup.LFT] = self.mesh.select_nodes_on_edge(0, 0, 0, 
                                                                         0, W, 0,
                                                                         eps)
        
        self.node_groups[NodeGroup.RGT] = self.mesh.select_nodes_on_edge(L, 0, 0, 
                                                                         L, W, 0,
                                                                         eps)
        
        self.node_groups[NodeGroup.TOP] = self.mesh.select_nodes_on_edge(0, W, 0, 
                                                                         L, W, 0,
                                                                         eps)
        
        self.node_groups[NodeGroup.BOT] = self.mesh.select_nodes_on_edge(0, 0, 0, 
                                                                         L, 0, 0,
                                                                         eps)
        
    def do_mesh(self):
        v = self.validate_before_meshing()
        if not v.is_ok():
            raise Exception(v)
        
        L = self.length
        W = self.width
        q = Quad.new_by_coord(0, 0, 0,
                              0, W, 0,
                              L, W, 0,
                              L, 0, 0)
        
        n_len = self.__get_n_elems_on_edge(self.length)
        n_wid = self.__get_n_elems_on_edge(self.width)
        self.mesh = q.mesh_tria(n_len, n_wid, 1)
        self.__create_node_groups()

    def set_constraint(self, node_group: NodeGroup, constraint: Constraint):
        self.constraints[node_group] = constraint

    def set_force(self, node_group: NodeGroup, force: ForceVector):
        self.forces[node_group] = force

    def compute(self):
        self.__create_finite_elements()
        self.__create_stiffeness_matrix()
        self.__apply_constraints_to_stiffeness_matrix()
        self.__create_force_vector()
        self.__solve_disp()

    def __apply_constraints_to_nodes(self):
        for node_group_key in self.constraints:

        
    def __create_finite_elements(self):
        self.fin_elems.clear()
        for elem in self.mesh.elements:
            fe3_elem = Fe3(elem.i, elem.j, elem.k)
            fe3_elem.set_material(self.material)
            fe3_elem.compute()
            self.fin_elems.append(fe3_elem)

    def __create_stiffeness_matrix(self):
        n_elems = len(self.fin_elems)
        max_arr_size = n_elems * 36
        n_indeces = 3
        block_size = 2
        self.cb = CooBuilder(max_arr_size, n_indeces, block_size)
        for elem in self.fin_elems:
            self.cb.accept_matrix(elem.k_mbr_6x6)

    def __apply_constraints_to_stiffeness_matrix(self):
        pass


    def __get_fixed_dofs(self):
        constrained_nodes = []
        for node_group_key in self.constraints:
            constrained_nodes.extend(self.node_groups[node_group_key])

        for node in constrained_nodes:

        fixed_dofs = []

    def __create_force_vector(self):
        pass

    def __solve_disp(self):
        pass
    
    def show_just_mesh(self):
        assert self.mesh != None
        nodes = self.mesh.nodes
        n_nodes = self.mesh.get_n_nodes()
        points = np.zeros((n_nodes, 3), dtype=float)
        for i in range(n_nodes):
            node = nodes[i]
            index = node.index
            points[index, 0] = node.x
            points[index, 1] = node.y
            points[index, 2] = node.z

        n_elems = self.mesh.get_n_elements()
        cells = np.zeros((n_elems, 4), dtype=int)
        for i in range(n_elems):
            elem = self.mesh.elements[i]
            cells[i, 0] = 3
            cells[i, 1] = elem.i.index
            cells[i, 2] = elem.j.index
            cells[i, 3] = elem.k.index
        
        pv_mesh = pyvista.PolyData(points, cells)
        pl = pyvista.Plotter()
        pl.add_mesh(pv_mesh, show_edges=True, line_width=1)
        pl.camera_position = 'xy'
        pl.show_bounds()
        pl.show()

    
    def show(self):
        assert self.mesh != None
        nodes = self.mesh.nodes
        n_nodes = self.mesh.get_n_nodes()
        points = np.zeros((n_nodes, 3), dtype=float)
        labels = [None] * n_nodes
        for i in range(n_nodes):
            node = nodes[i]
            index = node.index
            points[index, 0] = node.x
            points[index, 1] = node.y
            points[index, 2] = node.z
            labels[index] = str(index)

        n_elems = self.mesh.get_n_elements()
        cells = np.zeros((n_elems, 4), dtype=int)
        for i in range(n_elems):
            elem = self.mesh.elements[i]
            cells[i, 0] = 3
            cells[i, 1] = elem.i.index
            cells[i, 2] = elem.j.index
            cells[i, 3] = elem.k.index
        
        pv_mesh = pyvista.PolyData(points, cells)
        # pv_mesh.plot(show_bounds=True, cpos='xy', show_edges=True)

        pl = pyvista.Plotter()
        pl.add_mesh(pv_mesh, show_edges=True, line_width=1)
        # pl.add_point_labels(pv_mesh.points, labels, font_size=10, point_size=10)
        pl.camera_position = 'xy'
        pl.show_bounds()
        pl.show()
    