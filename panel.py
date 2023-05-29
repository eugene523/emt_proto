from enum import Enum
import math
from math_utils import SpBuilder
from shellmat import ShellMaterial
from meshing import *
from validation import *
from fea import Fe3
import numpy as np
from scipy.sparse.linalg import spsolve
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
        self.dof = 2

        # --- fea data and results --- #
        self.fin_elems:   list[Fe3] = None
        self.fixed_dofs:  list[int] = None
        self.sp_builder:  SpBuilder = None

        # [K] * {F} = {D}
        self.k_glob = None # stiffeness matrix of the whole construction
        self.f_glob: np.ndarray = None # force vector
        self.d_glob: np.ndarray = None # displacement vector
        self.disp_mag: np.ndarray = None # displacement magnitudes

    def set_material(self, material: ShellMaterial):
        self.material = material

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
        self.__apply_constraints_to_nodes()
        self.__apply_forces_to_nodes()
        self.__create_fixed_dofs_list()
        self.__create_finite_elements()
        self.__create_global_stiffeness_matrix()
        self.__create_global_force_vector()
        self.__solve_disp()
        self.__compute_disp_magnitudes()

    def __apply_constraints_to_nodes(self):
        self.mesh.clear_constraints()
        for node_group_key in self.constraints:
            constraint = self.constraints[node_group_key]
            nodes = self.node_groups[node_group_key]
            for node in nodes:
                node.constraint_superposition(constraint)

    def __apply_forces_to_nodes(self):
        self.mesh.clear_forces()
        for node_group_key in self.forces:
            force = self.forces[node_group_key]
            nodes = self.node_groups[node_group_key]
            for node in nodes:
                node.force_superposition(force)
        
    def __create_fixed_dofs_list(self):
        fixed_dofs = []
        for node in self.mesh.nodes:
            if not node.is_free():
                i = node.index
                p = self.dof * i
                c = node.constraint_vector

                if c.tx == Constraint.FIXED:
                    fixed_dofs.append(p)

                if c.ty == Constraint.FIXED:
                    fixed_dofs.append(p + 1)

                """ if c.tz == Constraint.FIXED:
                    fixed_dofs.append(p + 2)
                    
                if c.rx == Constraint.FIXED:
                    fixed_dofs.append(p + 3)

                if c.ry == Constraint.FIXED:
                    fixed_dofs.append(p + 4)
                    
                if c.rz == Constraint.FIXED:
                    fixed_dofs.append(p + 5) """

        self.fixed_dofs = fixed_dofs


    def __create_finite_elements(self):
        assert(self.material != None)
        fin_elems = []
        for elem in self.mesh.elements:
            fe3_elem = Fe3(elem.i, elem.j, elem.k)
            fe3_elem.set_material(self.material)
            fe3_elem.compute()
            fin_elems.append(fe3_elem)

        self.fin_elems = fin_elems


    def __create_global_stiffeness_matrix(self):
        # Building global stiffeness matrix
        n_elems = len(self.fin_elems)
        n_indeces = 3 # number of nodes per element
        block_size = 2
        mat_size = 6
        mat_n_entries = mat_size ** 2
        max_arr_size = n_elems * mat_n_entries
        n_nodes = self.mesh.get_n_nodes()
        sp_size = n_nodes * self.dof
        sp_builder = SpBuilder(max_arr_size, n_indeces, block_size, sp_size)
        for elem in self.fin_elems:
            sp_builder.accept_matrix(elem.k_mbr_6x6, elem.get_index_vector())

        # Applying constraints to global stiffeness matrix.
        # All fixed degrees of freedom are indeces of rows and columns.
        # We zero out this rows and columns and place 1.0 at position k[i, i],
        # where k is global stiffeness matrix, i is index.
        for i in self.fixed_dofs:
            sp_builder.set_row_col_to_zero_and_place_1(i)

        self.k_glob = sp_builder.get_csr()


    def __create_global_force_vector(self):
        # Building global force vector
        n_nodes = self.mesh.get_n_nodes()
        v_size = n_nodes * self.dof
        f_glob = np.zeros((v_size, 1), dtype=float)
        for node in self.mesh.nodes:
            p = node.index * self.dof
            f = node.force_vector

            if f.fx != 0:
                f_glob[p, 0] += f.fx

            if f.fy != 0:
                f_glob[p + 1, 0] += f.fy

            """ if f.fz != 0:
                f_glob[p + 2, 0] += f.fz

            if f.mx != 0:
                f_glob[p + 3, 0] += f.mx

            if f.my != 0:
                f_glob[p + 4, 0] += f.my
            
            if f.mz != 0:
                f_glob[p + 5, 0] += f.mz """
        
        # Applying constraints to global force vector.
        # All fixed degrees of freedom are indeces of components in force vector
        # We zero out this components.
        for i in self.fixed_dofs:
            f_glob[i, 0] = 0.0

        self.f_glob = f_glob


    def __solve_disp(self):
        self.d_glob = spsolve(self.k_glob, self.f_glob)

    def __compute_disp_magnitudes(self):
        disp_mag = np.zeros((2, 2), dtype=float)

        for node in self.mesh.nodes:
            p = node.index * self.dof
            for i in range(self.dof):
                v = self.d_glob[p + i]
                disp_mag[i, 0] = min(disp_mag[i, 0], v)
                disp_mag[i, 1] = max(disp_mag[i, 1], v)

        self.disp_mag = disp_mag


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

    def show_static_deform_mesh(self):
        assert self.mesh != None
        nodes = self.mesh.nodes
        n_nodes = self.mesh.get_n_nodes()
        points = np.zeros((n_nodes, 3), dtype=float)
        for i in range(n_nodes):
            node = nodes[i]
            index = node.index
            p = index * self.dof
            points[index, 0] = node.x + self.d_glob[p + 0]
            points[index, 1] = node.y + self.d_glob[p + 1]
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
    