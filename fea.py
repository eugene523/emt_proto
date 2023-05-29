from enum import Enum
import math
import math_utils
import shellmat
from boundary import *
from meshing import *
import numpy as np


class Fe3:
    def __init__(self, i: Node, j: Node, k: Node):
        self.i = i
        self.j = j
        self.k = k

        self.i_glb = i.get_coord_vect()
        self.j_glb = j.get_coord_vect()
        self.k_glb = k.get_coord_vect()

        # Матрицы поворота (направляющих косинусов) элемента 
        # [3x3], [9x9], [18x18], [24x24].
        self.r_3x3:   np.ndarray = None
        self.r_9x9:   np.ndarray = None
        self.r_18x18: np.ndarray = None
        self.r_24x24: np.ndarray = None

        # Локальные координаты узлов.
        self.i_loc: np.ndarray = None
        self.j_loc: np.ndarray = None
        self.k_loc: np.ndarray = None

        self.area: float = 0.0
        '''Площадь элемента.'''

        self.phi: float = 0.0
        '''Угол наклона для перевода матриц упругости (в радианах).'''

        # Матрицы для поворота векторов напряжений и деформаций.
        self.t1_3x3:    np.ndarray = None
        self.t1_tr_3x3: np.ndarray = None
        self.t2_3x3:    np.ndarray = None
        self.t2_tr_3x3: np.ndarray = None

        # ---------------------------------------------------------------------

        self.material: shellmat.ShellMaterial = None

        self.mbr_glb_3x3: np.ndarray = None
        '''Матрица упругости материала [3x3] для мембранной компоненты в ГСК.'''

        self.mbr_loc_3x3: np.ndarray = None
        '''Матрица упругости материала [3x3] для мембранной компоненты в ЛСК.'''

        self.bnd_glb_3x3: np.ndarray = None
        '''Матрица упругости материала [3x3] для изгибной компоненты в ГСК.'''

        self.bnd_loc_3x3: np.ndarray = None
        '''Матрица упругости материала [3x3] для изгибной компоненты в ЛСК.'''

        # ---------------------------------------------------------------------

        self.b_mbr_3x6: np.ndarray = None
        '''Матрица градиентов [3x6] для мембранной компоненты.'''

        self.k_mbr_6x6: np.ndarray = None
        '''Матрица жесткости  [6x6] для мембранной компоненты в ЛСК.'''

        self.b_bnd_3x9_list: list[np.ndarray] = None
        '''Матрицы градиентов [3x9] для задачи изгиба.'''
    
    def get_index_vector(self) -> list[int]:
        return [self.i.index, self.j.index, self.k.index]

    def set_material(self, material: shellmat.ShellMaterial):
        self.material = material
        self.mbr_glb_3x3 = material.get_mbr_3x3()
        self.bnd_glb_3x3 = material.get_bnd_3x3()

    def compute(self):
        self.__compute_area()
        assert self.area != 0.0

        self.__compute_rotation_matrices()
        self.__compute_local_coordinates()
        self.__compute_phi()
        self.__compute_t_matrices_3x3()
        self.__compute_mbr_bnd_loc_3x3()
        self.__compute_b_mbr_3x6()
        self.__compute_k_mbr_6x6()

    def __compute_area(self):
        self.area = math_utils.get_tria_area(self.i.x, self.i.y, self.i.z,
                                             self.j.x, self.j.y, self.j.z,
                                             self.k.x, self.k.y, self.k.z)

    def __compute_rotation_matrices(self):
        v_ij = self.j_glb - self.i_glb
        v_ik = self.k_glb - self.i_glb
        
        v_x = v_ij / np.linalg.norm(v_ij)

        v_z = np.cross(v_ij, v_ik)
        v_z = v_z / np.linalg.norm(v_z)

        v_y = np.cross(v_z, v_x)

        self.r_3x3   = np.array([v_x, v_y, v_z], dtype = float)
        self.r_9x9   = np.kron(np.eye(3), self.r_3x3)
        self.r_18x18 = np.kron(np.eye(6), self.r_3x3)
        self.r_24x24 = np.kron(np.eye(8), self.r_3x3)

    def __compute_local_coordinates(self):
        self.i_loc = np.matmul(self.r_3x3, self.i_glb)
        self.j_loc = np.matmul(self.r_3x3, self.j_glb)
        self.k_loc = np.matmul(self.r_3x3, self.k_glb)

    def __compute_phi(self):
        v_ij = self.j_glb - self.i_glb
        ax = v_ij[0]
        ay = v_ij[1]
        self.phi = math.atan2(ay, ax)
        
    def __compute_t_matrices_3x3(self):
        s = math.sin(self.phi)
        c = math.cos(self.phi)

        c2    = c ** 2
        s2    = s ** 2
        sc    = s * c
        _2sc  = 2 * sc
        c2_s2 = c2 - s2

        t1_3x3 = np.array([
            [c2,  s2, -_2sc ],
            [s2,  c2,  _2sc ],
            [sc, -sc,  c2_s2]],
            dtype=float)
        
        t2_3x3 = np.array([
            [c2,    s2,  -sc   ],
            [s2,    c2,   sc   ],
            [_2sc, -_2sc, c2_s2]],
            dtype=float)
        
        self.t1_3x3    = t1_3x3
        self.t1_tr_3x3 = t1_3x3.transpose()
        self.t2_3x3    = t2_3x3
        self.t2_tr_3x3 = t2_3x3.transpose()
        
    def __compute_mbr_bnd_loc_3x3(self):
        self.mbr_loc_3x3 = math_utils.matmul_abc(self.t1_3x3, 
                                                 self.mbr_glb_3x3, 
                                                 self.t1_tr_3x3)
        
        self.bnd_loc_3x3 = math_utils.matmul_abc(self.t1_3x3, 
                                                 self.bnd_glb_3x3, 
                                                 self.t1_tr_3x3)

    def __compute_b_mbr_3x6(self):
        ix = self.i_loc[0]
        iy = self.i_loc[1]

        jx = self.j_loc[0]
        jy = self.j_loc[1]

        kx = self.k_loc[0]
        ky = self.k_loc[1]
        
        bi = jy - ky
        bj = ky - iy
        bk = iy - jy

        ci = kx - jx
        cj = ix - kx
        ck = jx - ix

        b_mbr_3x6 = np.array([
            [bi, 0,  bj, 0,  bk, 0 ],
            [0,  ci, 0,  cj, 0,  ck],
            [ci, bi, cj, bj, ck, bk]],
            dtype=float)
        
        coeff = 1 / (2 * self.area)
        b_mbr_3x6 *= coeff
        self.b_mbr_3x6 = b_mbr_3x6

    def __compute_k_mbr_6x6(self):
        b_mbr_tr_6x3 = self.b_mbr_3x6.transpose()

        k_mbr_6x6 = math_utils.matmul_abc(b_mbr_tr_6x3,
                                          self.mbr_loc_3x3,
                                          self.b_mbr_3x6)

        k_mbr_6x6 *= self.area
        self.k_mbr_6x6 = k_mbr_6x6
