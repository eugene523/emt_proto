import math
from orth2d import *
from validation import *
import numpy as np


class PlyStress:
    def __init__(self, eps12: np.ndarray, sig12: np.ndarray, criteria: list[Criterion]):
        self.eps12 = eps12
        self.sig12 = sig12
        self.criteria = criteria


class Ply:
    def __init__(self, material: Orth2d, thickness: float, angle_radian: float):
        self.material = material
        self.thickness = thickness
        self.angle_radian = angle_radian
        self.nlayers = 0
        self.ztop = 0.0
        self.zbot = 0.0

        # Матрицы трансформации для перевода напряжений и деформаций из
        # глобальной системы координат (XY) в собственную (12).
        self.t1 = np.zeros((3, 3), dtype=float)
        self.t2 = np.zeros((3, 3), dtype=float)

    @classmethod
    def new_with_nlayers(cls, material: Orth2d, nlayers: int, angle_radian: float):
        if material.prep_h == 0.0:
            raise Exception("Can't create ply with zero prepreg thickness.")
        thickness = nlayers * material.prep_h
        ply = cls(material, thickness, angle_radian)
        ply.nlayers = nlayers
        return ply

    @classmethod
    def new_from_another(cls, another_ply: 'Ply'):
        ply = cls(another_ply.material, another_ply.thickness, another_ply.angle_radian)
        ply.nlayers = another_ply.nlayers
        return ply

    def compute(self):
        self.compute_t1_and_t2()
        self.compute_stiffness_matrix()

    def compute_t1_and_t2(self):
        c = math.cos(self.angle_radian)
        s = math.sin(self.angle_radian)

        c_sq = c ** 2
        s_sq = s ** 2
        sc   = s * c

        # -----------------------------

        self.t1.fill(0.0)

        self.t1[0, 0] = c_sq
        self.t1[0, 1] = s_sq
        self.t1[0, 2] = -2 * sc

        self.t1[1, 0] = s_sq
        self.t1[1, 1] = c_sq
        self.t1[1, 2] = 2 * sc

        self.t1[2, 0] = sc
        self.t1[2, 1] = -sc
        self.t1[2, 2] = c_sq - s_sq

        # -----------------------------

        self.t2.fill(0.0)

        self.t2[0, 0] = c_sq
        self.t2[0, 1] = s_sq
        self.t2[0, 2] = -sc

        self.t2[1, 0] = s_sq
        self.t2[1, 1] = c_sq
        self.t2[1, 2] = sc

        self.t2[2, 0] = 2 * sc
        self.t2[2, 1] = -2 * sc
        self.t2[2, 2] = c_sq - s_sq

    def compute_stiffness_matrix(self):
        qxy = np.matmul(self.t1, self.material.q12)
        qxy = np.matmul(qxy, self.t1.transpose())

        coeff = self.ztop - self.zbot
        a = qxy * coeff

        coeff = (self.ztop ** 3 - self.zbot ** 3) / 3
        b = qxy * coeff

        coeff = (self.ztop ** 2 - self.zbot ** 2) / 2
        c = qxy * coeff

        self.dxy = np.block([[a, c], 
                             [c, b]])
        
    def get_ply_stress_data(self, eps_shellmat_xy: np.ndarray):
        # Вектор мембранных деформаций пакета в координатах XY
        eps_membr_xy = eps_shellmat_xy[0:3]

        # Вектор изгибных деформаций пакета в координатах XY
        eps_bend_xy  = eps_shellmat_xy[3:6]

        # Матрица трансформации (XY) -> (12)
        t1_t = self.t1.transpose()

        # Координата z срединной поверхности
        h_mid = (self.zbot + self.ztop) / 2

        # Вектор мембранных деформаций в координатах XY
        eps_xy = eps_membr_xy + eps_bend_xy * h_mid

        # Вектор мембранных деформаций на серединной поверхности слоя в координатах 12
        eps_12 = np.matmul(t1_t, eps_xy)

        # Вектор напряжений на серединной поверхности в координатах 12
        sig_12 = self.material.get_sig12(eps_12)

        # Вычисляем критерии
        criteria = self.material.get_criteria(sig_12)

        return PlyStress(eps_12, sig_12, criteria)


class ShellMaterialStress:
    def __init__(self, ply_stress_list: list[PlyStress]):
        self.ply_stress_list = ply_stress_list


class ShellMaterial:
    def __init__(self):
        # Слои композиционного пакета.
        self.plies: list[Ply] = []

        # Толщина композиционного пакета.
        self.thickness = 0.0

        # Распределенный вес пакета.
        self.area_density = 0.0

        # Матрица [6x6] упругости композиционного пакета.
        self.dxy: np.ndarray = None

        # Подматрицы матрицы упругости композиционного пакета.
        self.b_mat: np.ndarray = None
        self.c_mat: np.ndarray = None
        self.d_mat: np.ndarray = None

        # Инвертированная матрица упругости - матрица податливости.
        self.dxy_inv: np.ndarray = None

        # Упругие технические постоянные
        self.ex = 0.0
        self.ey = 0.0
        self.gxy = 0.0
        self.nu_xy = 0.0
        self.nu_yx = 0.0

    def add_ply(self, material: Orth2d, ply_thickness: float, angle_degree: float):
        angle_radian = math.radians(angle_degree)
        ply = Ply(material, ply_thickness, angle_radian)
        self.plies.append(ply)

    def add_ply_nlayers(self, material: Orth2d, nlayers: int, angle_degree: float):
        angle_radian = math.radians(angle_degree)
        ply = Ply.new_with_nlayers(material, nlayers, angle_radian)
        self.plies.append(ply)

    def remove_ply(self, ply_index: int):
        self.plies.pop(ply_index)

    def nplies(self) -> int:
        return len(self.plies)

    def make_plies_symmetric(self):
        nplies = len(self.plies)
        for i in range(self.nplies() - 1, -1, -1):
            ply = self.plies[i]
            self.plies.append(Ply.new_from_another(ply))

    def repeat_plies(self, nrepetitions):
        nplies = self.nplies
        for _ in range(0, nrepetitions):
            for j in range(0, nplies):
                new_ply = Ply.new_from_another(self.plies[j])
                self.plies.append(new_ply)

    def validate(self) -> ValidationResult:
        v = ValidationResult()

        # Проверяем наличие слоев.
        if len(self.plies) == 0:
            v.add_issue("В композитном пакете отсутствуют слои.")
        
        # Проверяем у всех ли слоев задан материал.
        for i in range(0, len(self.plies)):
            if self.plies[i].material == None:
                v.add_issue(f"Отсутствует материал у слоя №{i + 1}.")

        # Проверка того, что все слои имеют ненулевую толщину.
        for i in range(0, len(self.plies)):
            if self.plies[i].thickness == 0.0:
                v.add_issue(f"Нулевые толщины слоев не допускаются. Слой №{i + 1}.")
        
        return v
    
    def compute(self):
        v = self.validate()
        if not v.is_ok():
            raise Exception(v)

    def __compute_materials(self):
        used_materials = set(self.plies)
        for m in used_materials:
            m.compute()

    def __compute_area_density(self):
        self.area_density = 0.0
        for ply in self.plies:
            self.area_density += ply.material.density * ply.thickness
    
    def __compute_thickness(self):
        self.thickness = 0.0
        for ply in self.plies:
            self.thickness += ply.thickness

    def __compute_matrices(self):
        zbot = -self.thickness / 2
        for i in range(len(self.plies) - 1, -1, -1):
            ztop = zbot + self.plies[i].thickness
            self.plies[i].ztop = ztop
            self.plies[i].zbot = zbot
            zbot = ztop

        for ply in self.plies:
            ply.compute()

        # Вычисляем матрицу упругости композитного пакета.
        self.dxy = np.zeros((6, 6), dtype=float)
        for ply in self.plies:
            self.dxy += ply.dxy
        
        self.dxy_inv = np.linalg.inv(self.dxy)

        self.b_mat = math_utils.extract_submatrix(self.dxy, 0, 0, 3, 3)
        self.c_mat = math_utils.extract_submatrix(self.dxy, 0, 3, 3, 3)
        self.d_mat = math_utils.extract_submatrix(self.dxy, 3, 3, 3, 3)
        
        

    