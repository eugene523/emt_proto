from math import sqrt
import numpy as np


def nearly_equal(a: float, b: float, eps: float) -> bool:
    if a == b:
        return True
    
    a_abs = abs(a)
    b_abs = abs(b)

    if (a > 0 and b > 0) or (a < 0 and b < 0):
        min_val = min(a_abs, b_abs)
        max_val = max(a_abs, b_abs)

        ab_eps = (max_val - min_val) / min_val
        if ab_eps < eps:
            return True
    
    half_eps = eps / 2
    return a_abs < half_eps and b_abs < half_eps


def get_tria_area(ax: float, ay: float, az: float,
                  bx: float, by: float, bz: float,
                  cx: float, cy: float, cz: float):

	# Площадь треугольника ABC вычисляем как половину от 
	# модуля вектора U, который является векторным произведением U = AB x BC.
	# Векторное произведение вычисляется как определитель:
	#
	#     | i  j  k  |
	# U = | x1 y1 z1 |
	#     | x2 y2 z2 |
	#
	# i,  j,  k  - единичные векторы базиса (в коде не используються)
	# x1, y1, z1 - компоненты вектора AB
	# x2, y2, z2 - компоненты вектора AC

	# Вычисляем компоненты вектора AB
	x1 = bx - ax
	y1 = by - ay
	z1 = bz - az

	# Вычисляем компоненты вектора AC
	x2 = cx - ax
	y2 = cy - ay
	z2 = cz - az

	# Вычисляем компоненты вектора U
	u1 =   y1 * z2 - y2 * z1
	u2 = -(x1 * z2 - x2 * z1)
	u3 =   x1 * y2 - x2 * y1

	# Вычисляем модуль вектора U
	u = sqrt(u1 * u1 + u2 * u2 + u3 * u3)

	# Половина от модуля вектора U - есть площадь треугольника ABC
	return 0.5 * u


def extract_submatrix(matrix: np.ndarray, 
                      row_from: int, col_from:int, 
                      nrows: int, ncols: int) -> np.ndarray:
    extracted = np.ndarray((nrows, ncols), dtype=float)
    for r in range(0, nrows):
        for c in range(0, ncols):
            val = matrix[row_from + r, col_from + c]
            extracted[r, c] = val
    return extracted


def matmul_abc(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    ab = np.matmul(a, b)
    abc = np.matmul(ab, c)
    return abc


class CooBuilder:
    def __init__(self, max_arr_size: int, n_indeces: int, block_size: int):
        self.rows: np.ndarray = np.zeros((1, max_arr_size), dtype=int)
        self.cols: np.ndarray = np.zeros((1, max_arr_size), dtype=int)
        self.vals: np.ndarray = np.zeros((1, max_arr_size), dtype=float)
        self.place: int = 0
        self.n_indeces: int = n_indeces
        self.block_size: int = block_size

    def accept_matrix(self, m: np.ndarray, indeces: list[int]):
        assert(len(indeces) == self.n_indeces)
        assert((m.shape[0] / self.n_indeces) == self.block_size)
        for a in range(self.n_indeces):
            _r = indeces[a]
            for b in range(self.n_indeces):
                _c = indeces[b]
                for r in range(self.block_size):
                    for c in range(self.block_size):
                        i = a * self.block_size + r
                        j = b * self.block_size + c
                        v = m[i, j]
                        _i = _r * self.block_size + r
                        _j = _c * self.block_size + c
                        self.rows[self.place] = _i
                        self.cols[self.place] = _j
                        self.vals[self.place] = v
                        self.place += 1

    def set_row_zero(self, row: int):
        for i in range(self.place):
            if self.rows[i] == row:
                self.vals[i] == 0.0

    def set_col_zero(self, col: int):
        for i in range(self.place):
            if self.cols[i] == col:
                self.vals[i] == 0.0

    def set_val(self, row: int, col: int, val: float):
        for i in range(self.place):
            if self.rows[i] == row and self.cols[i] == col:
                self.vals[i] = val
                return
        raise Exception(f"No such element with row={row}, col={col}.")

    def set_row_col_to_zero_and_place_1(self, rc: int):
        self.set_row_zero(rc)
        self.set_col_zero(rc)
        self.set_val(rc, rc, 1.0)