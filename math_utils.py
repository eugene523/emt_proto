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


def extract_submatrix(matrix: np.ndarray, 
                      row_from: int, col_from:int, 
                      nrows: int, ncols: int) -> np.ndarray:
    extracted = np.ndarray((nrows, ncols), dtype=float)
    for r in range(0, nrows):
        for c in range(0, ncols):
            val = matrix[row_from + r, col_from + c]
            extracted[r, c] = val
    return extracted


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