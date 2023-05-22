import math
import math_utils
from boundary import *
import numpy as np


class Node:
    def __init__(self, index: int, x: float, y: float, z: float):
        self.index = index
        '''Index of the node in a mesh.'''

        self.x = x
        self.y = y
        self.z = z
        self.constraint_vector: ConstraintVector = None
        self.force_vector: np.ndarray = np.zeros((1, DOF))

    def append_constraints(self, constraints: list[Constraint]):
        assert len(constraints) == DOF
        for i in range(0, DOF):
            self.constraints[i] |= constraints[i]

    def set_force(self, force: np.ndarray):
        assert force.shape == (1, DOF)
        self.force = force

    def get_coord_vect(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)
    
    def is_equal(self, other: 'Node', eps: float) -> bool:
        if self is other:
            return True
        
        return (math.isclose(self.x, other.x, rel_tol=eps) and
                math.isclose(self.y, other.y, rel_tol=eps) and
                math.isclose(self.z, other.z, rel_tol=eps))

    def is_near_point(self, px: float, py: float, pz: float, eps: float) -> bool:
        dx = self.x - px
        dy = self.y - py
        dz = self.z - pz
        r = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        return r < eps
    
    def is_on_edge(self, 
                   x1: float, y1: float, z1: float,
                   x2: float, y2: float, z2: float,
                   eps: float) -> bool:
        # Компоненты вектора AB, образованного парой точек.
        dx = x2 - x1
        dy = y2 - y1
        dz = z2 - z1

        # Длина вектора AB.
        L = math.sqrt(dx ** 2 + dy ** 2 + dz ** 2)

        # Компоненты единичного вектора E_AB, сонаправленного с AB.
        dx /= L
        dy /= L
        dz /= L

        # Компоненты вектора AN, с началом в первой точке A и
        # с концом в узловой точке N.
        n_dx = self.x - x1
        n_dy = self.y - y1
        n_dz = self.z - z1

        # Скалярное произведение s = (E_AB x AN).
        # Важно! Скалярное произведение равно длине проекции AN на AB.
        s = n_dx * dx + n_dy * dy + n_dz * dz

        # Умножая скалярное произведение s (она же длина проекции)
        # на компоненты единичного вектора, получаем компоненты вектора
        # проекции AN на AB.
        dx *= s
        dy *= s
        dz *= s

        # Вычитаем из вектора AN вектор проекции, и получаем
        # вектор перпендикулярный прямой AB с концом в точке N.
        n_dx -= dx
        n_dy -= dy
        n_dz -= dz

        # Вычисляем длину перпендикулярного вектора.
        L = math.sqrt(n_dx ** 2 + n_dy ** 2 + n_dz ** 2)

        # Если длина L меньше delta, то считаем что точка лежит на прямой.
        return L < eps


class Elem:
    def __init__(self, i: Node, j: Node, k: Node):
        self.i: Node = i
        self.j: Node = j
        self.k: Node = k
        self.tag: int = 0

    def area(self) -> float:
        return math_utils.get_tria_area(self.i.x, self.i.y, self.i.z,
                                        self.j.x, self.j.y, self.j.z,
                                        self.k.x, self.k.y, self.k.z)
    
    def set_tag(self, tag: int):
        self.tag = tag


class Mesh:
    def __init__(self):
        self.nodes: list[Node] = []
        self.elements: list[Elem] = []

    def add_node(self, node: Node):
        self.nodes.append(node)

    def add_element(self, i_index: int, j_index: int, k_index: int):
        elem = Elem(self.nodes[i_index],
                    self.nodes[j_index],
                    self.nodes[k_index])
        
        self.elements.append(elem)

    def set_tag_to_all_elements(self, tag: int):
        for e in self.elements:
            e.set_tag(tag)

    def area(self) -> float:
        a = 0
        for elem in self.elements:
            a += elem.area()
        return a

    def get_n_nodes(self) -> int:
        return len(self.nodes)

    def get_n_elements(self) -> int:
        return len(self.elements)
    
    def select_nodes_near_point(self, 
                                px: float, py: float, pz: float, 
                                eps: float) -> list[Node]:
        selected = []
        for node in self.nodes:
            if node.is_near_point(px, py, pz, eps):
                selected.append(node)
        return selected
    
    def select_nodes_on_edge(self, 
                             x1: float, y1: float, z1: float,
                             x2: float, y2: float, z2: float,
                             eps: float) -> list[Node]:
        selected = []
        for node in self.nodes:
            if node.is_on_edge(x1, y1, z1, x2, y2, z2, eps):
                selected.append(node)
        return selected
    

# -----------------------------------------------------------------------------


class Point:
    def __init__(self, x: float, y: float, z: float):
        self.x: float = x
        self.y: float = y
        self.z: float = z

    @classmethod
    def new_from_another(cls, another: 'Point') -> 'Point':
        return cls(another.x, another.y, another.z)

    def translate(self, dx: float, dy: float, dz: float):
        return Point(self.x + dx, self.y + dy, self.z + dz)
    
    def translate_mod(self, dx: float, dy: float, dz: float):
        self.x += dx
        self.y += dy
        self.z += dz


class Quad:

    # b ----- c
    # |       |
    # a ----- d

    def __init__(self, a: Point, b: Point, c: Point, d: Point):
        self.a = a
        self.b = b
        self.c = c
        self.d = d

    @classmethod
    def new_by_translating_side(cls, a: Point, b: Point, 
                                dx: float, dy: float, dz: float) -> 'Quad':
        _a = Point.new_from_another(a)
        _b = Point.new_from_another(b)
        c = b.translate(dx, dy, dz)
        d = a.translate(dx, dy, dz)
        return cls(_a, _b, c, d)
    
    @classmethod
    def new_by_coord(cls,
                     ax: float, ay: float, az: float,
                     bx: float, by: float, bz: float,
                     cx: float, cy: float, cz: float,
                     dx: float, dy: float, dz: float) -> 'Quad':
        a = Point(ax, ay, az)
        b = Point(bx, by, bz)
        c = Point(cx, cy, cz)
        d = Point(dx, dy, dz)
        q = cls(a, b, c, d)
        return q
    
    def area(self) -> float:
        s1 = math_utils.get_tria_area(self.a.x, self.a.y, self.a.z,
                                      self.b.x, self.b.y, self.b.z,
                                      self.c.x, self.c.y, self.c.z)
        
        s2 = math_utils.get_tria_area(self.a.x, self.a.y, self.a.z,
                                      self.d.x, self.d.y, self.d.z,
                                      self.c.x, self.c.y, self.c.z)
        return s1 + s2
    
    def translate(self, dx: float, dy: float, dz: float) -> 'Quad':
        a = self.a.translate(dx, dy, dz)
        b = self.b.translate(dx, dy, dz)
        c = self.c.translate(dx, dy, dz)
        d = self.d.translate(dx, dy, dz)
        return Quad(a, b, c, d)
    
    def translate_mod(self, dx: float, dy: float, dz: float):
        self.a.translate_mod(dx, dy, dz)
        self.b.translate_mod(dx, dy, dz)
        self.b.translate_mod(dx, dy, dz)
        self.d.translate_mod(dx, dy, dz)

    def __create_nodes(self, mesh: Mesh, n_len: int, n_wid: int):
        a = self.a
        b = self.b
        d = self.d

        dx1 = (d.x - a.x) / n_len
        dy1 = (d.y - a.y) / n_len
        dz1 = (d.z - a.z) / n_len

        dx2 = (b.x - a.x) / n_wid
        dy2 = (b.y - a.y) / n_wid
        dz2 = (b.z - a.z) / n_wid

        for j in range(0, n_wid + 1):
            for i in range(0, n_len + 1):
                n = (n_len + 1) * j + i
                x = a.x + (i * dx1 + j * dx2)
                y = a.y + (i * dy1 + j * dy2)
                z = a.z + (i * dz1 + j * dz2)
                node = Node(n, x, y, z)
                mesh.add_node(node)
    
    def __create_elements_tria(self, mesh: Mesh, n_len: int, n_wid: int, start_variant: int):
        assert start_variant == 1 or start_variant == -1
        for j in range(0, n_wid):
            variant = start_variant
            for i in range(0, n_len):

                # n01 --- n11
                #  |       |
                #  |       |
                # n00 --- n10

                n00 = j * (n_len + 1) + i
                n10 = n00 + 1
                n01 = (j + 1) * (n_len + 1) + i
                n11 = n01 + 1

                if variant == 1:
                    mesh.add_element(n00, n10, n11)
                    mesh.add_element(n11, n01, n00)
                else:
                    mesh.add_element(n00, n10, n01)
                    mesh.add_element(n11, n01, n10)
                
                variant *= -1
            start_variant *= -1
    
    def mesh_tria(self, n_len: int, n_wid: int, start_variant: int):
        m = Mesh()
        self.__create_nodes(m, n_len, n_wid)
        self.__create_elements_tria(m, n_len, n_wid, start_variant)
        return m
