import math_utils
from fea_common import *
import numpy as np


class Node:
    def __init__(self, index: int, x: float, y: float, z: float):
        self.index = index
        '''Index of the node in a mesh.'''

        self.x = x
        self.y = y
        self.z = z
        self.constraints = [Constraint.Free] * DOF
        self.force: np.ndarray = np.zeros((1, DOF))

    def append_constraints(self, constraints: list[Constraint]):
        assert len(constraints) == DOF
        for i in range(0, DOF):
            self.constraints[i] |= constraints[i]

    def set_force(self, force: np.ndarray):
        assert force.shape == (1, DOF)
        self.force = force

    def get_coord_vect(self) -> np.ndarray:
        return np.array([self.x, self.y, self.z], dtype=float)


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

    def get_nnodes(self) -> int:
        return len(self.nodes)

    def get_nelements(self) -> int:
        return len(self.elements)
    

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

    def __create_nodes(self, mesh: Mesh, nlen: int, nwid: int):
        a = self.a
        b = self.b
        d = self.d

        dx1 = (d.x - a.x) / nlen
        dy1 = (d.y - a.y) / nlen
        dz1 = (d.z - a.z) / nlen

        dx2 = (b.x - a.x) / nwid
        dy2 = (b.y - a.y) / nwid
        dz2 = (b.z - a.z) / nwid

        for j in range(0, nwid + 1):
            for i in range(0, nlen + 1):
                n = (nlen + 1) * j + i
                x = a.x + (i * dx1 + j * dx2)
                y = a.y + (i * dy1 + j * dy2)
                z = a.z + (i * dz1 + j * dz2)
                node = Node(n, x, y, z)
                mesh.add_node(node)
    
    def __create_elements_tria(self, mesh: Mesh, nlen: int, nwid: int, start_variant: int):
        assert start_variant == 1 or start_variant == -1
        for j in range(0, nwid):
            variant = start_variant
            for i in range(0, nlen):

                # n01 --- n11
                #  |       |
                #  |       |
                # n00 --- n10

                n00 = j * (nlen + 1) + i
                n10 = n00 + 1
                n01 = (j + 1) * (nlen + 1) + i
                n11 = n01 + 1

                if variant == 1:
                    mesh.add_element(n00, n10, n11)
                    mesh.add_element(n11, n01, n00)
                else:
                    mesh.add_element(n00, n10, n01)
                    mesh.add_element(n11, n01, n10)
                
                variant *= -1
            start_variant *= -1
    
    def mesh_tria(self, nlen: int, nwid: int, start_variant: int):
        m = Mesh()
        self.__create_nodes(m, nlen, nwid)
        self.__create_elements_tria(m, nlen, nwid, start_variant)
        return m


