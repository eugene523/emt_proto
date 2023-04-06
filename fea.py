from enum import Enum
import numpy as np


DOF = 6
'''Degrees of freedom.'''


class Constraint(Enum):
    Free = 0
    Fixed = 1


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


class Elem:
    def __init__(self, i: Node, j: Node, k: Node):
        self.i = i
        self.j = j
        self.k = k
        self.material_tag = -1

    def area(self) -> float:
        pass

