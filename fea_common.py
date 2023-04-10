from enum import Enum

DOF = 6
'''Degrees of freedom.'''


class Constraint(Enum):
    Free = 0
    Fixed = 1