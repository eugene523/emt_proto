from enum import Enum


class NodeGroup(Enum):
    LEFT   = 0
    RIGHT  = 1
    TOP    = 2
    BOTTOM = 3
    N00    = 4
    N01    = 5
    N10    = 6
    N11    = 7


class Panel:
    def __init__(self, length: float, width: float):
        self.length = length
        self.width = width
