from enum import Enum
from shellmat import ShellMaterial
from meshing import *
from validation import *


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
    ERR_LENGTH_NOT_SET      = "Panel length is not setted."
    ERR_WIDTH_NOT_SET       = "Panel width is not setted."
    ERR_ELEM_LENGTH_NOT_SET = "Mesh element length is not setted."

    def __init__(self, length: float, width: float):
        self.length: float = length
        self.width: float = width
        self.material: ShellMaterial = None
        self.elem_length: float = 0.0
        self.mesh: Mesh = None

    def validate_before_meshing(self) -> ValidationResult:
        v = ValidationResult()

        if self.length == 0:
            v.add_issue(self.ERR_LENGTH_NOT_SET)

        if self.width == 0:
            v.add_issue(self.ERR_WIDTH_NOT_SET)

        if self.elem_length == 0:
            v.add_issue(self.ERR_ELEM_LENGTH_NOT_SET)

        return v
        
    def do_mesh(self):
        v = self.validate_before_meshing()
        if not v.is_ok():
            raise Exception(v)
        
        m = Mesh()
        




    