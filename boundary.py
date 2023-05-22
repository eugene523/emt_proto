from enum import Enum

DOF = 6
'''Degrees of freedom.'''


class DofType(Enum):
    TX = 0
    TY = 1
    TZ = 2
    RX = 3
    RY = 4
    RZ = 5


class Constraint(Enum):
    FREE  = 0
    FIXED = 1


class ConstraintVector:
    def __init__(self):
        self.tx: Constraint = Constraint.FREE
        self.ty: Constraint = Constraint.FREE
        self.tz: Constraint = Constraint.FREE
        self.rx: Constraint = Constraint.FREE
        self.ry: Constraint = Constraint.FREE
        self.rz: Constraint = Constraint.FREE

    @classmethod
    def new_fixed(cls) -> 'ConstraintVector':
        c = cls()
        c.tx = Constraint.FIXED
        c.ty = Constraint.FIXED
        c.tz = Constraint.FIXED
        c.rx = Constraint.FIXED
        c.ry = Constraint.FIXED
        c.rz = Constraint.FIXED

    def set_dof(self, dof_type: DofType, constraint: Constraint):
        match dof_type:
            case DofType.TX:
                self.tx = constraint

            case DofType.TY:
                self.ty = constraint

            case DofType.TZ:
                self.tz = constraint

            case DofType.RX:
                self.rx = constraint

            case DofType.RY:
                self.ry = constraint

            case DofType.RZ:
                self.rz = constraint

    def is_totally_free(self) -> bool:
        return (self.tx == Constraint.Free and
                self.ty == Constraint.Free and
                self.tz == Constraint.Free and
                self.rx == Constraint.Free and
                self.ry == Constraint.Free and
                self.rz == Constraint.Free)
    

class ForceVector:
    def __init__(self, 
                 fx: float, fy: float, fz: float, 
                 mx: float, my: float, mz: float):
        self.fx: float = fx
        self.fy: float = fy
        self.fz: float = fz
        self.mx: float = mx
        self.my: float = my
        self.mz: float = mz