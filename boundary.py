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
        return c

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

    def superposition(self, other: 'ConstraintVector'):
        if other.tx == Constraint.FIXED:
            self.tx = Constraint.FIXED

        if other.ty == Constraint.FIXED:
            self.ty = Constraint.FIXED

        if other.tz == Constraint.FIXED:
            self.tz = Constraint.FIXED

        if other.rx == Constraint.FIXED:
            self.rx = Constraint.FIXED

        if other.ry == Constraint.FIXED:
            self.ry = Constraint.FIXED

        if other.rz == Constraint.FIXED:
            self.rz = Constraint.FIXED

    def is_free(self) -> bool:
        return (self.tx == Constraint.FREE and
                self.ty == Constraint.FREE and
                self.tz == Constraint.FREE and
                self.rx == Constraint.FREE and
                self.ry == Constraint.FREE and
                self.rz == Constraint.FREE)
    

class ForceVector:
    def __init__(self):
        self.fx: float = 0.0
        self.fy: float = 0.0
        self.fz: float = 0.0
        self.mx: float = 0.0
        self.my: float = 0.0
        self.mz: float = 0.0

    @classmethod
    def new(cls,
            fx: float, fy: float, fz: float, 
            mx: float, my: float, mz: float) -> 'ForceVector':
        f = cls()
        f.fx: float = fx
        f.fy: float = fy
        f.fz: float = fz
        f.mx: float = mx
        f.my: float = my
        f.mz: float = mz
        return f

    def superposition(self, other: 'ForceVector'):
        self.fx += other.fx
        self.fy += other.fy
        self.fz += other.fz
        self.mx += other.mx
        self.my += other.my
        self.mz += other.mz