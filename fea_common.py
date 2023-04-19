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
        self.tx = Constraint.FREE
        self.ty = Constraint.FREE
        self.tz = Constraint.FREE
        self.rx = Constraint.FREE
        self.ry = Constraint.FREE
        self.rz = Constraint.FREE

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