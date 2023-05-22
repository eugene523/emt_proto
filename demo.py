from panel import Panel, NodeGroup
from boundary import *


def demo_1x1():
    p = Panel(1, 1)
    p.elem_length = 0.1
    p.do_mesh()
    p.set_constraint(NodeGroup.LFT, ConstraintVector.new_fixed)
    p.set_force(NodeGroup.RGT, ForceVector(1, 0, 0, 0, 0, 0))
    p.show()