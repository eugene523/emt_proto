import orth2d
import shellmat
from panel import Panel, NodeGroup
from boundary import *


def demo_1x1():
    material = orth2d.get_orth2d_mock(shellmat.Orth2dMockKind.KMU4)

    shell_material = shellmat.ShellMaterial()
    shell_material.add_ply(material, 1e-3, 0)
    shell_material.compute()

    p = Panel(length=1.0, width=1.0)
    p.material = shell_material
    p.elem_length = 0.1
    p.do_mesh()
    p.set_constraint(NodeGroup.LFT, ConstraintVector.new_fixed())
    p.set_force(NodeGroup.RGT, ForceVector.new(1, 0, 0, 0, 0, 0))
    p.compute()
    p.show()