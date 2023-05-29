import orth2d
import shellmat
import material_mock
from panel import Panel, NodeGroup
from boundary import *


def demo_1x1():
    material = material_mock.get_material_mock(material_mock.MaterialMockKind.D16)

    shell_material = shellmat.ShellMaterial()
    shell_material.add_ply(material, 1e-3, 0)
    shell_material.compute()

    p = Panel(length=1.0, width=1.0)
    p.material = shell_material
    p.elem_length = 0.1
    p.do_mesh()
    p.set_constraint(NodeGroup.LFT, ConstraintVector.new_fixed())
    p.set_force(NodeGroup.RGT, ForceVector.new(0, 1e+5, 0, 0, 0, 0))
    p.compute()
    p.show_static_deform_mesh()
    print(p.disp_mag)

if __name__ == '__main__':
    demo_1x1()