from panel import Panel


def demo_1x1():
    p = Panel(1, 1)
    p.elem_length = 0.1
    p.do_mesh()
    p.show()