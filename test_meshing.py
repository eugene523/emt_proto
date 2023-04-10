from meshing import Mesh, Point, Quad
import unittest


class TestGeom(unittest.TestCase):
    def test_quad(self):
        dx = 2.5
        dy = 1.5
        test_area = dx * dy

        a = Point(0,  0,  0)
        b = Point(0,  dy, 0)
        c = Point(dx, dy, 0)
        d = Point(dx, 0,  0)

        q1 = Quad(a, b, c, d)
        self.assertAlmostEqual(q1.area(), test_area)

        q2 = Quad.new_by_translating_side(a, b, dx, 0, 0)
        self.assertAlmostEqual(q2.area(), test_area)

        q3 = Quad.new_by_coord(0,  0,  0,
                               0,  dy, 0,
                               dx, dy, 0,
                               dx, 0,  0)
        
        self.assertAlmostEqual(q3.area(), test_area)

    def test_quad_translate(self):
        # If we translate geometric figure, then its area shouldn't change.
        dx = 2.5
        dy = 1.5
        test_area = dx * dy

        tx = 0.1
        ty = 2.7
        tz = 3.5

        q = Quad.new_by_coord(0,  0,  0,
                              0,  dy, 0,
                              dx, dy, 0,
                              dx, 0,  0)
        
        qt = q.translate(tx, ty, tz)
        
        self.assertAlmostEqual(qt.area(), test_area)
        self.assertAlmostEqual(q.area(), qt.area())

    def test_quad_mesh_tria(self):
        dx = 2.5
        dy = 1.5
        test_area = dx * dy

        q = Quad.new_by_coord(0,  0,  0,
                              0,  dy, 0,
                              dx, dy, 0,
                              dx, 0,  0)
        
        nlen = 10
        nwid = 20
        m = q.mesh_tria(nlen, nwid, 1)
        self.assertEqual(m.get_nnodes(), (nlen + 1) * (nwid + 1))
        self.assertEqual(m.get_nelements(), 2 * nlen * nwid)
        self.assertAlmostEqual(m.area(), test_area)


'''
class TestMesh(unittest.TestCase):
    def test_quad_mesh(self):
        dx = 2.5
        dy = 1.5
        q = Quad.new_by_coord(0,  0,  0,
                              0,  dy, 0,
                              dx, dy, 0,
                              dx, 0,  0)
'''

if __name__ == '__main__':
    unittest.main()