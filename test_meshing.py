from meshing import Mesh, Node, Point, Quad
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
        self.assertEqual(m.get_n_nodes(), (nlen + 1) * (nwid + 1))
        self.assertEqual(m.get_n_elements(), 2 * nlen * nwid)
        self.assertAlmostEqual(m.area(), test_area)


class TestMesh(unittest.TestCase):
    def test_node_is_equal(self):
        eps = 1e-6
        n1 = Node(0, 0.1, 0.2, 0.3)
        n2 = n1
        self.assertTrue(n1.is_equal(n2, eps))

        n2 = Node(0, 0.1 + 1e-7, 0.2 + 1e-7, 0.3 + 1e-7)
        self.assertTrue(n1.is_equal(n2, eps))

        n2 = Node(0, 0.1 + 1e-5, 0.2 + 1e-5, 0.3 + 1e-5)
        self.assertFalse(n1.is_equal(n2, eps))

    def test_node_is_on_edge(self):
        pass
        # Todo


if __name__ == '__main__':
    unittest.main()