from meshing import Node, Point, Quad
import unittest


def get_mock_quad(dx: float, dy: float) -> Quad:
    q = Quad.new_by_coord(0,  0,  0,
                          0,  dy, 0,
                          dx, dy, 0,
                          dx, 0,  0)
    return q


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

        q3 = get_mock_quad(dx, dy)
        
        self.assertAlmostEqual(q3.area(), test_area)

    def test_quad_translate(self):
        # If we translate geometric figure, then its area shouldn't change.
        dx = 2.5
        dy = 1.5
        test_area = dx * dy

        tx = 0.1
        ty = 2.7
        tz = 3.5

        q = get_mock_quad(dx, dy)  
        qt = q.translate(tx, ty, tz)

        self.assertAlmostEqual(qt.area(), test_area)
        self.assertAlmostEqual(q.area(), qt.area())

    def test_quad_mesh_tria(self):
        dx = 2.5
        dy = 1.5
        test_area = dx * dy

        q = get_mock_quad(dx, dy)
        nlen = 10
        nwid = 20
        m = q.mesh_tria(nlen, nwid, 1)
        self.assertEqual(m.get_n_nodes(), (nlen + 1) * (nwid + 1))
        self.assertEqual(m.get_n_elements(), 2 * nlen * nwid)
        self.assertAlmostEqual(m.area(), test_area)


class TestMesh(unittest.TestCase):
    def test_is_near_point(self):
        eps = 1e-6
        n1 = Node(0, 1, 1, 1)
        self.assertTrue(n1.is_near_point(1, 1, 1, eps))

        n2 = Node(0, 1 + 1e-5, 1, 1)
        self.assertFalse(n2.is_near_point(1, 1, 1, eps))

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
        eps = 1e-6
        n1 = Node(0, 0.5, 0.5, 0.5)
        self.assertTrue(n1.is_on_edge(0, 0, 0, 1, 1, 1, eps))
        
        n2 = Node(0, 0.5 + eps, 0.5 - eps, 0.5)
        self.assertFalse(n2.is_on_edge(0, 0, 0, 1, 1, 1, eps))

    def test_select_nodes(self):
        eps = 1e-6
        dx = 200
        dy = 100
        n_len = 20
        n_wid = 10
        n_nodes_dx_test = n_len + 1
        n_nodes_dy_test = n_wid + 1

        q = get_mock_quad(dx, dy)
        m = q.mesh_tria(n_len, n_wid, 1)

        # -----------------------------------------------------------

        n00 = m.select_nodes_near_point(0, 0, 0, eps)
        self.assertEqual(len(n00), 1)

        n10 = m.select_nodes_near_point(dx, 0, 0, eps)
        self.assertEqual(len(n10), 1)

        n01 = m.select_nodes_near_point(0, dy, 0, eps)
        self.assertEqual(len(n01), 1)

        n11 = m.select_nodes_near_point(dx, dy, 0, eps)
        self.assertEqual(len(n11), 1)

        # -----------------------------------------------------------

        bottom_nodes = m.select_nodes_on_edge(0, 0, 0,
                                              1, 0, 0,
                                              eps)
        self.assertEqual(len(bottom_nodes), n_nodes_dx_test)

        # -----------------------------------------------------------

        top_nodes = m.select_nodes_on_edge(0, dy, 0,
                                           1, dy, 0,
                                           eps)
        self.assertEqual(len(top_nodes), n_nodes_dx_test)

        # -----------------------------------------------------------

        left_nodes = m.select_nodes_on_edge(0, 0, 0,
                                            0, 1, 0,
                                            eps)
        self.assertEqual(len(left_nodes), n_nodes_dy_test)

        # -----------------------------------------------------------

        right_nodes = m.select_nodes_on_edge(dx, 0, 0,
                                             dx, 1, 0,
                                             eps)
        self.assertEqual(len(right_nodes), n_nodes_dy_test)


if __name__ == '__main__':
    unittest.main()