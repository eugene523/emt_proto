import unittest
import orth2d
import shellmat

class TestShellMaterial(unittest.TestCase):
    def test_add_ply(self):
        kmu4 = orth2d.get_orth2d_mock(orth2d.Orth2dMockKind.KMU4)

        sm = shellmat.ShellMaterial()
        sm.add_ply(kmu4, 1.5,  0)
        sm.add_ply(kmu4, 1.5,  45)
        sm.add_ply(kmu4, 1.5, -45)
        sm.add_ply(kmu4, 1.5,  90)
        self.assertEqual(sm.nplies(), 4)

        sm.compute()
        self.assertAlmostEqual(sm.thickness, 6)

    def __test_ply(self, 
                   sm: shellmat.ShellMaterial, 
                   ply_index: int,
                   mat: orth2d.Orth2d, 
                   thickness: float,
                   angle_degree: float):
        self.assertIs(sm.get_ply_material(ply_index), mat)
        self.assertAlmostEqual(sm.get_ply_thickness(ply_index), thickness)
        self.assertAlmostEqual(sm.get_ply_angle_deg(ply_index), angle_degree)

    def test_make_symmetric(self):
        kmu4  = orth2d.get_orth2d_mock(orth2d.Orth2dMockKind.KMU4)
        vku25 = orth2d.get_orth2d_mock(orth2d.Orth2dMockKind.VKU25)

        # materials, ply thicknesses and agles
        m = [kmu4, vku25, kmu4, vku25]
        t = [1.0,  2.35,  3.0,  3.5]
        a = [0,    45,   -45,   90]

        sm = shellmat.ShellMaterial()
        sm.add_ply(m[0], t[0], a[0])
        sm.add_ply(m[1], t[1], a[1])
        sm.add_ply(m[2], t[2], a[2])
        sm.add_ply(m[3], t[3], a[3])
        sm.make_symmetric()
        self.assertEqual(sm.nplies(), 8)
        
        self.__test_ply(sm, 0, m[0], t[0], a[0])
        self.__test_ply(sm, 1, m[1], t[1], a[1])
        self.__test_ply(sm, 2, m[2], t[2], a[2])
        self.__test_ply(sm, 3, m[3], t[3], a[3])
        self.__test_ply(sm, 4, m[3], t[3], a[3])
        self.__test_ply(sm, 5, m[2], t[2], a[2])
        self.__test_ply(sm, 6, m[1], t[1], a[1])
        self.__test_ply(sm, 7, m[0], t[0], a[0])

    def test_repeat_plies(self):
        kmu4  = orth2d.get_orth2d_mock(orth2d.Orth2dMockKind.KMU4)
        vku25 = orth2d.get_orth2d_mock(orth2d.Orth2dMockKind.VKU25)

        # materials, ply thicknesses and agles
        m = [kmu4, vku25, kmu4]
        t = [1.0,  2.35,  3.0]
        a = [0,    45,   -45]

        sm = shellmat.ShellMaterial()
        sm.add_ply(m[0], t[0], a[0])
        sm.add_ply(m[1], t[1], a[1])
        sm.add_ply(m[2], t[2], a[2])
        sm.repeat_plies(3)
        self.assertEqual(sm.nplies(), 9)

        self.__test_ply(sm, 0, m[0], t[0], a[0])
        self.__test_ply(sm, 1, m[1], t[1], a[1])
        self.__test_ply(sm, 2, m[2], t[2], a[2])
        self.__test_ply(sm, 3, m[0], t[0], a[0])
        self.__test_ply(sm, 4, m[1], t[1], a[1])
        self.__test_ply(sm, 5, m[2], t[2], a[2])
        self.__test_ply(sm, 6, m[0], t[0], a[0])
        self.__test_ply(sm, 7, m[1], t[1], a[1])
        self.__test_ply(sm, 8, m[2], t[2], a[2])

    def test_compute(self):
        pass # todo




if __name__ == '__main__':
    unittest.main()