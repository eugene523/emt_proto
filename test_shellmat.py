import math
import unittest
import orth2d
import shellmat
import numpy as np

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
        tol = 1e-2
        kmu4 = orth2d.get_orth2d_mock(orth2d.Orth2dMockKind.KMU4)

        sm = shellmat.ShellMaterial()
        sm.add_ply_nlayers(kmu4, 1,  0.0)
        sm.add_ply_nlayers(kmu4, 2,  90.0)
        sm.add_ply_nlayers(kmu4, 3,  45.0)
        sm.add_ply_nlayers(kmu4, 4, -45.0)
        sm.compute()

        self.assertTrue(math.isclose(sm.thickness,    2.3e-3, rel_tol=tol))
        self.assertTrue(math.isclose(sm.area_density, 4.094,  rel_tol=tol))
        self.assertEqual(sm.nlayers(), 10)

        # Checking engineering constants
        self.assertTrue(math.isclose(sm.ex,    3.35e+10, rel_tol=tol))
        self.assertTrue(math.isclose(sm.ey,    4.29e+10, rel_tol=tol))
        self.assertTrue(math.isclose(sm.gxy,   2.44e+10, rel_tol=tol))
        self.assertTrue(math.isclose(sm.nu_xy, 4.15e-1,  rel_tol=tol))
        self.assertTrue(math.isclose(sm.nu_yx, 5.32e-1,  rel_tol=tol))

        # Checking stiffness matrix
        dxy_test = np.array([
            [ 9.88e+7,  5.25e+7, -6.94e+6,  1.09e+4, -1.57e+4,  2.15e+4],
            [ 5.25e+7,  1.27e+8, -6.94e+6, -1.57e+4,  2.05e+4,  2.15e+4],
            [-6.94e+6, -6.94e+6,  5.61e+7,  2.15e+4,  2.15e+4, -1.57e+4],
            [ 1.09e+4, -1.57e+4,  2.15e+4,  5.57e+1,  1.83e+1, -1.41e+1],
            [-1.57e+4,  2.05e+4,  2.15e+4,  1.83e+1,  5.33e+1, -1.41e+1],
            [ 2.15e+4,  2.15e+4, -1.57e+4, -1.41e+1, -1.41e+1,  1.99e+1]], 
            dtype=float)
        
        self.assertTrue(np.allclose(sm.dxy, dxy_test, rtol=tol))

        # -----------------------------------------------------------

        b_test = np.array([
            [ 9.88e+7,  5.25e+7, -6.94e+6],
            [ 5.25e+7,  1.27e+8, -6.94e+6],
            [-6.94e+6, -6.94e+6,  5.61e+7]], 
            dtype=float)
        
        self.assertTrue(np.allclose(sm.b_mat, b_test, rtol=tol))

        # -----------------------------------------------------------

        c_test = np.array([
            [ 1.09e+4, -1.57e+4,  2.15e+4],
            [-1.57e+4,  2.05e+4,  2.15e+4],
            [ 2.15e+4,  2.15e+4, -1.57e+4]],
            dtype=float)
        
        self.assertTrue(np.allclose(sm.c_mat, c_test, rtol=tol))
        
        # -----------------------------------------------------------

        d_test = np.array([
            [ 5.57e+1,  1.83e+1, -1.41e+1],
            [ 1.83e+1,  5.33e+1, -1.41e+1],
            [-1.41e+1, -1.41e+1,  1.99e+1]],
            dtype=float)
        
        self.assertTrue(np.allclose(sm.d_mat, d_test, rtol=tol))

    def test_compute2(self):
        tol = 1e-2
        kmu4 = orth2d.get_orth2d_mock(orth2d.Orth2dMockKind.KMU4)

        sm = shellmat.ShellMaterial()
        sm.add_ply_nlayers(kmu4, 2,  0.0)
        sm.add_ply_nlayers(kmu4, 3,  45.0)
        sm.add_ply_nlayers(kmu4, 3, -45.0)
        sm.add_ply_nlayers(kmu4, 2,  90.0)
        sm.make_symmetric()
        sm.compute()

        self.assertTrue(math.isclose(sm.thickness,    4.6e-3, rel_tol=tol))
        self.assertTrue(math.isclose(sm.area_density, 8.19,   rel_tol=tol))
        self.assertEqual(sm.nlayers(), 20)

        # Checking engineering constants
        self.assertTrue(math.isclose(sm.ex,    4.41e+10, rel_tol=tol))
        self.assertTrue(math.isclose(sm.ey,    4.41e+10, rel_tol=tol))
        self.assertTrue(math.isclose(sm.gxy,   2.16e+10, rel_tol=tol))
        self.assertTrue(math.isclose(sm.nu_xy, 3.86e-1,  rel_tol=tol))
        self.assertTrue(math.isclose(sm.nu_yx, 3.86e-1,  rel_tol=tol))

        # Checking stiffness matrix
        dxy_test = np.array([
            [2.38e+8,  9.2e+7,   0.0,      4.37e-11, 1.18e-11, 7.28e-12],
            [9.2e+7,   2.38e+8,  7.45e-9,  1.18e-11, 1.64e-11, 7.28e-12],
            [0.0,      7.45e-9,  9.92e+7,  7.28e-12, 3.64e-12, 1.0e-11 ],
            [4.37e-11, 1.18e-11, 7.28e-12, 6.77e+2,  1.4e+2,   6.6e+1  ],
            [1.18e-11, 1.64e-11, 7.28e-12, 1.4e+2,   2.08e+2,  6.6e+1  ],
            [7.28e-12, 3.64e-12, 1.0e-11,  6.6e+1,   6.6e+1,   1.53e+2 ]], 
            dtype=float)
        
        self.assertTrue(np.allclose(sm.dxy, dxy_test, rtol=tol))

        # Computing shell material stress
        dist_load = np.array([100, 200, 300, 400, 500, 600], dtype=float)
        sm_stress = sm.get_stress(dist_load)

        # Checking sig12 vectors
        sig12_table_test = np.array([
        #    sig1       sig2      tau12
            [ 4.17e+6,  2.37e+7,  3.19e+7],
            [ 4.48e+8, -2.04e+6,  9.44e+6],
            [-9.87e+7,  1.35e+7, -5.08e+6],
            [ 4.03e+7,  9.22e+5, -3.55e+6],
            [-4.01e+7, -9.15e+5,  3.53e+6],
            [ 9.84e+7, -1.35e+7,  5.08e+6],
            [-4.47e+8,  2.03e+6, -9.43e+6],
            [-4.13e+6, -2.37e+7, -3.19e+7]], 
            dtype=float)
        
        sig12_table = sm_stress.get_sig12_table()
        self.assertTrue(np.allclose(sig12_table, sig12_table_test, rtol=tol))

        # Checking criteria (failure index)
        fi_table_test = np.array([
        #   MaxStress Hill      TsaiWu
            [5.14e-1, 5.07e-1,  6.77e-1],
            [5.46e-1, 3.22e-1,  3.50e-1],
            [2.82e-1, 9.72e-2,  2.31e-1],
            [5.73e-2, 6.02e-3,  2.68e-2],
            [5.69e-2, 4.83e-3, -1.69e-2],
            [1.20e-1, 3.05e-2, -1.09e-1],
            [4.47e-1, 2.26e-1,  2.10e-1],
            [5.14e-1, 2.89e-1,  4.34e-3]],
            dtype=float)
        
        # Checking criteria (factor of safety)
        fos_table_test = np.array([
            [1.94,    1.4,     1.29],
            [1.83,    1.76,    1.77],
            [3.55,    3.21,    2.88],
            [1.74e+1, 1.29e+1, 1.22e+1],
            [1.76e+1, 1.44e+1, 1.68e+1],
            [8.33,    5.72,    5.66],
            [2.24,    2.1,     2.02],
            [1.95,    1.86,    2.28]],
            dtype=float)
        
        fos_table = sm_stress.get_crit_table(orth2d.CriterionValueType.FACTOR_OF_SAFETY)
        self.assertTrue(np.allclose(fos_table, fos_table_test, rtol=tol))


if __name__ == '__main__':
    unittest.main()