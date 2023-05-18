from enum import Enum
from math import sqrt
import math_utils
import numpy as np


class CriterionType(Enum):
    MAX_STRESS = 0
    HILL = 1
    TSAI_WU = 2
    HOFFMAN = 3


class CriterionValueType(Enum):
    FAILURE_INDEX = 0
    FACTOR_OF_SAFETY = 1
    MARGIN_OF_SAFETY = 2


class Criterion:
    def __init__(self, 
                 criterion_type: CriterionType, 
                 failure_index: float, 
                 factor_of_safety: float, 
                 margin_of_safety: float):
        self.criterion_type   = criterion_type
        self.failure_index    = failure_index
        self.factor_of_safety = factor_of_safety
        self.margin_of_safety = margin_of_safety

    def get_criterion_value(self, criterion_value_type: CriterionValueType) -> float:
        match criterion_value_type:
            case CriterionValueType.FAILURE_INDEX:
                return self.failure_index
            
            case CriterionValueType.FACTOR_OF_SAFETY:
                return self.factor_of_safety
            
            case CriterionValueType.MARGIN_OF_SAFETY:
                return self.margin_of_safety
            
            case _:
                raise Exception(f"Unknown criterion value type {criterion_value_type}")


class Orth2d:
    '''Orthotropic prepreg material'''

    '''
    Принятые обозначения:
    (12) - Главная система координат:
    1 - направление вдоль волокон.
    2 - направление поперек волокон.
    T (Tension)     - растяжение материала.
    C (Compression) - сжатие материала.
    '''

    def __init__(self):
        self.name = "No name"

        self.e1   = 0.0   # Модуль упругости вдоль волокон.
        self.e2   = 0.0   # Модуль упругости поперек волокон.
        self.g12  = 0.0   # Модуль сдвига.
        self.nu12 = 0.0   # Коэффициент Пуассона при растяжении вдоль волокон.
        self.nu21 = 0.0   # Коэффициент Пуассона при растяжении поперек волокон.

        self.sig1t   = 0.0   # Предел прочности на растяжение вдоль волокон.
        self.sig1c   = 0.0   # Предел прочности на сжатие вдоль волокон.
        self.sig2t   = 0.0   # Предел прочности на расяжение поперек волокон.
        self.sig2c   = 0.0   # Предел прочности на сжатие поперек волокон.
        self.tau_max = 0.0   # Предел прочности на сдвиг.

        self.eps1t = 0.0      # Максимальные деформации при растяжении вдоль волокон.
        self.eps1c = 0.0      # Максимальные деформации при сжатии вдоль волокон.
        self.eps2t = 0.0      # Максимальные деформации при растяжении поперек волокон.
        self.eps2c = 0.0      # Максимальные деформации при сжатии поперек волокон.
        self.gamma_max = 0.0  # Максимальные сдвиговые деформации.

        self.density = 0.0;  # Плотность материала.
        self.prep_h  = 0.0;  # Толщина препрега.
        
        # Матрица упругости (мембранная)
        self.q12 = np.zeros((3, 3), dtype=float)

    def compute(self):
        '''Подготовить перед использованием в прочностных расчетах'''

        if self.e1 == 0.0:
            raise Exception("Parameter e1 of the material is not setted.")
        
        if self.e2 == 0.0:
            raise Exception("Parameter e2 of the material is not setted.")
        
        if self.nu12 == 0.0:
            raise Exception("Parameter nu12 of material is not setted.")
        
        if self.g12 == 0.0:
            self.g12 = self.e1 / (2 * (1 + self.nu12))

        self.nu21 = self.e2 * self.nu12 / self.e1

        self.q12.fill(0.0)
        k = 1 - self.nu12 * self.nu21
        self.q12[0, 0] = self.e1 / k
        self.q12[0, 1] = self.e1 * self.nu21 / k
        self.q12[1, 0] = self.e2 * self.nu12 / k
        self.q12[1, 1] = self.e2 / k
        self.q12[2, 2] = self.g12

    def get_sig12(self, eps12: np.ndarray) -> np.ndarray:
        return self.q12.dot(eps12)
    
    def get_criterion_max_stress(self, sig12: np.ndarray) -> Criterion:
        sig1  = sig12[0]
        sig2  = sig12[1]
        tau12 = sig12[2]
        
        v1 = abs(sig1 / self.sig1t) if sig1 >= 0.0 else abs(sig1 / self.sig1c)
        v2 = abs(sig2 / self.sig2t) if sig2 >= 0.0 else abs(sig2 / self.sig2c)
        v3 = abs(tau12 / self.tau_max)

        fi  = max(v1, v2, v3)
        fos = 1 / fi
        mos = fos - 1

        return Criterion(CriterionType.MAX_STRESS, fi, fos, mos)
    
    def get_criterion_hill(self, sig12: np.ndarray) -> Criterion:
        sig1  = sig12[0]
        sig2  = sig12[1]
        tau12 = sig12[2]

        sig1t_sq = self.sig1t ** 2
        sig1c_sq = self.sig1c ** 2
        sig2t_sq = self.sig2t ** 2
        sig2c_sq = self.sig2c ** 2
        
        fxx = (1 / sig1t_sq) if sig1 >= 0.0 else (1 / sig1c_sq)

        fyy = (1 / sig2t_sq) if sig2 >= 0.0 else (1 / sig2c_sq)

        fxy = (-1 / sig1t_sq) if sig1 * sig2 >= 0.0 else (-1 / sig1c_sq)

        fss = 1 / self.tau_max ** 2

        fi = (fxx * sig1 ** 2 +
              fyy * sig2 ** 2 +
              fxy * sig1 * sig2 +
              fss * tau12 ** 2)
        
        fos = 1 / sqrt(fi)

        mos = fos - 1

        return Criterion(CriterionType.HILL, fi, fos, mos)
    
    def get_criterion_tsai_wu(self, sig12: np.ndarray) -> Criterion:
        sig1  = sig12[0]
        sig2  = sig12[1]
        tau12 = sig12[2]

        fx  = (1 / self.sig1t) - (1 / self.sig1c)
        fy  = (1 / self.sig2t) - (1 / self.sig2c)
        fxx = 1 / (self.sig1t * self.sig1c)
        fyy = 1 / (self.sig2t * self.sig2c)
        fss = 1 / self.tau_max ** 2

        # coefficient of interaction
        ixy = 1

        fxy = -(ixy * sqrt(fxx * fyy))

        a = (fxx * (sig1 ** 2) +
             fyy * (sig2 ** 2) +
             fxy * (sig1 * sig2) +
             fss * (tau12 ** 2))
        
        b = fx * sig1 + fy * sig2

        fi = a + b

        d = sqrt(b ** 2 + 4 * a)
        fos = (-b + d) / (2 * a)

        mos = fos - 1

        return Criterion(CriterionType.TSAI_WU, fi, fos, mos)
    
    def get_criterion_hoffman(self, sig12: np.ndarray) -> Criterion:
        sig1  = sig12[0]
        sig2  = sig12[1]
        tau12 = sig12[2]

        fx  = (1 / self.sig1t) - (1 / self.sig1c)
        fy  = (1 / self.sig2t) - (1 / self.sig2c)
        fxx =  1 / (self.sig1t * self.sig1c)
        fyy =  1 / (self.sig2t * self.sig2c)
        fxy = -1 / (self.sig1t * self.sig1c)
        fss =  1 / (self.tau_max ** 2)

        a = (fxx * (sig1 ** 2) +
             fyy * (sig2 ** 2) +
             fxy * (sig1 * sig2) +
             fss * (tau12 ** 2))
        
        b = fx * sig1 + fy * sig2

        fi = a + b

        d = sqrt(b ** 2 + 4 * a)
        fos = (-b + d) / (2 * a)

        mos = fos - 1

        return Criterion(CriterionType.HOFFMAN, fi, fos, mos)
    
    def get_criteria(self, sig12: np.ndarray) -> list[Criterion]:
        criteria = [
            self.get_criterion_max_stress(sig12),
            self.get_criterion_hill(sig12),
            self.get_criterion_tsai_wu(sig12),
            self.get_criterion_hoffman(sig12)
        ]
        return criteria
    
    def is_isotropic(self) -> bool:
        math_utils.nearly_equal(self.e1, self.e2, 1e-3)


class Orth2dMockKind(Enum):
    KMU4 = 1
    VKU25 = 2

def get_orth2d_mock(material_kind: Orth2dMockKind) -> Orth2d:
    match material_kind:
        case Orth2dMockKind.KMU4:
            m = Orth2d()
            m.name    = "KMU4"
            m.e1      = 1.28e+11
            m.e2      = 8.4e+9
            m.g12     = 4.6e+9
            m.nu12    = 0.36
            m.sig1t   = 8.2e+8
            m.sig1c   = 1.0e+9
            m.sig2t   = 4.8e+7
            m.sig2c   = 1.5e+8
            m.tau_max = 6.2e+7
            m.density = 1.78e+3
            m.prep_h  = 2.3e-4
            m.compute()
            return m
        
        case Orth2dMockKind.VKU25:
            m = Orth2d()
            m.name    = "VKU25"
            m.e1      = 1.1e+11
            m.e2      = 8.21e+9
            m.g12     = 4.08e+9
            m.nu12    = 0.25
            m.sig1t   = 2.07e+9
            m.sig1c   = 1.14e+9
            m.sig2t   = 4.1e+7
            m.sig2c   = 1.63e+8
            m.tau_max = 9.0e+7
            m.density = 1.57e+3
            m.prep_h  = 2.15e-4
            m.compute()
            return m
        
        case _:
            raise Exception(f"Unknown Orth2dMockKind ({material_kind}).")