from enum import Enum
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


class Material:
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
            raise Exception("Параметр E1 материала не задан.")
        
        if self.e2 == 0.0:
            raise Exception("Параметр E2 материала не задан.")
        
        if self.g12 == 0.0:
            self.g12 = self.e1 / (2 * (1 + self.nu12))

        self.nu21 = self.e2 * self.nu12 / self.e1

        self.q12.fill(0)
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
        fos = 1.0 / fi
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
        
        fxx = 1.0 / sig1t_sq if sig1 >= 0.0 else 1.0 / sig1c_sq

        fyy = 1.0 / sig2t_sq if sig2 >= 0.0 else 1.0 / sig2c_sq

        fxy = -1.0 / sig1t_sq if sig1 * sig2 >= 0.0 else -1.0 / sig1c_sq

        fss = 1.0 / self.tau_max ** 2

        fi = (fxx * sig1 ** 2 +
              fyy * sig2 ** 2 +
              fxy * sig1 * sig2 +
              fss * tau12 ** 2)
        
        fos = 1.0 / fi ** 0.5

        mos = fos - 1

        return Criterion(CriterionType.HILL, fi, fos, mos)
        