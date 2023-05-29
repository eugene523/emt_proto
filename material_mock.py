from enum import Enum
from orth2d import Orth2d

class MaterialMockKind(Enum):
    STEEL = 2
    D16   = 1
    KMU4  = 3
    VKU25 = 4

def get_material_mock(material_kind: MaterialMockKind) -> Orth2d:
    match material_kind:
        case MaterialMockKind.STEEL:
            m = Orth2d()
            m.name    = "STEEL"
            m.e1      = 2.1e+11
            m.e2      = 2.1e+11
            m.g12     = 8.1e+10
            m.nu12    = 3.0e-1
            m.nu21    = 3.0e-1
            m.sig1t   = 4.0e+8
            m.sig1c   = 4.0e+8
            m.sig2t   = 4.0e+8
            m.sig2c   = 4.0e+8
            m.tau_max = 2.3e+8
            m.density = 0.0 # not setted
            m.prep_h  = 0.0 # not setted
            m.compute()
            return m
        
        case MaterialMockKind.D16:
            m = Orth2d()
            m.name    = "D16"
            m.e1      = 7.2e+10
            m.e2      = 7.2e+10
            m.g12     = 2.76e+10
            m.nu12    = 3.0e-1
            m.nu21    = 3.0e-1
            m.sig1t   = 4.0e+8
            m.sig1c   = 4.0e+8
            m.sig2t   = 4.0e+8
            m.sig2c   = 4.0e+8
            m.tau_max = 2.3e+8
            m.density = 2.73e+3
            m.prep_h  = 0.0 # not setted
            m.compute()
            return m
        
        case MaterialMockKind.KMU4:
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
        
        case MaterialMockKind.VKU25:
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