def nearly_equal(a: float, b: float, eps: float) -> bool:
    if a == b:
        return True
    
    a_abs = abs(a)
    b_abs = abs(b)

    if (a > 0 and b > 0) or (a < 0 and b < 0):
        min_val = min(a_abs, b_abs)
        max_val = max(a_abs, b_abs)

        ab_eps = (max_val - min_val) / min_val
        if ab_eps < eps:
            return True
    
    half_eps = eps / 2
    return a_abs < half_eps and b_abs < half_eps