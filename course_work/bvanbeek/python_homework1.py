"""Write a function to compute a definite integral of an aribtrary function (also taking an arbitrary number of parameters) on a given interval according to one of the known quadrature rules, such as Simpson rule: https://en.wikipedia.org/wiki/Simpson%27s_rule.

The more general the function (e.g. for functions of multiple dimensions), the better.

"""
import unittest
from itertools import combinations_with_replacement


def get_simspon_integral(func, a, b, *args, **kwargs):
    """Calculate the intergral of `func` over the `[a, b]` interval using Simpson's rule."""
    if a > b:
        raise ValueError("`a` must be smaller than `b`")
    elif a == b:
        return 0.0

    f = lambda x: func(x, *args, **kwargs)
    ret = (b - a) / 6
    ret *= f(a) + 4 * f((a + b) / 2) + f(b)
    return ret


class TestGetSimpsonIntegral(unittest.TestCase):
    def test_linear(self, n=5):
        func = lambda x: x
        for start, stop in combinations_with_replacement(range(n), 2):
            integral = get_simspon_integral(func, start, stop)
            integral_ref = (stop**2 / 2) - (start**2 / 2)
            self.assertAlmostEqual(integral, integral_ref)


unittest.main(argv=[''], verbosity=2, exit=False)
