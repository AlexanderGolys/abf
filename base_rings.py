from abstract import *
import functools, itertools
import numpy as np


class ZRing(BaseRing):
    def __init__(self):
        super().__init__(name='Z', euclidean=True, integral=True, pid=True, ufd=True, normed=True)

    def __call__(self, integer):
        return BaseElement(self, integer)

    def add(self, a, b):
        return self(a.value + b.value)

    def mul(self, a, b):
        return self(a.value * b.value)

    @property
    def one(self):
        return self(1)

    @property
    def zero(self):
        return self(0)

    def neg(self, a):
        return self(-a.value)

    def eq(self, a, b):
        return a.value == b.value

    def mod(self, a, b):
        raise NotImplementedError("Divisibility check is not implemented.")

    def maybe_unit_check(self, element):
        return element.value == 1 or element.value == -1

    def maybe_zero_check(self, element):
        return element.value == 0

    def force_zero_check(self, element):
        return element.value == 0

    def force_unit_check(self, element):
        return element.value == 1 or element.value == -1

    def truediv(self, a, b):
        if self.force_zero_check(a) and not self.force_zero_check(b):
            return self.zero
        if self.force_zero_check(b) and self.force_zero_check(a):
            raise ValueError("0/0")
        if a.value % b.value != 0:
            raise ValueError(f"{b} not divisible by {a}")
        return self(a.value // b.value)

    def longdiv(self, a, b):
        if self.force_zero_check(a) and not self.force_zero_check(b):
            return self.zero, self.zero
        if self.force_zero_check(b) and self.force_zero_check(a):
            raise ValueError("0/0")
        return self(a.value // b.value), self(a.value % b.value)

    def from_canonical_subring(self, element):
        return element

    def abs(self, a):
        return self(abs(a.value))

    @staticmethod
    def element_str(element):
        return str(element.value)


class QField(Field):
    def __init__(self):
        super().__init__([ZRing()], 'Q', characteristics=0, normed=True)

    def __call__(self, numerator, denominator=1):
        if denominator == 0:
            raise ValueError("Denominator cannot be zero.")
        if numerator == 0:
            return BaseElement(self, (0, 1))
        if numerator*denominator > 0:
            numerator = abs(numerator)
            denominator = abs(denominator)
        else:
            numerator = -abs(numerator)
            denominator = abs(denominator)
        gcd = np.gcd(numerator, denominator)
        return BaseElement(self, (numerator//gcd, denominator//gcd))

    def add(self, a, b):
        return self(a.value[0] * b.value[1] + b.value[0] * a.value[1], a.value[1] * b.value[1])

    def mul(self, a, b):
        return self(a.value[0] * b.value[0], a.value[1] * b.value[1])

    @property
    def one(self):
        return self(1)

    @property
    def zero(self):
        return self(0)

    def neg(self, a):
        return self(-a.value[0], a.value[1])

    def eq(self, a, b):
        return a.value == b.value

    def maybe_zero_check(self, element):
        return element.value[0] == 0

    def force_zero_check(self, element):
        return element.value[0] == 0

    def truediv(self, a, b):
        if self.force_zero_check(b):
            raise ZeroDivisionError
        return self(a.value[0] * b.value[1], a.value[1] * b.value[0])

    def from_canonical_subring(self, element):
        if element.ring == ZRing:
            return self(element.value)
        if element.ring == QField:
            return element
        raise ValueError("Only canonical subrings are Z and Q.")

    def abs(self, element):
        return abs(element.value[0]/element.value[1])

    @staticmethod
    def element_str(element):
        if element.value[1] == 1:
            return str(element.value[0])
        return str(element.value[0]) + '/' + str(element.value[1])


class RFloating(Field):

    def __init__(self):
        super().__init__([QField(), ZRing()], 'R', exact_values=False, characteristics=0, normed=True)

    def __call__(self, value, denominator=1):
        return BaseElement(self, value/denominator)

    def add(self, a, b):
        return self(a.value + b.value)

    def mul(self, a, b):
        return self(a.value * b.value)

    def neg(self, a):
        return self(-a.value)

    def maybe_zero_check(self, element):
        return np.isclose(element.value, 0)

    def force_zero_check(self, element):
        return np.isclose(element.value, 0)

    def truediv(self, a, b):
        return self(a.value / b.value)

    def from_canonical_subring(self, element):
        if element.ring == ZRing():
            return self(element.value)
        if element.ring == QField():
            return self(element.value[0]/element.value[1])
        if element.ring == RFloating():
            return element
        raise ValueError("Only canonical subrings are Z, Q and R.")

    def abs(self, element):
        return abs(element.value)

    @property
    def one(self):
        return self(1)

    @property
    def zero(self):
        return self(0)

    def eq(self, a, b):
        return np.isclose(a.value, b.value)

    @staticmethod
    def element_str(element):
        return str(element.value)


class CFloating(RFloating):
    def __init__(self):
        super().__init__()
        self.canonical_subrings = [ZRing(), QField(), RFloating()]
        self.name = 'C'
        self.properties['normal'] = True

    def __call__(self, real, imaginary=0, denominator=1):
        return BaseElement(self, complex(real, imaginary)/denominator)











