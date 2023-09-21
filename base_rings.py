from abstract import *
import functools, itertools
import numpy as np


class BaseElement:
    def __init__(self, ring, value, name=None):
        self.ring = ring
        self.value = value
        self.name = name or str(value)

    def __str__(self):
        return self.name

    def __add__(self, other):
        return self.ring.add(self, other)

    def __mul__(self, other):
        return self.ring.mul(self, other)

    def __neg__(self):
        return self.ring.neg(self)

    def __eq__(self, other):
        return self.ring.eq(self, other)

    def __ne__(self, other):
        return not self.ring.eq(self, other)

    def __sub__(self, other):
        return self.ring.add(self, -other)

    def __truediv__(self, other):
        return self.ring.truediv(self, other)

    def __floordiv__(self, other):
        return self.ring.longdiv(self, other)

    def __mod__(self, other):
        return self.ring.mod(self, other)

    def __pow__(self, other):
        return functools.reduce(lambda x, y: x * y, [self] * other)

    def __abs__(self):
        return self.ring.abs(self)

    @functools.cached_property
    def grad(self):
        return self.ring.grad(self)


class Z(BaseRing):
    name = "Z"
    pid = True
    euclidean = True
    integral = True
    ufd = True
    normed = True

    def __call__(self, integer):
        return BaseElement(self, integer)

    def add(self, a, b):
        return self(a.value + b.value)

    def mul(self, a, b):
        return self(a.value * b.value)

    def one(self):
        return self(1)

    def zero(self):
        return self(0)

    def neg(self, a):
        return self(-a.value)

    def eq(self, a, b):
        return a.value == b.value

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
            return self.zero()
        if self.force_zero_check(b) and self.force_zero_check(a):
            raise ValueError("0/0")
        if a % b != 0:
            raise ValueError(f"{b} not divisible by {a}")
        return self(a.value // b.value)

    def longdiv(self, a, b):
        if self.force_zero_check(a) and not self.force_zero_check(b):
            return self.zero(), self.zero()
        if self.force_zero_check(b) and self.force_zero_check(a):
            raise ValueError("0/0")
        return self(a.value // b.value), self(a.value % b.value)

    def from_canonical_subring(self, element):
        return element

    def abs(self, a):
        return self(abs(a.value))


class Q(Field):
    field = True
    name = "Q"
    euclidean = True
    integral = True
    pid = True
    ufd = True
    local = True
    artinian = True
    normed = True

    def __init__(self):
        super().__init__([Z])

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

    def one(self):
        return self(1)

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
        if element.ring == Z:
            return self(element.value)
        if element.ring == Q:
            return element
        raise ValueError("Only canonical subrings are Z and Q.")

    def abs(self, element):
        return abs(element.value[0]/element.value[1])


class RFloating(Field):
    name = "R"

    def __init__(self):
        super().__init__([Q, Z])

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
        if element.ring == Z:
            return self(element.value)
        if element.ring == Q:
            return self(element.value[0]/element.value[1])
        if element.ring == RFloating:
            return element
        raise ValueError("Only canonical subrings are Z, Q and R.")

    def abs(self, element):
        return abs(element.value)

    def one(self):
        return self(1)

    def zero(self):
        return self(0)

    def eq(self, a, b):
        return np.isclose(a.value, b.value)


class CFloating(RFloating):
    name = "C"

    def __call__(self, real, imaginary=0, denominator=1):
        return BaseElement(self, complex(real, imaginary)/denominator)











