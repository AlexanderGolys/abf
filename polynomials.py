import functools

import numpy as np
from base_rings import *


class MonomialOrder:
    def __init__(self, name, comparator):
        self.name = name
        self.comparator = comparator

    def __call__(self, a, b):
        return self.comparator(a.exponent_index, b.exponent_index)

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return self.name == other.name

    def sort(self, monomials):
        return sorted(monomials, key=functools.cmp_to_key(lambda x, y: self(x, y)), reverse=True)


def lex_comp(l1, l2):
    diff = l1 - l2
    for d in diff:
        if d > 0:
            return 1
        if d < 0:
            return -1
    return 0


def grlex_comp(l1, l2):
    if sum(l1) > sum(l2):
        return 1
    if sum(l1) < sum(l2):
        return -1
    return lex_comp(l1, l2)


def grevlex_comp(l1, l2):
    if sum(l1) > sum(l2):
        return 1
    if sum(l1) < sum(l2):
        return -1
    diff = l1 - l2
    for d in diff[::-1]:
        if d > 0:
            return -1
        if d < 0:
            return 1
    return 0


grlex = MonomialOrder('grlex', grlex_comp)
lex = MonomialOrder('lex', lex_comp)
grevlex = MonomialOrder('grevlex', grevlex_comp)


class Monomial:
    def __init__(self, coefficient, generator_powers, var_names=None, infinite_variables=False):
        self.ring = coefficient.ring
        self.coefficient = coefficient
        self.exponent_index = np.array(generator_powers).astype(int)
        self.no_variables = len(generator_powers) if not infinite_variables else None
        self.degree = sum(self.exponent_index)
        self.infinite_variables = infinite_variables
        if coefficient == self.ring.zero:
            self.degree = -1
            self.exponent_index = np.array([0] * self.no_variables)
        if var_names is None:
            if self.no_variables < 5:
                var_names = ['x', 'y', 'z', 'w'][:self.no_variables]
            else:
                var_names = [f'x_{i}' for i in range(self.no_variables)]
        self.var_names = var_names

    def __str__(self):
        if self.degree < 1:
            return str(self.coefficient)
        if self.coefficient == self.ring.one:
            return f'{"".join([f"{var}^{e}" if e > 1 else (f"{var}" if e != 0 else "") for var, e in zip(self.var_names, self.exponent_index)])}'
        elif self.coefficient == -self.ring.one:
            return f'-{"".join([f"{var}^{e}" if e > 1 else (f"{var}" if e != 0 else "") for var, e in zip(self.var_names, self.exponent_index)])}'
        elif self.coefficient == self.ring.zero:
            return '0'
        return f'{self.coefficient if self.coefficient != self.ring.one else ""}{"".join([f"{var}^{e}" if e > 1  else (f"{var}" if e != 0 else "") for var, e in zip(self.var_names, self.exponent_index)])}'

    def __eq__(self, other):
        if self.ring != other.ring:
            return False
        if self.coefficient != other.coefficient:
            return False
        if self.no_variables is None or other.no_variables is None:
            return list(self.exponent_index) + [self.ring.zero]*len(other.exponent_index) == list(other.exponent_index) + [self.ring.zero]*len(self.exponent_index)
        return self.exponent_index == other.exponent_index

    def __mul__(self, other):
        assert isinstance(other, (Monomial, BaseElement))
        if isinstance(other, BaseElement):
            return Monomial(self.coefficient * other, self.exponent_index, self.var_names)
        if self.ring != other.ring:
            raise ValueError('Monomials must be in the same coefficient ring')
        if self.no_variables is None or other.no_variables is None:
            return Monomial(self.coefficient * other.coefficient,
                            np.array(list(self.exponent_index) + [self.ring.zero]*max(0, (len(other.exponent_index) - len(self.exponent_index))))
                            + np.array(list(other.exponent_index) + [self.ring.zero]*max(0, (len(self.exponent_index) - len(other.exponent_index)))),
                            self.var_names)

        return Monomial(self.coefficient*other.coefficient,
                        self.exponent_index + other.exponent_index,
                        self.var_names)

    def __call__(self, *args):
        if self.no_variables is not None and len(args) != self.no_variables:
            raise ValueError('Wrong number of arguments')
        if self.no_variables is None:
            exponent_index = list(self.exponent_index) + [self.ring.zero]*max(0, len(args)-len(self.exponent_index))
            return self.coefficient * np.prod([arg ** exp for arg, exp in zip(args, exponent_index)])
        return self.coefficient * np.prod([arg ** exp for arg, exp in zip(args, self.exponent_index)])

    @staticmethod
    def constant(coefficient, no_variables):
        if no_variables is None:
            return Monomial(coefficient, [0], infinite_variables=True)
        return Monomial(coefficient, [0] * no_variables)

    def compatibilise_infinite_variables(self, other):
        if self.infinite_variables and len(self.exponent_index) < len(other.exponent_index):
            self.exponent_index = np.array(list(self.exponent_index) + [self.ring.zero]*(len(other.exponent_index) - len(self.exponent_index)))
        if other.infinite_variables and len(other.exponent_index) < len(self.exponent_index):
            other.exponent_index = np.array(list(other.exponent_index) + [other.ring.zero]*(len(self.exponent_index) - len(other.exponent_index)))

    def compare(self, other, order):
        self.compatibilise_infinite_variables(other)
        return order(self, other)

    def divisible_by(self, other):
        self.compatibilise_infinite_variables(other)
        if self.ring != other.ring:
            raise ValueError('Monomials must be in the same coefficient ring')
        if self.no_variables != other.no_variables and not (self.infinite_variables or other.infinite_variables):
            raise ValueError('Monomials must have the same number of variables')
        diff = self.exponent_index - other.exponent_index
        if np.any(diff < 0):
            return False
        return True

    def divides(self, other):
        self.compatibilise_infinite_variables(other)
        return other.divisible_by(self)

    @staticmethod
    def create_one(no_variables, coef_ring):
        if no_variables is None:
            return Monomial(coef_ring.one, [0], infinite_variables=True)
        return Monomial(coef_ring.one, [0] * no_variables)

    @staticmethod
    def create_zero(no_variables, coef_ring):
        if no_variables is None:
            return Monomial(coef_ring.zero, [0], infinite_variables=True)
        return Monomial(coef_ring.zero, [0] * no_variables)

    def __truediv__(self, other):
        if self.ring != other.ring:
            raise ValueError('Monomials must be in the same coefficient ring')
        if self.no_variables != other.no_variables and not (self.infinite_variables or other.infinite_variables):
            raise ValueError('Monomials must have the same number of variables')
        if not self.divisible_by(other):
            raise ValueError('Divisor must divide dividend')
        if other.coefficient == self.ring.zero:
            raise ZeroDivisionError('Coefficient is zero')
        return Monomial(self.coefficient/other.coefficient, self.exponent_index - other.exponent_index, self.var_names)

    def __pow__(self, other):
        if not isinstance(other, int) or other < 0:
            raise ValueError('Exponent must be a non-negative integer')
        if other == 0 and self.degree == -1:
            raise ValueError('0^0 is undefined')
        if other == 0:
            return Monomial.create_one(self.no_variables, self.ring)
        return Monomial(self.coefficient ** other, self.exponent_index * other, self.var_names)


class Polynomial:
    def __init__(self, *monomials, order=grlex):
        self.monomials = list(filter(lambda x: x.degree >= 0 and not x.coefficient.maybe_is_zero(), monomials))
        if len(self.monomials) == 0:
            self.monomials.append(Monomial(monomials[0].ring.zero, [0] * monomials[0].no_variables))
        self.no_variables = self.monomials[0].no_variables
        self.order = order
        self.monomials = sorted(self.monomials, key=functools.cmp_to_key(lambda x, y: x.compare(y, order=self.order)), reverse=True)
        self.infinite_variables = self.monomials[0].infinite_variables
        for i in range(len(self.monomials) - 1, 0, -1):
            if np.array_equal(self.monomials[i].exponent_index, self.monomials[i - 1].exponent_index):
                self.monomials[i - 1].coefficient += self.monomials[i].coefficient
                self.monomials.pop(i)
        self.value_dict = {tuple(m.exponent_index): m.coefficient for m in self.monomials}
        self.ring = self.monomials[0].ring

    @property
    def degree(self):
        return max([m.degree for m in self.monomials if not m.coefficient.maybe_is_zero()] + [-1])

    def __str__(self):
        return ' + '.join([str(m) for m in self.monomials])

    def multidegree(self):
        return self.monomials[0].degree

    def leading_coefficient(self):
        return self.monomials[0].coefficient

    def leading_monomial(self):
        return self.monomials[0]

    def constant_term(self):
        if tuple([0]*self.no_variables) in self.value_dict:
            return self.value_dict[tuple([0]*self.no_variables)]
        return self.monomials[0].ring.zero

    def __add__(self, other):
        assert isinstance(other, (Polynomial, Monomial))
        other = other if not isinstance(other, Monomial) else Polynomial.from_monomial(other)
        if self.monomials[0].ring != other.monomials[0].ring:
            raise ValueError('Polynomials must be in the same coefficient ring')
        monomials = self.monomials
        for m in other.monomials:
            duplicate = False
            for i, included in enumerate(monomials):
                if np.array_equal(m.exponent_index, included.exponent_index):
                    monomials[i] = Monomial(included.coefficient + m.coefficient, m.exponent_index, m.var_names)
                    duplicate = True
                    break
            if not duplicate:
                monomials.append(m)
        return Polynomial(*monomials, order=self.order)

    @staticmethod
    def from_monomial(monomial, order=grlex):
        return Polynomial(monomial, order=order)

    @staticmethod
    def from_constant(coefficient, no_variables, order=grlex):
        return Polynomial.from_monomial(Monomial.constant(coefficient, no_variables), order)

    def __mul__(self, other):
        if isinstance(other, BaseElement) and other.ring == self.ring:
            return Polynomial(*[m * other for m in self.monomials], order=self.order)

        other = other if not isinstance(other, Monomial) else Polynomial(other)

        if self.ring != other.ring:
            raise ValueError('Polynomials must be in the same coefficient ring')
        return Polynomial(*[m1 * m2 for m1 in self.monomials for m2 in other.monomials], order=self.order)

    def __neg__(self):
        return self * -self.ring.one

    def __eq__(self, other):
        return all([m1 == m2 for m1, m2 in zip(self.monomials, other.monomials)]) or self.degree == other.degree == -1

    def __sub__(self, other):
        return self + (-other)

    def __truediv__(self, other):
        if isinstance(other, Monomial):
            return Polynomial(*[m/other for m in self.monomials], order=self.order)
        if isinstance(other, Polynomial) and len(other.monomials) == 1:
            return Polynomial(*[m/other.monomials[0] for m in self.monomials], order=self.order)
        raise NotImplementedError('Division of polynomials is implemented only by monomials')

    def __call__(self, *args):
        if len(args) != self.no_variables:
            raise ValueError('Wrong number of arguments')
        return sum([m(*args) for m in self.monomials])

    @staticmethod
    def create_one(no_variables, coef_ring, order=grlex):
        return Polynomial.from_monomial(Monomial.create_one(no_variables, coef_ring), order)

    @staticmethod
    def create_zero(no_variables, coef_ring, order=grlex):
        return Polynomial.from_monomial(Monomial.create_zero(no_variables, coef_ring), order)

    @property
    def monomial(self):
        return len(self.monomials) == 1







