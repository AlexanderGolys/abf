import numpy as np


class Generator:
    def __init__(self, name, id=None, grading=0):
        self.id = id or name
        self.name = name
        self.grading = grading

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        return self.name


class BaseRing:
    def __init__(self, name, one, zero, valid_check=None, units_check=None):
        self.name = name
        self.add = lambda x, y: x + y
        self.multiply = lambda x, y: x * y
        self.negate = lambda x: -x
        self.one = one
        self.zero = zero
        if valid_check is not None:
            self.valid_check = lambda x: valid_check(x) or x.of == name
        else:
            self.valid_check = lambda x: x.of == name

        if units_check is not None:
            self.units_check = units_check
        else:
            self.units_check = lambda x: False if x == zero else None

    def __str__(self):
        return self.name

    def isunit(self):
        return self.units_check(self)



class Field(BaseRing):
    def __init__(self, name, one, zero, valid_check):
        super().__init__(name, one, zero, valid_check, lambda x: x != zero)


Z = BaseRing('Z', 1, 0, lambda x: isinstance(x, int), lambda x: abs(x) == 1)
R = Field('R', 1, 0, lambda x: isinstance(x, (int, float)))
Q = Field('Q', 1, 0, lambda x: isinstance(x, int))
C = Field('C', 1, 0, lambda x: isinstance(x, (int, float)))


class AlgebraElement:
    def __init__(self, algebra, polynomial):
        self.algebra = algebra
        self.polynomial = polynomial

    def __add__(self, other):
        if self.algebra != other.algebra:
            raise ValueError('Elements must be in the same algebra')
        return AlgebraElement(self.algebra, self.polynomial + other.polynomial)

    def __sub__(self, other):
        if self.algebra != other.algebra:
            raise ValueError('Elements must be in the same algebra')
        return AlgebraElement(self.algebra, self.polynomial - other.polynomial)

    def __mul__(self, other):
        if self.algebra != other.algebra:
            raise ValueError('Elements must be in the same algebra')
        return AlgebraElement(self.algebra, self.polynomial * other.polynomial)

    def __eq__(self, other):
        if self.algebra != other.algebra:
            if self.algebra == other.algebra.base_ring:
                return self == other.polynomial.constant_term
            if other.algebra == self.algebra.base_ring:
                return self.polynomial.constant_term == other
            return False
        return self.polynomial == other.polynomial

    def __str__(self):
        return str(self.polynomial)


class Monomial:
    def __init__(self, coefficient, generator_powers=None, coef_ring=C, var_names=None):
        self.coef_ring = coef_ring
        self.coefficient = coefficient
        self.exponent_index = np.array(generator_powers)
        self.no_variables = len(generator_powers)
        self.degree = sum(self.exponent_index)
        if coefficient == coef_ring.zero:
            self.degree = -1
        if var_names is None:
            if self.no_variables == 1:
                var_names = ['x']
            elif self.no_variables == 2:
                var_names = ['x', 'y']
            elif self.no_variables == 3:
                var_names = ['x', 'y', 'z']
            elif self.no_variables == 4:
                var_names = ['x', 'y', 'z', 'w']
            else:
                var_names = [f'x_{i}' for i in range(self.no_variables)]
        self.var_names = var_names

    def __str__(self):
        if self.coefficient == self.coef_ring.one:
            return f'{"".join([f"{var}^{e}" if e != 1 else f"{var}" for var, e in zip(self.var_names, self.exponent_index)])}'
        elif self.coefficient == -self.coef_ring.one:
            return f'-{"*".join([f"{var}^{e}" if e != 1 else f"{var}" for var, e in zip(self.var_names, self.exponent_index)])}'
        elif self.coefficient == self.coef_ring.zero:
            return '0'
        return f'{self.coefficient}{"*" if self.coefficient != self.coef_ring.one else ""}{"".join([f"{var}^{e}" if e != 1 else f"{var}" for var, e in zip(self.var_names, self.exponent_index)])}'

    def __eq__(self, other):
        if self.coef_ring != other.coef_ring:
            return False
        if self.coefficient != other.coefficient:
            return False
        if self.exponent_index != other.exponent_index:
            return False
        return True

    def __mul__(self, other):
        if self.coef_ring != other.coef_ring:
            raise ValueError('Monomials must be in the same coefficient ring')
        return Monomial(self.coefficient.multiply(other.coefficient),
                        self.exponent_index + other.exponent_index,
                        self.coef_ring,
                        self.var_names)

    def __hash__(self):
        return hash((self.coef_ring, self.coefficient, tuple(self.exponent_index)))

    def __call__(self, *args):
        if len(args) != self.no_variables:
            raise ValueError('Wrong number of arguments')
        return self.coefficient * np.prod([arg ** exp for arg, exp in zip(args, self.exponent_index)])

    @staticmethod
    def constant(coef_ring, coefficient, no_variables):
        return Monomial(coefficient, [0] * no_variables, coef_ring)

    def compare(self, other, order):
        # TODO: NEW CLASS FOR ORDER!!
        if order == 'lex':
            diff = self.exponent_index - other.exponent_index
            for d in diff:
                if d > 0:
                    return 1
                elif d < 0:
                    return -1
            return 0
        if order == 'grlex':
            if self.degree > other.degree:
                return 1
            elif self.degree < other.degree:
                return -1
            diff = self.exponent_index - other.exponent_index
            for d in diff:
                if d > 0:
                    return 1
                elif d < 0:
                    return -1
            return 0

        if order == 'grevlex':
            if self.degree > other.degree:
                return 1
            elif self.degree < other.degree:
                return -1
            diff = self.exponent_index - other.exponent_index
            for d in diff[::-1]:
                if d < 0:
                    return 1
                elif d > 0:
                    return -1
            return 0

    def check_divisibility(self, other):
        if self.coef_ring != other.coef_ring:
            raise ValueError('Monomials must be in the same coefficient ring')
        if self.no_variables != other.no_variables:
            raise ValueError('Monomials must have the same number of variables')
        diff = self.exponent_index - other.exponent_index
        if np.any(diff < 0):
            return False
        return True

    @staticmethod
    def create_one(no_variables, coef_ring=C):
        return Monomial(coef_ring.one, [0] * no_variables, coef_ring)

    @staticmethod
    def create_zero(no_variables, coef_ring=C):
        return Monomial(coef_ring.zero, [0] * no_variables, coef_ring)

    def __truediv__(self, other):
        if self.coef_ring != other.coef_ring:
            raise ValueError('Monomials must be in the same coefficient ring')
        if self.no_variables != other.no_variables:
            raise ValueError('Monomials must have the same number of variables')
        if not self.check_divisibility(other):
            raise ValueError('Divisor must divide dividend')
        return Monomial(self.coefficient/other.coefficient, self.exponent_index - other.exponent_index, self.coef_ring, self.var_names)

    def __pow__(self, other):
        if not isinstance(other, int) or other < 0:
            raise ValueError('Exponent must be a non-negative integer')
        if other == 0 and self.degree == -1:
            raise ValueError('0^0 is undefined')
        return Monomial(self.coefficient ** other, self.exponent_index * other, self.coef_ring, self.var_names)


class Polynomial:
    def __init__(self, monomials, order='grlex'):
        self.monomials = list(filter(lambda x: x.degree >= 0, monomials))
        self.no_variables = self.monomials[0].no_variables
        self.order = order
        self.sort()
        for i in range(len(self.monomials) - 1, 0, -1):
            if self.monomials[i].exponent_index == self.monomials[i - 1].exponent_index:
                self.monomials[i - 1].coefficient = self.monomials[i - 1].coefficient.add(self.monomials[i].coefficient)
                self.monomials.pop(i)
        self.value_dict = {tuple(m.exponent_index): m.coefficient for m in self.monomials}

    def sort(self):
        self.monomials.sort(key=lambda x: x.compare(x, self.order), reverse=True)

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
        return self.monomials[0].coef_ring.zero

    def __add__(self, other):
        if self.monomials[0].coef_ring != other.monomials[0].coef_ring:
            raise ValueError('Polynomials must be in the same coefficient ring')
        monomials = self.monomials
        for m in other.monomials:
            for included in monomials:
                if m.exponent_index == included.exponent_index:
                    included.coefficient += m.coefficient
                    break
            else:
                monomials.append(m)
        return Polynomial(monomials, self.order)

    @staticmethod
    def from_monomial(monomial, order='grlex'):
        return Polynomial([monomial], order)

    @staticmethod
    def from_constant(coefficient, no_variables, order='grlex'):
        return Polynomial.from_monomial(Monomial.constant(coefficient.coef_ring, coefficient, no_variables), order)

    def __mul__(self, other):
        if self.monomials[0].coef_ring != other.monomials[0].coef_ring:
            raise ValueError('Polynomials must be in the same coefficient ring')
        return Polynomial([m1 * m2 for m1 in self.monomials for m2 in other.monomials], self.order)

    def __call__(self, *args):
        if len(args) != self.no_variables:
            raise ValueError('Wrong number of arguments')
        return sum([m(*args) for m in self.monomials])








