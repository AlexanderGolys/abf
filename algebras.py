from base_rings import *
from polynomials import *
from modules import *
import modules


class AlgebraFG(BaseRing):
    def __init__(self, name, base_ring, generators, order=grlex, **properties):
        canonical_subrings = base_ring.canonical_subrings + [base_ring]
        super().__init__(canonical_subrings, name, **properties)
        self.base_ring = base_ring
        self.generators = generators
        self.no_generators = len(generators)
        self.order = order

    def __call__(self, polynomial):
        if isinstance(polynomial, Monomial):
            return BaseElement(self, Polynomial(polynomial, order=self.order))
        polynomial.order = self.order
        return BaseElement(self, polynomial)

    def add(self, a, b):
        return self(a.value + b.value)

    def mul(self, a, b):
        return self(a.value * b.value)

    @property
    def one(self):
        return self(Polynomial.create_one(self.no_generators, self.base_ring, self.order))

    @property
    def zero(self):
        return self(Polynomial.create_zero(self.no_generators, self.base_ring, self.order))

    def neg(self, a):
        return self(-a.value)

    def eq(self, a, b):
        return a.value == b.value

    @staticmethod
    def element_str(element):
        return str(element.value)

    def from_canonical_subring(self, element):
        if element.ring in self.base_ring.canonical_subrings:
            element = self.base_ring.from_canonical_subring(element)
        if element.ring == self:
            return element
        return self(Polynomial.from_constant(element, self.no_generators, self.order))

    def info(self):
        print(self.name + ': finitely generated ' + str(self.base_ring) + '-algebra')
        print(str(self.no_generators) + ' generators')


class AlgebraFP(AlgebraFG):
    def __init__(self, name, base_ring, generators, relations, order=grlex, **properties):
        super().__init__(name, base_ring, generators, order, **properties)
        self.relations = relations
        self.no_relations = len(relations)

    def ideal(self, name=None):
        return modules.Ideal(self.relations, name)

    def info(self):
        print(self.name + ': finitely presented ' + str(self.base_ring) + '-algebra')
        print(str(self.no_generators) + ' generators, ' + str(len(self.relations)) + ' relations')


class PolynomialAlgebra(AlgebraFP):
    def __init__(self, base_ring, no_variables, variable_names=None, order=grlex):
        if variable_names is None:
            if no_variables < 5:
                variable_names = ['x', 'y', 'z', 'w'][:no_variables]
            else:
                variable_names = [f'x_{i}' for i in range(no_variables)]
        name = str(base_ring) + '[' + ', '.join(variable_names) + ']'
        AlgebraFP.__init__(self, name, base_ring, variable_names, [], order)
        if self.no_generators == 1 and self.base_ring.properties['pid']:
            self.properties['pid'] = True
        if self.base_ring.properties['ufd']:
            self.properties['ufd'] = True
        if self.base_ring.properties['integral']:
            self.properties['integral'] = True
        if self.base_ring.properties['noetherian']:
            self.properties['noetherian'] = True

    def grad(self, element):
        return element.value.degree

    def element_str(self, element):
        return str(element.value)

    def truediv(self, f, g):
        return self(f.value/g.value)

    def lcm_lead(self, f, g):
        leading_f = f.value.leading_monomial()
        leading_g = g.value.leading_monomial()
        lcm_index = [max(e1, e2) for e1, e2 in zip(leading_f.exponent_index, leading_g.exponent_index)]
        return self(Polynomial(Monomial(self.base_ring.one, lcm_index, leading_f.var_names), order=self.order))

    def S_polynomial(self, f, g):
        if f == self.zero or g == self.zero:
            raise ValueError("S-polynomial of zero polynomials is not defined.")
        lcm = self.lcm_lead(f, g)
        return lcm/self(f.value.leading_monomial()) * f - lcm/self(g.value.leading_monomial()) * g


class KAlgebraFP(AlgebraFP):
    def __init__(self, name, base_ring, generators, relations, order=grlex, **properties):
        assert base_ring.properties['field']
        AlgebraFP.__init__(self, name, base_ring, generators, relations, order, **properties)


class KPolynomialAlgebra(PolynomialAlgebra, KAlgebraFP):
    def __init__(self, field, no_variables, variable_names=None, order=grlex):
        PolynomialAlgebra.__init__(self, field, no_variables, variable_names, order)
        KAlgebraFP.__init__(self, self.name, field, self.generators, [], order)

    def multi_long_div(self, f, g_list):
        a = [copy.copy(self.zero)]*len(g_list)
        r = copy.copy(self.zero)
        p = self(Polynomial(*f.value.monomials))
        while not p == self.zero:
            i = 0
            occurred = False
            while i < len(g_list) and not occurred:
                if g_list[i] != self.zero and p.value.leading_monomial().divisible_by(g_list[i].value.leading_monomial()):
                    occurred = True
                    quotient = self(p.value.leading_monomial() / g_list[i].value.leading_monomial())
                    a[i] = a[i] + quotient
                    p = p - g_list[i] * quotient
                else:
                    i += 1
            if not occurred:
                r = r + self(Polynomial(p.value.leading_monomial()))
                p = p - self(Polynomial(p.value.leading_monomial()))
        return a, r


class Matrix:
    def __init__(self, coefficients):
        self.coefficients = coefficients
        self.no_rows = len(coefficients)
        self.no_columns = len(coefficients[0])
        self.ring = coefficients[0][0].ring

    def __add__(self, other):
        assert self.no_rows == other.no_rows and self.no_columns == other.no_columns
        return Matrix([[self.coefficients[i][j] + other.coefficients[i][j] for j in range(self.no_columns)] for i in range(self.no_rows)])

    def __mul__(self, scalar):
        return Matrix([[scalar * self.coefficients[i][j] for j in range(self.no_columns)] for i in range(self.no_rows)])

    def __matmul__(self, other):
        assert self.no_columns == other.no_rows
        return Matrix([[sum([self.coefficients[i][k] * other.coefficients[k][j] for k in range(self.no_columns)]) for j in range(other.no_columns)] for i in range(self.no_rows)])

    def vanishing_ideal_of_coefs(self):
        return modules.Ideal([self.coefficients[i][j] for i in range(self.no_rows) for j in range(self.no_columns)], name='vanishing ideal of coefficients')

    def submatrix(self, rows, columns, keep=False):
        if isinstance(rows, int):
            rows = [rows]
        if isinstance(columns, int):
            columns = [columns]
        if keep:
            return Matrix([[self.coefficients[i][j] for j in range(self.no_columns) if j in columns] for i in range(self.no_rows) if i in rows])
        return Matrix([[self.coefficients[i][j] for j in range(self.no_columns) if j not in columns] for i in range(self.no_rows) if i not in rows])

    def det(self):
        assert self.no_rows == self.no_columns
        if self.no_rows == 1:
            return self.coefficients[0][0]
        if self.no_rows == 2:
            return self.coefficients[0][0] * self.coefficients[1][1] - self.coefficients[0][1] * self.coefficients[1][0]
        return sum([(-1)**i*self.ring.one * self.coefficients[0][i] * self.submatrix(0, i).det() for i in range(self.no_columns)])

    def __str__(self):
        return '\n'.join(['|' + ', '.join([str(self.coefficients[i][j]) for j in range(self.no_columns)]) + '|' for i in range(self.no_rows)])

    def __pow__(self, power):
        assert self.no_rows == self.no_columns
        if power == 0:
            return Matrix([[self.ring.one if i == j else self.ring.zero for j in range(self.no_columns)] for i in range(self.no_rows)])
        if power == 1:
            return self
        if power % 2 == 0:
            return (self@self)**(power//2)
        return self@(self@self)**((power-1)//2)

    def zero_matrix(self, no_rows, no_columns):
        return Matrix([[self.ring.zero for j in range(no_columns)] for i in range(no_rows)])

    def identity_matrix(self, no_rows):
        return Matrix([[self.ring.one if i == j else self.ring.zero for j in range(no_rows)] for i in range(no_rows)])

    def __eq__(self, other):
        if other == 0:
            return self == self.zero_matrix(self.no_rows, self.no_columns)
        if other == 1:
            return self.is_square and self == self.identity_matrix(self.no_rows)
        if self.no_rows != other.no_rows or self.no_columns != other.no_columns:
            return False
        return all([self.coefficients[i][j] == other.coefficients[i][j] for i in range(self.no_rows) for j in range(self.no_columns)])

    def __ne__(self, other):
        return not self == other

    def __neg__(self):
        return self * -self.ring.one

    def __sub__(self, other):
        return self + (-other)

    def transpose(self):
        return Matrix([[self.coefficients[j][i] for j in range(self.no_rows)] for i in range(self.no_columns)])

    @property
    def is_square(self):
        return self.no_rows == self.no_columns

    def trace(self):
        assert self.is_square
        return sum([self.coefficients[i][i] for i in range(self.no_rows)])


class QuotientAlgebra(AlgebraFP):
    def __init__(self, ideal, name=None):
        if name is None:
            name = str(ideal.ring) + '/' + ideal.name
        super().__init__(name, ideal.ring, ideal.ring.generators, ideal.generators)
        self.ideal = ideal

    @staticmethod
    def element_str(element):
        return str(element.value) + ' + ' + str(element.ideal)


class QuotientKAlgebra(QuotientAlgebra):
    def __call__(self, polynomial):
        if isinstance(polynomial, Monomial):
            return BaseElement(self, Polynomial(polynomial, order=self.order))
        polynomial.order = self.order
        if self.ideal.groebner_basis is not None:
            polynomial = self.ideal.groebner_reminder(polynomial)
        else:
            polynomial = self.ideal.reminder_wrt_family(polynomial, self.ideal.generators)
        return BaseElement(self, polynomial)







