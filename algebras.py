from polynomials import *
import modules
import arrows


class Algebra(BaseRing):
    def __init__(self, name, base_ring, order=grlex, **properties):
        canonical_subrings = base_ring.canonical_subrings + [base_ring]
        super().__init__(canonical_subrings, name, **properties)
        self.base_ring = base_ring
        self.order = order
        self.no_generators = None

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
        return self(Polynomial.create_one(None, self.base_ring, self.order))

    @property
    def zero(self):
        return self(Polynomial.create_zero(None, self.base_ring, self.order))

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
        return self(Polynomial.from_constant(element, None, self.order))

    def info(self):
        print(self.name + ': ' + str(self.base_ring) + '-algebra')

    def localisation_at_element(self, element):
        if element == self.zero:
            return AlgebraFG('0', self.base_ring, [])
        poly = Polynomial(Monomial(element, [1], ['t']), Monomial(-self.one, [0], ['t']))
        free_alg = PolynomialAlgebra(self, 1, ['t'])
        return modules.Ideal([free_alg(poly)]).quotient_algebra()

    def localisation_morphism_at_element(self, element):
        localisation = self.localisation_at_element(element)
        return arrows.AlgebraMorphism(self, localisation, lambda x: x)


class AlgebraFG(Algebra):
    def __init__(self, name, base_ring, generators, order=grlex, **properties):
        super().__init__(name, base_ring, order, **properties)
        self.generators = generators
        self.no_generators = len(generators)

    def info(self):
        print(self.name + ': finitely generated ' + str(self.base_ring) + '-algebra')
        print(str(self.no_generators) + ' generators')

    @property
    def one(self):
        return self(Polynomial.create_one(self.no_generators, self.base_ring, self.order))

    @property
    def zero(self):
        return self(Polynomial.create_zero(self.no_generators, self.base_ring, self.order))

    @property
    def generator_elements(self):
        return [self(Polynomial(Monomial(self.base_ring.one, [1 if j == i else 0 for j in range(self.no_generators)], self.generators))) for i, g in enumerate(self.generators)]


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

    def to_quotient_algebra(self):
        return QuotientAlgebra(self.ideal(), name=self.name)


class KAlgebraFP(AlgebraFP):
    def __init__(self, name, base_ring, generators, relations, order=grlex, **properties):
        assert base_ring.properties['field']
        AlgebraFP.__init__(self, name, base_ring, generators, relations, order, **properties)

    def to_quotient_algebra(self):
        return QuotientKAlgebra(self.ideal(), name=self.name)


class QuotientAlgebra(AlgebraFP):
    def __init__(self, ideal, name=None):
        if name is None:
            name = ideal.ring.name + '/' + ideal.name
        super().__init__(name, ideal.ring.base_ring, ideal.ring.generators, ideal.generators, order=ideal.ring.order)
        self.ideal = ideal

    @staticmethod
    def element_str(element):
        return str(element.value) + ' + ' + str(element.ring.ideal)


class QuotientKAlgebra(QuotientAlgebra, KAlgebraFP):
    def __init__(self, ideal, name=None):
        QuotientAlgebra.__init__(self, ideal, name)
        KAlgebraFP.__init__(self, name, ideal.ring, ideal.ring.generators, ideal.generators, order=ideal.order)

    def __call__(self, polynomial):
        if isinstance(polynomial, Monomial):
            return BaseElement(self, Polynomial(polynomial, order=self.order))
        polynomial.order = self.order
        if self.ideal.groebner_basis is not None:
            polynomial = self.ideal.groebner_reminder(polynomial)
        else:
            polynomial = self.ideal.reminder_wrt_family(polynomial, self.ideal.generators)
        return BaseElement(self, polynomial)


class PolynomialAlgebra(AlgebraFP):
    def __init__(self, base_ring, no_variables, variable_names=None, order=grlex):
        if variable_names is None:
            if no_variables < 5:
                variable_names = ['x', 'y', 'z', 'w'][:no_variables]
            else:
                variable_names = [f'x_{i}' for i in range(no_variables)]
        name = str(base_ring) + '[' + ', '.join(variable_names) + ']'

        if isinstance(base_ring, PolynomialAlgebra):
            self.__init__(base_ring.base_ring, no_variables + base_ring.no_generators,
                          variable_names + base_ring.generators, order=order)
            return

        super().__init__(name, base_ring, variable_names, [], order)
        if self.no_generators == 1 and self.base_ring.properties['field']:
            self.properties['pid'] = True
        if self.base_ring.properties['ufd']:
            self.properties['ufd'] = True
        if self.base_ring.properties['integral']:
            self.properties['integral'] = True
        if self.base_ring.properties['noetherian']:
            self.properties['noetherian'] = True

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


class KPolynomialAlgebra(PolynomialAlgebra, QuotientKAlgebra):
    def __init__(self, field, no_variables, variable_names=None, order=grlex):
        PolynomialAlgebra.__init__(self, field, no_variables, variable_names, order)
        QuotientKAlgebra.__init__(self, modules.Ideal.zero_ideal(self), name=self.name)

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
        return Matrix([[sum([self.coefficients[i][k] * other.coefficients[k][j] for k in range(self.no_columns)], self.ring.one) for j in range(other.no_columns)] for i in range(self.no_rows)])

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

    @staticmethod
    def fill_free_coefficients(ring, *shapes, order=grlex):
        matrices = []
        shapes = [(n, n) if isinstance(n, int) else n for n in shapes]
        no_vars = sum([s[0]*s[1] for s in shapes])
        varnames = []
        if len(shapes) < 5:
            names = ['x', 'y', 'z', 'w'][:len(shapes)]
            for shape, name in zip(shapes, names):
                if shape[0] == shape[1] == 1:
                    varnames += [name]
                elif shape[0] == 1 or shape[1] == 1:
                    varnames += [name + '_' + str(i) for i in range(max(shape[0], shape[1]))]
                else:
                    varnames += [name + '_' + str(i)+str(j) for i in range(shapes[0][0]) for j in range(shapes[0][1])]
        else:
            for k, shape in enumerate(shapes):
                name = 'x_' + str(k)
                if shape[0] == shape[1] == 1:
                    varnames += [name]
                elif shape[0] == 1 or shape[1] == 1:
                    varnames += [name + '_' + str(i) for i in range(max(shape[0], shape[1]))]
                else:
                    varnames += [name + '_' + str(i)+str(j) for i in range(shapes[0][0]) for j in range(shapes[0][1])]
        algebra = PolynomialAlgebra(ring, no_vars, varnames, order=order)

        used_vars = 0
        for shape in shapes:
            coefficients = []
            for i in range(shape[1]):
                row = algebra.generator_elements[used_vars + i*shape[1]:used_vars + (i+1)*shape[1]]
                coefficients.append(row)
            used_vars += shape[0]*shape[1]
            matrices.append(Matrix(coefficients))
        return matrices

    def vanishing_ideal(self, name=None):
        coefs = sum(self.coefficients, [])
        return modules.Ideal(coefs, name=name)

    def vanishing_algebra(self, name=None):
        return self.vanishing_ideal().quotient_algebra(name)










