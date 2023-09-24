from abstract import *
from polynomials import *
from algebras import *
from base_rings import *


class Ideal(FGModule):
    def __init__(self, generators, name="", **properties):
        for g in generators:
            assert isinstance(g, BaseElement)

        generators = list(filter(lambda x: not x == generators[0].ring.zero, generators))

        if name == "":
            name = '(' + ', '.join([str(g) for g in generators]) + ')'

        properties['ideal'] = True
        super().__init__(generators[0].ring, generators, name=name, **properties)

    def scalar_mul(self, element, scalar):
        return self([scalar * c for c in element.coefficients])

    def as_ring_element(self, element):
        return sum([c * g for c, g in zip(element.coefficients, self.generators)], start=self.ring.zero)


class PolynomialIdeal(Ideal):
    def __init__(self, generators, name="", groebner=None, **properties):
        super().__init__(generators, name=name, **properties)
        self.groebner_basis = groebner

    def info(self):
        print(self.name + ': finitely generated ' + str(self.ring) + '-module')
        print(str(self.no_generators) + f' generator{"s" if self.no_generators > 1 else ""}')
        print('<' + ', '.join([str(g) for g in self.generators]) + '>')
        print(', '.join(filter(lambda x: self.properties[x], self.properties().keys())))
        if self.groebner_basis is not None:
            print('Groebner basis: (' + ', '.join([str(g) for g in self.groebner_basis]) + ')')
        print()


class KPolynomialIdeal(PolynomialIdeal):
    def __init__(self, generators, name="", groebner=None, find_groebner=True, **properties):
        super().__init__(generators, name=name, groebner=groebner, **properties)
        if find_groebner:
            self.to_groebner()

    @property
    def monomial_ideal(self):
        return all([g.monomial for g in self.generators])

    def S_polynomial(self, f, g):
        f_ = f if not isinstance(f, ModuleElement) else f()
        g_ = g if not isinstance(g, ModuleElement) else g()
        return self.ring.S_polynomial(f_, g_)

    def recalculate_groebner(self):
        self.groebner_basis = None
        self.to_groebner()

    def to_groebner(self):
        if self.groebner_basis is not None:
            return

        if self.no_generators == 1:
            self.groebner_basis = copy.copy(self.generators)
            return

        basis = copy.copy(self.generators)
        pairs = list(itertools.combinations(basis, 2))
        while pairs:
            f, g = pairs.pop()
            h = self.S_polynomial(f, g)
            r = self.reminder_wrt_family(h, basis)
            if r != self.ring.zero:
                pairs += [(r, b) for b in basis]
                basis.append(r)
        self.groebner_basis = basis

    def convert_basis_to_groebner(self):
        if self.groebner_basis is None:
            self.to_groebner()
        return self.__class__(self.groebner_basis, name=self.name, groebner=self.groebner_basis, find_groebner=False, **self.properties)

    def check_if_basis_groebner(self):
        for f, g in zip(self.generators):
            if f != g and self.S_polynomial(f, g) != self.ring.zero:
                return False
        return True

    def groebner_reminder(self, f):
        if self.groebner_basis is None:
            self.to_groebner()
        return self.ring.multi_long_div(f, self.groebner_basis)[1]

    def reminder_wrt_family(self, f, F):
        return self.ring.multi_long_div(f, F)[1]

    def check_if_belongs(self, element):
        self.to_groebner()
        return self.groebner_reminder(element) == self.ring.zero

    def _reduce_basis(self):
        if self.no_generators < 2:
            return
        for i, b in enumerate(self.generators):
            complementary_basis = [copy.copy(g) for g in self.generators if g != b]
            reduced_ideal = KPolynomialIdeal(complementary_basis)

            if reduced_ideal.check_if_belongs(b):
                return reduced_ideal._reduce_basis()
        return self

    def reduce_basis(self):
        reduced = self._reduce_basis()
        reduced.recalculate_groebner()
        return reduced

    def leading_monomial_ideal(self):
        self.to_groebner()
        return KMonomialIdeal(*[g.value.leading_monomial() for g in self.groebner_basis])

    def in_ideal(self):
        return self.leading_monomial_ideal().to_ideal()


class KMonomialIdeal:
    def __init__(self, *monomials):
        generators = []
        for m in monomials:
            if not isinstance(m, Polynomial) and m.monomial:
                generators.append(m.monomials[0])
            else:
                generators.append(m)
        self.monomials = generators

    def to_ideal(self):
        return KPolynomialIdeal(*[Polynomial(m) for m in self.monomials])





