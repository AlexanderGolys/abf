import warnings

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

    def quotient_algebra(self):
        return QuotientAlgebra(self.ring, self)

    @staticmethod
    def zero_deal(ring):
        return Ideal([ring.zero], name='(0)')


class PolynomialIdeal(Ideal):
    def __init__(self, generators, name="", groebner=None, **properties):
        super().__init__(generators, name=name, **properties)
        self.groebner_basis = groebner
        self.order = self.ring.order

    def info(self):
        print(self.name + ': finitely generated ' + str(self.ring) + '-module')
        print(str(self.no_generators) + f' generator{"s" if self.no_generators > 1 else ""}')
        print('<' + ', '.join([str(g) for g in self.generators]) + '>')
        print(', '.join(filter(lambda x: self.properties[x], self.properties().keys())))
        if self.groebner_basis is not None:
            print('Groebner basis: (' + ', '.join([str(g) for g in self.groebner_basis]) + ')')
        print()



class KPolynomialIdeal(PolynomialIdeal):
    def __init__(self, generators, name="", groebner=None, find_groebner=False, **properties):
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

    def check_if_groebner_minimal(self):
        if self.groebner_basis is None:
            return False
        for i, p in enumerate(self.groebner_basis):
            if p.value.leading_coefficient() != 1:
                return False
            for g in self.groebner_basis[i+1:]:
                if g.value.leading_monomial().divisible_by(p.value.leading_monomial()):
                    return False
                if p.value.leading_monomial().divisible_by(g.value.leading_monomial()):
                    return False
        return True

    def check_if_groebner_reduced(self):
        if self.groebner_basis is None:
            return False
        for i, p in enumerate(self.groebner_basis):
            leading_monomial = p.value.leading_monomial()
            if leading_monomial.coefficient != 1:
                return False
            for g in self.groebner_basis[:i] + self.groebner_basis[i+1:]:
                for m in g.value.monomials:
                    if m.divisible_by(leading_monomial):
                        return False
        return True

    def reduce_wrt_family(self, f, g, F=None):
        if F is None:
            F = self.groebner_basis
        lm_f, lm_g = f.value.leading_monomial(), g.value.leading_monomial()
        if lm_f.lcm(lm_g) == lm_f*lm_g:
            return self.ring.zero
        return self.ring.multi_long_div(self.S_polynomial(f, g), F)[1]

    def groebner_to_minimal(self):
        minimal = []
        for p in self.groebner_basis:
            complement = [g for g in self.groebner_basis if g != p]
            if not p != 0 and self.reminder_wrt_family(p, complement) == self.ring.zero:
                minimal.append(p/p.value.leading_coefficient())
        self.groebner_basis = minimal

    def groebner_to_reduced(self):
        self.groebner_to_minimal()
        G = self.groebner_basis
        H = []
        for i, p in enumerate(G):
            h = self.reminder_wrt_family(p, H + G[:i+1])
            H.append(h)
        self.groebner_basis = H
        self.order.sort(self.groebner_basis)

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
            reminder = self.reduce_wrt_family(f, g, basis)
            if reminder != 0:
                pairs += [(reminder, b) for b in basis]
                basis.append(reminder)
        self.groebner_basis = basis
        self.groebner_to_reduced()

    def convert_basis_to_groebner(self):
        if self.groebner_basis is None:
            self.to_groebner()
        return self.__class__(self.groebner_basis, name=self.name, groebner=self.groebner_basis, find_groebner=False, **self.properties)

    def check_if_basis_groebner(self, basis=None):
        if basis is None:
            basis = self.generators
        pairs = list(itertools.combinations(basis, 2))
        while pairs:
            f, g = pairs.pop()
            h = self.S_polynomial(f, g)
            r = self.reminder_wrt_family(h, basis)
            if r != self.ring.zero:
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

    def __contains__(self, item):
        return self.check_if_belongs(item)

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

    def leading_ideal(self):
        return self.leading_monomial_ideal().to_ideal()

    def __eq__(self, other):
        if not isinstance(other, KPolynomialIdeal):
            return False
        if self.groebner_basis is None:
            self.to_groebner()
        if other.groebner_basis is None:
            other.to_groebner()
        if self.order != other.order:
            warnings.warn('Ideals are not equal due to different orders.')
            return False
        return self.groebner_basis == other.groebner_basis

    def quotient_algebra(self):
        return QuotientKAlgebra(self.ring, self)


class KMonomialIdeal:
    def __init__(self, *monomials):
        generators = []
        for m in monomials:
            if not isinstance(m, Polynomial) and m.monomial:
                generators.append(m.monomials[0])
            else:
                generators.append(m)
        self.monomials = generators
        self.monomials = [m/m.coefficient for m in self.monomials if m.coefficient != 0]

    def to_ideal(self):
        return KPolynomialIdeal(*[Polynomial(m) for m in self.monomials])

    def __contains__(self, item):
        return any([item.value.divisible_by(m) for m in self.monomials])







