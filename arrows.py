import types

from modules import *


class Morphism(ABC):
    def __init__(self, domain, codomain):
        self.domain = domain
        self.codomain = codomain

    @abstractmethod
    def __call__(self, element):
        pass


class AlgebraMorphism(Morphism):
    def __init__(self, domain, codomain, function_on_generators):
        super().__init__(domain, codomain)
        self.function_on_generators = function_on_generators
        self.images = [self.function_on_generators(g) for g in self.domain.generators]

    def __call__(self, element):
        return element.value(self.images)


class ModuleMorphism(Morphism):
    def __init__(self, domain, codomain, function_on_generators):
        super().__init__(domain, codomain)
        self.function_on_generators = function_on_generators
        self.images = [self.function_on_generators(g) for g in self.domain.generators]

    def __call__(self, element):
        return sum([c * g for c, g in zip(element.coefficients, self.images)], start=self.codomain.zero)

    def matrix(self):
        return [self(g) for g in self.domain.generators]


class CovariantFunctor(ABC):
    def __init__(self, domain, codomain, name='F'):
        self.domain = domain
        self.codomain = codomain
        self.name = name

    @abstractmethod
    def __call__(self, obj):
        pass

    @abstractmethod
    def induced_morphism(self, morphism):
        pass

    def __str__(self):
        return self.name + ': ' + str(self.domain) + '->' + str(self.codomain)


class IdentityFunctor(CovariantFunctor):
    def __init__(self, category):
        super().__init__(category, category, name='Id')

    def __call__(self, obj):
        return obj

    def induced_morphism(self, morphism):
        return morphism


class FreeAlgebra(CovariantFunctor):
    def __init__(self, base_ring, order=grlex):
        super().__init__('Set', Algebra, name='Free(' + str(base_ring) + ')')
        self.order = order
        self.base_ring = base_ring

    def __call__(self, finite_set):
        return PolynomialAlgebra(self.base_ring, list(map(str, finite_set)), order=self.order)

    def induced_morphism(self, function):
        return AlgebraMorphism(self(function.domain), self(function.codomain), function)


class ExtensionOfScalars(CovariantFunctor):
    def __init__(self, morphism_of_base_algebras):
        def extended_module(module):
            new_base = morphism_of_base_algebras.codomain


class TensorAlgebra(CovariantFunctor):
    def __init__(self, algebra):
        self.algebra = algebra
        self.base_ring = algebra.base_ring

        def function_on_objects(algebra):
            if isinstance(algebra, QuotientAlgebra) and algebra.base_ring == self.base_ring:
                pass  # TODO: need extension of scalars

        def function_on_morphisms(function):
            pass  # TODO: need extension of scalars




