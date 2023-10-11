from modules import *
from algebras import *


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

    @staticmethod
    def identity(algebra):
        return AlgebraMorphism(algebra, algebra, lambda x: x)


class ModuleMorphism(Morphism):
    def __init__(self, domain, codomain, function_on_generators):
        super().__init__(domain, codomain)
        self.function_on_generators = function_on_generators
        self.images = [self.function_on_generators(g) for g in self.domain.generators]

    def __call__(self, element):
        return sum([c * g for c, g in zip(element.coefficients, self.images)], start=self.codomain.zero)

    def matrix(self):
        return [self(g) for g in self.domain.generators]

    def __str__(self):
        return str(self.matrix())

    @staticmethod
    def identity(module):
        return ModuleMorphism(module, module, lambda x: x)


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
        return self.name + ': ' + str(self.domain) + ' -> ' + str(self.codomain)


class ContravariantFunctor(ABC):
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
        return self.name + ': ' + str(self.domain) + '^op -> ' + str(self.codomain)


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
        super().__init__(Module, Module, name='-(x)' + str(morphism_of_base_algebras.codomain))
        self.new_base = morphism_of_base_algebras.codomain
        self.morphism_of_base_algebras = morphism_of_base_algebras

    def __call__(self, module):
        return Module(self.codomain, module.generators, name=module.name,
                      custom_mul=lambda c, g: self.morphism_of_base_algebras(c), **module.properties)

    def induced_morphism(self, morphism):
        return ModuleMorphism(self(morphism.domain), self(morphism.codomain), morphism.function_on_generators)


class TensorAlgebraOverBase(CovariantFunctor):
    def __init__(self, new_algebra):
        super().__init__(Algebra, Algebra, name='-(x)' + str(new_algebra))
        self.new_algebra = new_algebra

    def __call__(self, algebra):
        if isinstance(algebra, PolynomialAlgebra):
            return PolynomialAlgebra(self.new_algebra, algebra.generators, order=algebra.order)
        if isinstance(self.new_algebra, PolynomialAlgebra):
            return PolynomialAlgebra(algebra, self.new_algebra.generators, order=self.new_algebra.order)
        if isinstance(algebra, QuotientAlgebra) and isinstance(self.new_algebra, QuotientAlgebra):
            bigger_free_algebra = PolynomialAlgebra(algebra.base_ring,
                                                    algebra.no_generators + self.new_algebra.no_generators,
                                                    variable_names=[f'x_{i}' for i in range(algebra.no_generators)] + [f'y_{i}' for i in range(self.new_algebra.no_generators)],
                                                    order=algebra.order)
            g1 = algebra.ideal.generators
            g2 = self.new_algebra.ideal.generators
            dim1 = algebra.no_generators
            dim2 = self.new_algebra.no_generators
            names1 = [f'x_{i}' for i in range(algebra.no_generators)]
            names2 = [f'y_{i}' for i in range(self.new_algebra.no_generators)]
            g = []
            for poly in g1:
                extended = list(map(lambda p: Polynomial(*list(map(lambda m: Monomial(m.coefficient, np.array(list(m.exponent_index) + [0]*dim2), names1), p.monomials))), poly.polynomials))
                g.append(extended)
            for poly in g2:
                extended = list(map(lambda p: Polynomial(*list(map(lambda m: Monomial(m.coefficient, np.array([0]*dim1 + list(m.exponent_index)), names2), p.monomials))), poly.polynomials))
                g.append(extended)
            return QuotientAlgebra(bigger_free_algebra, Ideal(g, name=str(algebra) + '(x)' + str(self.new_algebra)))
        raise NotImplementedError('Only tensor products of polynomial algebras and quotient algebras are implemented.')

    def induced_morphism(self, morphism):
        raise NotImplementedError






