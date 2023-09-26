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



