from modules import *
from algebras import *


class GradedAlgebraFG(AlgebraFG):
    def __init__(self, name, base_ring, generators, gradings, order=grlex):
        super().__init__(base_ring, generators, order=order, name=name)
        self.gradings = gradings
        self.gradation_dict = {g: gradings[i] for i, g in enumerate(self.generators)}

    def grade_of_monomial(self, monomial):
        return np.dot(monomial.exponents, self.gradings)

    def is_homo(self, element):
        return all([self.grade_of_monomial(m) == self.grade_of_monomial(self.generators[0]) for m in element.value.monomials])

    def grade(self, element):
        if not self.is_homo(element):
            raise ValueError("Element is not homogeneous.")
        return self.grade_of_monomial(element.value.monomials[0])

    def max_graded_component(self, element):
        return self(max(element.value.monomials, key=self.grade_of_monomial))

    def max_grade(self, element):
        return self.grade(self.max_graded_component(element))

    def homogeneous_component(self, element, grade):
        return self(sum([m for m in element.value.monomials if self.grade_of_monomial(m) == grade], start=Monomial.create_zero(self.no_generators, self)))

    def grade_n_module(self, n):
        base = list(filter(lambda x: self.gradation_dict[x] <= n, self.generators))
        min_grade = min([self.gradation_dict[g] for g in base])
        while min_grade < n:
            base.sort(key=lambda x: self.gradation_dict[x])
            x = base.pop(0)
            base.extend([x * g for g in self.generators if self.gradation_dict[g] + self.gradation_dict[x] <= n])
            if len(base) == 0:
                break
            min_grade = min([self.gradation_dict[g] for g in base])
        base = list(set(base))
        return FreeFGModule(self.base_ring, base, name=self.name + "_" + str(n))

    def irrelevant_deal(self):
        return HomogeneousIdeal(list(filter(lambda g: self.gradation_dict[g] > 0, self.generators)), name=self.name + "_+")

    def generated_in_degree(self):
        return max(self.gradings)


class GradedFGModule(FGModule, ABC):
    def __init__(self, base_ring, no_generators, generators_gradation, name="AbstractGradedModule", generator_names=None, **properties):
        super().__init__(base_ring, no_generators, name=name, generator_names=generator_names, **properties)
        self.generators_grades = generators_gradation

    @abstractmethod
    def grade_n_part(self, n):
        pass

    @abstractmethod
    def is_homo(self, element):
        pass

    @abstractmethod
    def grade(self, element):
        pass

    def generated_in_degree(self):
        return max(self.generators_grades)


class HomogeneousIdeal(Ideal, GradedFGModule):
    def __init__(self, generators, name="", **properties):
        GradedFGModule.__init__(self, generators[0].ring, [g.degree() for g in generators], name, generators, **properties)
        Ideal.__init__(self, generators, name=name, **properties)

    def scalar_mul(self, element, scalar):
        return self([scalar * c for c in element.coefficients])

    def as_ring_element(self, element):
        return sum([c * g for c, g in zip(element.coefficients, self.generators)], start=self.ring.zero)

    def grade_n_part(self, n):
        return self([g * x for g, x in zip(self.generators, self.ring.grade_n_part(n).generators) if (g*x).grade == n])

    def is_homo(self, element):
        return self.ring.is_homo(self.as_ring_element(element))

    def grade(self, element):
        return self.ring.grade(self.as_ring_element(element))


class PolynomialGradedAlgebra(GradedAlgebraFG, PolynomialAlgebra):
    def __init__(self, base_ring, no_variables, variable_names=None, order=grlex):
        PolynomialAlgebra.__init__(self, base_ring, no_variables, variable_names=variable_names, order=order)



