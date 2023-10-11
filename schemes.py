from arrows import *
from graded import *


class OpenOfAffine(OpenSet):
    def __init__(self, distinguished_union):
        super().__init__(set(distinguished_union))
        self.ring = list(distinguished_union)[0].ring
        self.distinguished_list = list(self.id)
        self.is_distinguished = len(self.id) == 1

    def union(self, other):
        return OpenOfAffine(self.id.union(other.id))

    def intersection(self, other):
        if self.is_distinguished:
            return OpenOfAffine([self.distinguished_list[0]*u for u in other.id])
        return self.distinguished_list[0]*other + self.distinguished_list[1:]*other


class DistinguishedOpen(OpenSet):
    def __init__(self, element):
        super().__init__(element)
        self.ring = element.ring

    def union(self, other):
        OpenOfAffine([self, other])

    def intersection(self, other):
        return DistinguishedOpen(self.id*other.id)

    def to_open_of_affine(self):
        return OpenOfAffine([self.id])


class AffineScheme(Scheme):
    def __init__(self, algebra, name=None, **properties):
        if name is None:
            name = 'Spec(' + str(algebra) + ')'
        super().__init__([0], name,
                         affine=True,
                         noetherian=algebra.properties['noetherian'],
                         finite=algebra.properties['artinian'],
                         integral=algebra.properties['integral'],
                         reduced=algebra.properties['integral'],
                         irreducible=algebra.properties['integral'],
                         separated=True,
                         quasi_projective=True,
                         connected=algebra.properties['integral'],
                         normal=algebra.properties['normal'],
                         **properties)
        self.algebra = algebra

    def __call__(self, i):
        if i != 0:
            raise ValueError("affine cover of Spec(R) has only one element with id 0.")
        return self.algebra

    def restriction(self, i, ij):
        if i != 0:
            raise ValueError("affine cover of Spec(R) has only one element with id 0.")
        if ij != 0:
            raise ValueError("affine cover of Spec(R) has only one element with id 0.")
        return AlgebraMorphism.identity(self.algebra)

    def sections_general(self, open_set):
        if not isinstance(open_set, DistinguishedOpen):
            if open_set.is_distinguished:
                return self.sections_general(DistinguishedOpen(open_set.distinguished_list[0]))
            raise ValueError("Only sections of distinguished open sets are implemented")
        return self.algebra.localisation(open_set.id)

    def transition_general(self, open1, open2):
        if not isinstance(open1, DistinguishedOpen):
            if open1.is_distinguished:
                return self.transition_general(DistinguishedOpen(open1.distinguished_list[0]), open2)
        if not isinstance(open2, DistinguishedOpen):
            if open2.is_distinguished:
                return self.transition_general(open1, DistinguishedOpen(open2.distinguished_list[0]))
        return AlgebraMorphism.identity(self(open1*open2))

    def restriction_general(self, open1, open2):
        return self(open1).localisation_morphism_at_element((open1*open2).id)

    def transition(self, open_set1, open_set2):
        return AlgebraMorphism.identity(self.algebra)

    def global_sections(self):
        return self.algebra


class Spec(ContravariantFunctor):
    def __init__(self, base_ring):
        super().__init__(Algebra, Scheme, name='Spec(-)')
        self.base_ring = base_ring

    def __call__(self, algebra):
        return AffineScheme(algebra, name='Spec(' + str(algebra) + ')')

    def induced_morphism(self, morphism):
        pass  # TODO: implement scheme morphisms


class AffineSpace(AffineScheme):
    def __init__(self, ring, n):
        algebra = PolynomialAlgebra(ring, n)
        super().__init__(algebra, name='A^' + str(n))


class ProjectiveSpace(Scheme):
    def __init__(self, algebra, n, **properties):
        name = 'P^' + str(n)
        self.algebra = algebra
        self.n = n
        super().__init__(list(range(n+1)), name,
                         affine=False,
                         projective=True,
                         noetherian=algebra.properties['noetherian'],
                         finite=algebra.properties['artinian'],
                         integral=algebra.properties['integral'],
                         reduced=algebra.properties['integral'],
                         irreducible=algebra.properties['integral'],
                         separated=True,
                         quasi_projective=True,
                         connected=algebra.properties['integral'],
                         normal=algebra.properties['normal'],
                         proper=True,
                         **properties)
        self.algebra = algebra

    def __call__(self, i):
        return PolynomialAlgebra(self.algebra, self.n, [f'x_{j}/x_{i}' for j in range(self.n) if j != i])

    def sections_general(self, open_set):
        raise NotImplementedError("general sections of projective schemes are not implemented")

    #  ------ TODO -------
    def transition(self, i, j):
        pass

    def transition_general(self, open1, open2):
        pass

    def restriction(self, i, ij):
        pass

    def restriction_general(self, open1, open2):
        pass

    def global_sections(self):
        pass





