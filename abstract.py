import functools
import random
from abc import ABC, abstractmethod
import base_rings
from properties import *

def convert_subrings(method):
    def wrapper(self, *args, **kwargs):
        args = list(args)
        for i, other in enumerate(args):
            if isinstance(other, int):
                args[i] = self.ring.from_canonical_subring(base_rings.ZRing()(other))
                continue
            if other.ring != self.ring and other.ring in self.ring.canonical_subrings:
                args[i] = self.ring.from_canonical_subring(other)
        return method(self, *args, **kwargs)
    return wrapper




class BaseRing(ABC):
    def __init__(self, canonical_proper_subrings=(), name=None, **properties):
        self.canonical_subrings = list(canonical_proper_subrings)
        self.name = name or 'A' + str(random.randint(0, 1000000))
        self.properties = AlgebraProperties(**properties)

    def __str__(self):
        return self.name

    def __eq__(self, other):
        return str(self) == str(other)

    @abstractmethod
    def __call__(self, *args, **kwargs):
        pass

    @abstractmethod
    def add(self, a, b):
        pass

    @abstractmethod
    def mul(self, a, b):
        pass

    @property
    @abstractmethod
    def one(self):
        pass

    @property
    @abstractmethod
    def zero(self):
        pass

    @abstractmethod
    def neg(self, a):
        pass

    @abstractmethod
    def eq(self, a, b):
        pass

    def maybe_unit_check(self, element):
        return element == self.one or element == -self.one

    def maybe_zero_check(self, element):
        return element == self.zero

    def force_zero_check(self, element):
        raise NotImplementedError("Deterministic algorithm for checking if element is zero is not implemented.")

    def force_unit_check(self, element):
        raise NotImplementedError("Deterministic algorithm for checking if element is unit is not implemented.")

    def truediv(self, a, b):
        raise NotImplementedError("Division is not implemented.")

    def longdiv(self, a, b):
        raise NotImplementedError("Long division is not implemented.")

    @abstractmethod
    def from_canonical_subring(self, element):
        pass

    def abs(self, element):
        raise NotImplementedError("Norm is not implemented.")

    def grad(self, element):
        raise NotImplementedError("Gradation is not implemented.")

    @staticmethod
    @abstractmethod
    def element_str(element):
        pass

    def info(self):
        print(self.name + ': commutative ring with 1')
        print(', '.join(filter(lambda x: self.properties[x] and x != 'exact_values', self.properties().keys())))
        print('Canonical subrings: ' + ', '.join([str(ring) for ring in self.canonical_subrings]))
        if not self.properties['exact_values']:
            print('Elements have no exact representation, but based on floating point arithmetic')
        print()


class Field(BaseRing):
    def __init__(self, canonical_proper_subrings=(), name=None,
                 exact_values=True, algebraically_closed=False, characteristics=None, normed=False):
        super().__init__(canonical_proper_subrings,
                         name,
                         field=True,
                         euclidean=True,
                         integral=True,
                         ufd=True,
                         local=True,
                         artinian=True,
                         pid=True,
                         exact_values=exact_values,
                         normal=algebraically_closed,
                         characteristics=characteristics,
                         normed=normed)

    def maybe_unit_check(self, element):
        return not self.force_zero_check(element)

    def force_unit_check(self, element):
        return not self.force_zero_check(element)

    @abstractmethod
    def truediv(self, a, b):
        pass

    def longdiv(self, a, b):
        return self.truediv(a, b), self.zero

    def info(self):
        if self.properties['characteristics'] is None:
            print(self.name + ': field')
        else:
            print(self.name + ': field of characteristic ' + str(self.properties['characteristics']))
        print(f'{"algebraically closed" if self.properties["normal"] else "not algebraically closed"}{", normed" if self.properties["normed"] else ""}')
        print('Canonical subrings: ' + ', '.join([str(ring) for ring in self.canonical_subrings]))
        if not self.properties['exact_values']:
            print('Elements have no exact representation, but based on floating point arithmetic')
        print()


class BaseElement:
    def __init__(self, ring, value):
        self.ring = ring
        self.value = value

    def __str__(self):
        return self.ring.element_str(self)

    @convert_subrings
    def __add__(self, other):
        return self.ring.add(self, other)

    @convert_subrings
    def __mul__(self, other):
        return self.ring.mul(self, other)

    @convert_subrings
    def __neg__(self):
        return self.ring.neg(self)

    @convert_subrings
    def __eq__(self, other):
        return self.ring.eq(self, other)

    def __ne__(self, other):
        return not self.ring.eq(self, other)

    @convert_subrings
    def __sub__(self, other):
        return self.ring.add(self, -other)

    @convert_subrings
    def __truediv__(self, other):
        return self.ring.truediv(self, other)

    @convert_subrings
    def __floordiv__(self, other):
        return self.ring.longdiv(self, other)[0]

    @convert_subrings
    def __mod__(self, other):
        return self.ring.longdiv(self, other)[1]

    def __pow__(self, other):
        return functools.reduce(lambda x, y: x * y, [self] * other, self.ring.one)

    def __abs__(self):
        return self.ring.abs(self)

    @functools.cached_property
    def grade(self):
        if not self.ring.properties['graded']:
            return None
        return self.ring.grade(self)

    def maybe_is_zero(self):
        return self.ring.maybe_zero_check(self)

    def maybe_is_unit(self):
        return self.ring.maybe_unit_check(self)


class GroupElement(ABC):
    def __init__(self, group, value, sequence_of_generators=None):
        self.group = group
        self.value = value
        self.sequence_of_generators = sequence_of_generators

    def __str__(self):
        return self.group.element_str(self)

    def __eq__(self, other):
        return self.group.eq(self, other)

    def __invert__(self):
        return self.group.inv(self)

    def __mul__(self, other):
        return self.group.mul(self, other)

    def __add__(self, other):
        return self.group.add(self, other)

    def __neg__(self, other):
        return self.group.neg(self, other)

    def __pow__(self, n):
        if n == 0:
            return self.group.one
        if n < 0:
            return self.group.inv(self) ** (-n)
        if n == 1:
            return self
        return self * (self ** (n - 1))


class Group(ABC):
    def __init__(self, name=None, **properties):
        self.name = name or 'G' + str(random.randint(0, 1000000))
        self.properties = properties

    def __str__(self):
        return self.name

    @abstractmethod
    def __call__(self, *args, **kwargs):
        pass

    @abstractmethod
    def mul(self, a, b):
        pass

    @property
    @abstractmethod
    def one(self):
        pass

    @abstractmethod
    def inv(self, a):
        pass

    @abstractmethod
    def eq(self, a, b):
        pass

    def add(self, a, b):
        raise NotImplementedError("Additive operations are reserved for abelian groups")

    def neg(self, a, b):
        raise NotImplementedError("Additive operations are reserved for abelian groups")

    @staticmethod
    def element_str(element):
        return str(element.value)

    def info(self):
        print(self.name + ': group')
        print(', '.join(filter(lambda x: self.properties[x], self.properties.keys())))
        print()


class FGGroup(Group, ABC):
    def __init__(self, generators, name=None, **properties):
        super().__init__(name, finitely_generated=True, **properties)
        self.generators = generators
        self.no_generators = len(generators)
        self.generator_names = [str(g) for g in generators]


class FPGroup(FGGroup, ABC):
    def __init__(self, generators, relations, name=None, **properties):
        super().__init__(generators, name, finitely_presented=True, **properties)
        self.relations = relations
        self.no_relations = len(relations)
        self.relation_names = [str(r) for r in relations]


class AbelianGroup(Group, ABC):
    def __init__(self, name=None, **properties):
        super().__init__(name, abelian=True, **properties)

    def add(self, a, b):
        return a*b

    def neg(self, a, b):
        return self.inv(a)

    @functools.cached_property
    def zero(self):
        return self.one


class FGAbelianGroup(AbelianGroup, FGGroup, ABC):
    def __init__(self, generators, name=None, **properties):
        AbelianGroup.__init__(self, name, **properties)
        FGGroup.__init__(self, generators, name, **properties)


class OpenSet(ABC):
    def __init__(self, id):
        self.id = id

    @abstractmethod
    def union(self, other):
        pass

    def __add__(self, other):
        return self.union(other)

    @abstractmethod
    def intersection(self, other):
        pass

    def __mul__(self, other):
        return self.intersection(other)


class Scheme(ABC):
    def __init__(self, affine_cover, name=None, **properties):
        self.affine_cover = affine_cover
        self.name = name or 'X' + str(random.randint(0, 1000000))
        self.properties = SchemeProperties(**properties)

    def __str__(self):
        return self.name

    @abstractmethod
    def __call__(self, i):
        pass

    @abstractmethod
    def sections_general(self, open_set):
        pass

    @abstractmethod
    def transition(self, i, j):
        pass

    @abstractmethod
    def transition_general(self, open1, open2):
        pass

    @abstractmethod
    def restriction(self, i, ij):
        pass

    @abstractmethod
    def restriction_general(self, open1, open2):
        pass

    @abstractmethod
    def global_sections(self):
        pass

    def info(self):
        print(self.name + ': scheme')
        print(', '.join(filter(lambda x: self.properties[x], self.properties.properties.keys())))
        print()











