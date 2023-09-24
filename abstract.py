import copy
import functools
from abc import ABC, abstractmethod, abstractproperty


class AlgebraProperties:
    defaultProperties = {'field': False,
                         'euclidean': False,
                         'integral': False,
                         'noetherian': True,
                         'pid': False,
                         'ufd': False,
                         'local': False,
                         'artinian': False,
                         'dvr': False,
                         'normed': False,
                         'graded': False,
                         'exact_values': True,
                         'normal': False}

    def __init__(self, **kwargs):
        self.properties = copy.deepcopy(self.defaultProperties)
        for key, value in kwargs.items():
            self.properties[key] = value

    def __getitem__(self, key):
        return self.properties[key]

    def __setitem__(self, key, value):
        self.properties[key] = value

    def __call__(self):
        return self.properties


class ModuleProperties:
    defaultProperties = {'vector space': False,
                         'free': False,
                         'finite': False,
                         'finitely generated': True,
                         'finitely presented': True,
                         'ideal': False,
                         'cyclic': False,
                         'graded': False}

    def __init__(self, **kwargs):
        self.properties = copy.deepcopy(self.defaultProperties)
        for key, value in kwargs.items():
            self.properties[key] = value

    def __getitem__(self, key):
        return self.properties[key]

    def __setitem__(self, key, value):
        self.properties[key] = value

    def __call__(self):
        return self.properties


class BaseRing(ABC):

    def __init__(self, canonical_proper_subrings=(), name="AbstractBaseRing", **properties):
        self.canonical_subrings = list(canonical_proper_subrings)
        self.name = name
        self.properties = AlgebraProperties(**properties)

    def __str__(self):
        return self.name


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
        return element == self.one() or element == -self.one()

    def maybe_zero_check(self, element):
        return element == self.zero()

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
    def __init__(self, canonical_proper_subrings=(), name="AbstractField",
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
        return self.truediv(a, b), self.zero()

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

    def __add__(self, other):
        return self.ring.add(self, other)

    def __mul__(self, other):
        return self.ring.mul(self, other)

    def __neg__(self):
        return self.ring.neg(self)

    def __eq__(self, other):
        return self.ring.eq(self, other)

    def __ne__(self, other):
        return not self.ring.eq(self, other)

    def __sub__(self, other):
        return self.ring.add(self, -other)

    def __truediv__(self, other):
        return self.ring.truediv(self, other)

    def __floordiv__(self, other):
        return self.ring.longdiv(self, other)

    def __mod__(self, other):
        return self.ring.longdiv(self, other)[1]

    def __pow__(self, other):
        return functools.reduce(lambda x, y: x * y, [self] * other)

    def __abs__(self):
        return self.ring.abs(self)

    @functools.cached_property
    def grad(self):
        if not self.ring.properties['graded']:
            return None
        return self.ring.grad(self)

    def maybe_is_zero(self):
        return self.ring.maybe_zero_check(self)

    def maybe_is_unit(self):
        return self.ring.maybe_unit_check(self)


class FGModule(ABC):
    def __init__(self, base_ring, generators, name="AbstractBaseRing", **properties):
        assert isinstance(base_ring, BaseRing)
        assert len(generators) > 0

        self.ring = base_ring
        self.name = name
        self.properties = ModuleProperties(**properties)
        self.generators = generators
        self.no_generators = len(generators)
        if self.no_generators == 1:
            self.properties['cyclic'] = True
        if self.ring.properties['field']:
            self.properties['vector space'] = True

    def __call__(self, coefficients):
        return ModuleElement(self, coefficients)

    @property
    def zero(self):
        return self([self.ring.zero] * self.no_generators)

    def __str__(self):
        return self.name

    def add(self, a, b):
        return self([a + b for a, b in zip(a.coefficients, b.coefficients)])

    @abstractmethod
    def scalar_mul(self, element, scalar):
        pass

    def grad(self, element):
        if not self.properties['graded']:
            return None
        raise NotImplementedError("Gradation is not implemented.")

    def to_ring_element(self):
        raise NotImplementedError("Conversion to ring element is possible only for ideals.")

    def info(self):
        print(self.name + ': finitely generated ' + str(self.ring) + '-module')
        print(str(self.no_generators) + f' generator{"s" if self.no_generators > 1 else ""}')
        print('<' + ', '.join([str(g) for g in self.generators]) + '>')
        print(', '.join(filter(lambda x: self.properties[x], self.properties().keys())))
        print()



class ModuleElement:
    def __init__(self, module, coefficients):
        self.ring = coefficients[0].ring
        self.coefficients = coefficients
        self.module = module

    def __str__(self):
        if self.module.properties['ideal']:
            return str(self.module.as_ring_element(self))
        return ' + '.join([str(c) + '*' + str(g) for c, g in zip(self.coefficients, self.module.generators)])

    def __add__(self, other):
        return self.module.add(self, other)

    def __mul__(self, scalar):
        assert isinstance(scalar, BaseElement)
        assert scalar.ring == self.ring

        return self.module.scalar_mul(self, scalar)

    def __neg__(self):
        return self.module.scalar_mul(self, -self.ring.one())

    def __eq__(self, other):
        return self.coefficients == other.coefficients

    def __ne__(self, other):
        return not self.coefficients == other.coefficients

    def __sub__(self, other):
        return self.module.add(self, -other)

    @functools.cached_property
    def grad(self):
        if not self.module.properties['graded']:
            return None
        return self.module.grad(self)

    def maybe_is_zero(self):
        return all([self.ring.maybe_zero_check(c) for c in self.coefficients])

    def force_is_zero(self):
        return all([self.ring.force_zero_check(c) for c in self.coefficients])

    def __call__(self):
        return self.module.to_ring_element(self)

