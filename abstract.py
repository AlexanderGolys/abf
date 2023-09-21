from abc import ABC, abstractmethod


class BaseRing(ABC):
    field = False
    euclidean = False
    integral = False
    noetherian = True
    pid = False
    ufd = False
    local = False
    artinian = False
    dvr = False
    name = "AbstractBaseRing"
    normed = False
    graded = False

    def __init__(self, canonical_proper_subrings=()):
        self.canonical_subrings = list(canonical_proper_subrings) + [self.__class__]
        self.name = self.__class__.name

    @classmethod
    def __str__(cls):
        return cls.name

    @abstractmethod
    def __call__(self, *args, **kwargs):
        pass

    @abstractmethod
    def add(self, a, b):
        pass

    @abstractmethod
    def mul(self, a, b):
        pass

    @abstractmethod
    def one(self):
        pass

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

    def mod(self, a, b):
        raise NotImplementedError("Divisibility check is not implemented.")

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


class Field(BaseRing):
    field = True
    euclidean = True
    integral = True
    noetherian = True
    pid = True
    ufd = True
    local = True
    artinian = True
    dvr = False
    name = "AbstractField"
    normed = False
    graded = False

    def maybe_unit_check(self, element):
        return not self.force_zero_check(element)

    def force_unit_check(self, element):
        return not self.force_zero_check(element)

    @abstractmethod
    def truediv(self, a, b):
        pass

    def longdiv(self, a, b):
        return self.truediv(a, b), self.zero()
