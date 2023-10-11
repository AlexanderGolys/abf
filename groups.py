from modules import *

class FiniteGroup(FPGroup):
    def __init__(self, generators, relations, name=None, **properties):
        super().__init__(generators, relations, name=name, **properties)

    def __call__(self, sequence_of_generators):
        sequence_of_generators = self.reduce_sequence(sequence_of_generators)
        x = GroupElement(self, sequence_of_generators, sequence_of_generators)
        return self.reduce(x)

    @staticmethod
    def reduce_sequence(sequence):
        if len(sequence) < 2:
            return sequence
        for i, a, b in enumerate(zip(sequence[:-1], sequence[1:])):
            x, y = a[0], b[0]
            if x == y:
                sequence = sequence[:i] + [(x, a[1] + b[1])] + sequence[i + 2:]
                return FiniteGroup.reduce_sequence(sequence)

    def reduce(self, element):
        length = len(element.sequence_of_generators)
        for relation in self.relations:
            rel_len = len(relation)
            for i in range(length - rel_len + 1):
                if element.sequence_of_generators[i:i + rel_len] == relation:
                    element.sequence_of_generators = element.sequence_of_generators[:i] + element.sequence_of_generators[i + rel_len:]
                    return self.reduce(element)

        return element

    @functools.lru_cache(maxsize=1000)
    def mul(self, a, b):
        return self(a.sequence_of_generators + b.sequence_of_generators)

    def eq(self, g, h):
        return g.sequence_of_generators == h.sequence_of_generators

    def generator_elements(self):
        return [self([g]) for g in self.generators]

    @functools.cached_property
    def elements(self):
        def add_multiples(list_of_elements):
            new_list = []
            for g in list_of_elements:
                for h in self.generators:
                    if g*h not in new_list:
                        new_list.append(g*h)
            return new_list

        list_of_elements = self.generator_elements()
        while True:
            new_list = add_multiples(list_of_elements)
            if len(new_list) == len(list_of_elements):
                return list_of_elements
            list_of_elements = new_list

    @functools.cached_property
    def table(self):
        return {(g, h): g*h for g, h in itertools.product(self.elements, repeat=2)}

    @functools.cached_property
    def order(self):
        return len(self.elements)

    @property
    def one(self):
        return self([])

    @functools.lru_cache(maxsize=1000)
    def inv(self, a):
        multiples = [(g, self.table[(a, g)]) for g in self.elements]
        for g, h in multiples:
            if h == self.one:
                return g

    @staticmethod
    def element_str(element):
        return "*".join(element.sequence_of_generators)


class FreeGroup(FPGroup):
    def __call__(self, sequence_of_generators):
        x = GroupElement(self, sequence_of_generators, sequence_of_generators)
        return x

    def mul(self, a, b):
        return self(a.sequence_of_generators + b.sequence_of_generators)

    def eq(self, g, h):
        return g.sequence_of_generators == h.sequence_of_generators

    @property
    def one(self):
        return self([])

    def inv(self, a):
        return self(list(map(lambda s: (s[0], -s[1]), a.sequence_of_generators[::-1])))


class FreeAbelianGroup(FGAbelianGroup, FPGroup):
    def __init__(self, generators, name=None, **properties):
        AbelianGroup.__init__(self, name, **properties)
        FPGroup.__init__(self, generators, [], name=name, abelian=True, **properties)
        self.generators.sort(key=lambda x: str(x))

    def __call__(self, sequence_of_generators):
        if isinstance(sequence_of_generators[0], tuple):
            appeared = set(map(lambda x: x[0], sequence_of_generators))
            rest = set(self.generators) - appeared
            sequence_of_generators.extend([(g, 0) for g in rest])
            seq = sorted(sequence_of_generators, key=lambda g: str(g[0]))
            x = GroupElement(self, np.array(list(map(lambda x: x[1], seq))).astype(int), seq)
            return x
        else:
            return GroupElement(self, np.array(sequence_of_generators).astype(int), [(g, n) for g, n in zip(self.generators, sequence_of_generators) if n != 0])

    def mul(self, a, b):
        return self(a.value + b.value)

    def eq(self, g, h):
        return g.sequence_of_generators == h.sequence_of_generators

    @property
    def one(self):
        return self([])

    def inv(self, a):
        return self(list(map(lambda s: (s[0], -s[1]), a.sequence_of_generators[::-1])))








