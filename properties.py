import copy


class Properties:
    defaultProperties = {}

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


class AlgebraProperties(Properties):
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


class ModuleProperties(Properties):
    defaultProperties = {'vector_space': False,
                         'free': False,
                         'finite': False,
                         'finitely_generated': True,
                         'finitely_presented': True,
                         'ideal': False,
                         'cyclic': False,
                         'graded': False,
                         'zero': False}


class GroupProperties(Properties):
    defaultProperties = {'finite': False,
                         'abelian': False,
                         'cyclic': False,
                         'free': False,
                         'finitely generated': True,
                         'finitely presented': True,
                         'solvable': False,
                         'nilpotent': False,
                         'simple': False}


class SchemeProperties(Properties):
    defaultProperties = {'affine': False,
                         'projective': False,
                         'finite': False,
                         'reduced': False,
                         'irreducible': False,
                         'integral': False,
                         'normal': False,
                         'regular': False,
                         'smooth': False,
                         'complete': False,
                         'connected': False,
                         'separated': False,
                         'proper': False,
                         'quasi_projective': False,
                         'quasi_compact': True,
                         'noetherian': True}
