import unittest
from graded import *


class TestAll(unittest.TestCase):
    def test_primitive_algebras(self):
        Z, Q, R, C = ZRing(), QField(), RFloating(), CFloating()
        self.assertEqual(Z(2) + Z(-3), Z(-1))
        self.assertEqual(Z(2)*Z(-3), Z(-6))
        a = Q(2)/Z(-3)
        self.assertEqual(Q(2)/Z(-3), Q(2, -3))
        self.assertEqual(R(2)/Q(-3), R(-2/3))
        self.assertEqual(C(0, 1)**2, C(-1))

        for ring in [Z, Q, R, C]:
            self.assertEqual(ring.zero, ring(0))
            self.assertEqual(ring.one, ring(1))
            self.assertEqual(-ring.one, ring(-1))
            self.assertEqual(ring.one + ring.zero, ring.one)
            self.assertTrue(ring.one.maybe_is_unit())
            self.assertFalse(ring.zero.maybe_is_unit())
            self.assertEqual((ring.one*4)/(ring.one*2), ring.one*2)

        self.assertEqual(Q(3)/Q(2), Q(3, 2))
        self.assertEqual(Q(3)//Q(2), Q(3, 2))
        self.assertEqual(Q(3) % Q(2), Q.zero)
        self.assertEqual(Z(5) % Z(2), Z(1))
        self.assertEqual(Z(5) // Z(2), Z(2))
        self.assertEqual(Z(5)**0, Z(1))
        self.assertEqual(abs(Z(-5)), Z(5))
        self.assertEqual(Z.one.grade, None)

    def test_polynomials(self):
        pass

    def test_algebras(self):
        pass

    def test_modules(self):
        pass

    def test_graded(self):
        pass

    def test_matrices(self):
        Q = QField()
        M, N = Matrix.fill_free_coefficients(Q, 2, 2)
        commutator = M@N - N@M
        commuting_ideal = commutator.vanishing_ideal()
        commuting_algebra = commutator.vanishing_algebra()

        self.assertEqual(commuting_ideal.ring.no_generators, 8)
        self.assertEqual(commuting_algebra.no_relations, 4)
        print(commuting_ideal)
        print(commuting_algebra.generator_elements[0])
        print(commuting_algebra)




