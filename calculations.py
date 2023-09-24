from polynomials import *
from formatting import *
from algebras import *
from base_rings import *
import math
from modules import *

def test_algebras():
    # base rings
    Z, Q, R, C = ZRing(), QField(), RFloating(), CFloating()

    # polynomial rings

    Zx = PolynomialAlgebra(Z, 1)
    f = Zx(Polynomial(Monomial(Z(3), [2]), Monomial(Z(3), [2]), Monomial(Z(1), [1]), Monomial(Z(1), [0])))

    print(f)
    print(f**2)
    print()
    Zx.info()


def test_modules():
    Z, Q = ZRing(), QField()
    I = Ideal([Z(2)], name='(2)')
    two = I([Z.one])
    print(two)
    print(f'{two} + {two} = {two + two}')
    print(f'3 * {two} = {two * Z(3)}')
    I.info()


def test_groebner():
    Q = QField()
    P = KPolynomialAlgebra(Q, 2)
    P.info()
    # x, y, z, w = P(Polynomial(Monomial(Q.one, [1, 0, 0, 0]))), P(Polynomial(Monomial(Q.one, [0, 1, 0, 0]))), P(Polynomial(Monomial(Q.one, [0, 0, 1, 0]))), P(Polynomial(Monomial(Q.one, [0, 0, 0, 1])))
    x, y, z = P(Polynomial(Monomial(Q.one, [1, 0, 0]))), P(Polynomial(Monomial(Q.one, [0, 1, 0]))), P(Polynomial(Monomial(Q.one, [0, 0, 1])))
    I = KPolynomialIdeal([x**2-y, x**3-z], name='I')
    I.info()
    # J = I.reduce_basis()
    # J.info()


if __name__ == '__main__':
    test_groebner()
