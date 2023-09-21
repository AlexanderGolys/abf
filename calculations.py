from polynomials import *
from formatting import *
import math


if __name__ == '__main__':
    dual_numbers = FPAlgebra(C, [Generator('ε')], [PolynomialTensorForm([0, 0, 1], 'ε')], 'D', C.one, C.zero)
    dual_numbers.show_info()
    x, y, z, t = Generator('x'), Generator('y'), Generator('z'), Generator('t')
    two_lines = FPAlgebra(C, [x, y], [PolynomialTensorForm([[0, 0], [0, 1]], [x, y])], 'A', C.one, C.zero)
    two_lines.show_info()