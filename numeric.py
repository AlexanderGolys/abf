import math
import scipy.special


class Binomials:
    @staticmethod
    def binomial_list_n(k, max_n):
        return [scipy.special.comb(n, k, exact=True) for n in range(max_n + 1)]

    @staticmethod
    def binomial_list_nk(max_n, max_k):
        return [[scipy.special.comb(n, k, exact=True) for k in range(max_k + 1)] for n in range(max_n + 1)]

    @staticmethod
    def binomial_list_k(n):
        return [scipy.special.comb(n, k, exact=True) for k in range(n+1)]

