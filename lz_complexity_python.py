# -*- coding: utf-8 -*-
"""
Created on Wed Jun 26 11:39:16 2019

"""

def lz_complexity(s):
    '''
    Lempel-Ziv complexity as described in Kaspar and Schuster, Phys. Rev. A.
    The input iterable (see below) does not have to be binary (2-element), but
    most applications of LZ complexity have used strings of 0s and 1s.

    INPUT:
        s : string, list, or tuple, required
          sequence to calculate complexity for

    '''
    i, k, l = 0, 1, 1
    k_max = 1
    n = len(s)-1
    lzc = 1
    while True:
        if s[i+k-1] == s[l+k-1]:
            k += 1
            if l + k >= n - 1:
                lzc += 1
                break
        else:
            if k > k_max:
               k_max = k
            i += 1
            if i == l:
                lzc += 1
                l += k_max
                if l + 1 > n:
                    break
                else:
                    i = 0
                    k = 1
                    k_max = 1
            else:
                k = 1
    return lzc