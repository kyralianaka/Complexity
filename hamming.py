import numpy as np
import matplotlib.pyplot as plt
import pickle


def hamming(s):
    """
    This function calculates the minimum hamming distance between the
    sequences in each pair of rows in an array. By minimum, it is meant that
    the number of differences are counted for each possible alignment of each
    pair of sequences. The output is a symmetric matrix of the minimum distances.

        INPUT: s, numpy int array of arrays

        OUTPUT: hmap, numpy float32 array of arrays

    """
    L = len(s)
    hmap = np.zeros([L, L], dtype=np.float32)

    for i in range(len(s)):

        j = i
        while j < len(s):

            counters = np.zeros(len(s[0]))

            i_startidx = 0
            j_startidx = 0  # j is the next node in s (s[i+1])

            # main loop that compares 2 sequences
            while j_startidx < len(s[i]):
                iidx = i_startidx
                jidx = j_startidx

                while iidx < len(s[i]):
                    if s[i][iidx] != s[j][jidx]:
                        counters[j_startidx] += 1
                    iidx += 1
                    jidx += 1
                    if jidx == len(s[i]):
                        jidx = 0
                j_startidx += 1

            hd = np.amin(counters) / len(s[i])
            hmap[i, j] = hd
            hmap[j, i] = hd
            j += 1

    return hmap
