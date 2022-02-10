import numpy as np

def lz_complexity(s):
    '''
    Lempel-Ziv complexity as described in Kaspar and Schuster, Phys. Rev. A.
    The input iterable (see below) does not have to be binary (2-element), but
    most applications of LZ complexity have used strings of 0s and 1s.

    INPUT:
        s : string, list, or tuple, required
          sequence to calculate complexity for

    '''
    i = 0 # where we are in the current block
    k = 1 # current block length
    k_max = 1 # maximum length of matching block in history
    l = 1 # where we are in the sequence
    n = len(s)-1 # length of the sequence
    lzc = 1 # minimum lzc

    while True:
        # compare current elements with history
        if s[i+k-1] == s[l+k-1]:
            # add to the block length if the same
            k += 1 
            # check if at the end of the sequence
            if l + k >= n - 1:
                lzc += 1
                break
        else:
            # if the current block length is the largest matching block length
            # in the history, set it to kmax
            if k > k_max:
               k_max = k
            # increment what starting index in the history we are at
            i += 1
            # if we're "starting" at where we are in the sequence, increment
            # the lzc, add the maximum block length to where we are in the 
            # sequence and move on 
            if i == l:
                lzc += 1
                l += k_max
                # if we're at the end of the sequence
                if l + 1 > n:
                    break
                else:
                    i = 0
                    k = 1
                    k_max = 1
            # increment k to start the first if statement again but starting
            # at the next index in the sequence
            else:
                k = 1
    return lzc

'''
This algorithm might be improved by making sure that we don't keep checking the 
history for a bigger block when we can't get a block any bigger because there 
isn't enough history (kmax > l-i).
'''

def random_lz_complexity(n, p=0.5):
    '''
    Computes the expected Lempel-Ziv complexity for a random sequence of length
    n and expected probability of generating a 1 = p.  Useful for normalizing
    the raw lz_complexity.  This function will behave poorly if p is identically
    0 or 1.  Therefore, it would be best to estimate p from real (finite length)
    strings using pseudocounts.

    INPUT:
        n : int, required
          length of the random sequence

        p : float, optional
          probability of seeing a 1 in the sequence
    '''
    # source entropy
    h = -p*np.log2(p) - (1-p)*np.log2(1-p)
    # expected LZ complexity of binary representations of real numbers
    bn = n/np.log2(n)
    return h*bn

def nlz_complexity(s):
    '''
    Compute the LZC of a sequence as shown above normalized by the random LZC.
    '''
    n = len(s)
    p = sum(s)/len(s)
    return lz_complexity(s)/random_lz_complexity(n, p)
