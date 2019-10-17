# -*- coding: utf-8 -*-
"""
Created on Sun Mar  3 19:54:35 2019

@author: Kyra Kadhim
"""
from numpy import unique,array
from collections import defaultdict
from .utilities import *

def etc(s,ent_est='MM'):
    '''
    This is code adapted from the MATLAB script written by Nithin Nagaraj to calculate Effort-To-Compress complexity.
    This method utilizes Non-Sequential Recursive Pair Substitution (NSRPS), described on pg 5 of
    "Dynamical Complexity of Short and Noisy Time Series" published in 2016. ETC is claimed to be the number of
    iterations necessary for the entropy of the sequence to reach 0. This function currently only works for binary
    sequences.

    INPUT:
        s: string, list, or array, required
            should consist of only 1's and 0's

    OUTPUT:
        H : float
             entropy of each step of the algorithm
        N : integer
             ETC measure, number of iterations until the entropy reaches 0 in the NSRPS algorithm

    '''
    # convert to an array if a string was input
    if type(s) == str:
        s = array([int(c) for c in s])
    # The main loop for NSRPS iteration
    N = 0 # The ETC measure 'N'
    H_new = natstobits(entropyd(s,est=ent_est))
    L = len(s)
    ent_list = [H_new*L]
    # loop until entropy of sequence is zero
    while (H_new.all() > 1e-6) and (L > 1):
        # finds the pair of symbols with the maximum frequency
        pair = find_pair(s)
        # substitutes the pair with a new symbol
        s_new = substitute(s,pair)
        # Shannon Entropy of the new sequence
        H_new = natstobits(entropyd(s_new,est=ent_est))
        # Shannon Entropy of the new sequence
        L = len(s_new) # length of the new sequence
        ent_list.append(H_new*L) # adds the new entropy value to the list
        N = N + 1 # ETC measure incremented by 1
        s = s_new
    return ent_list,N


def find_pair(s):
    "Finds the most frequent pair of symbols in the sequence"
    L = len(s)
    indx = 0
    pairs = {}
    while indx < (L-1):
        a = s[indx]
        b = s[indx+1]
        pair = (a,b)
        if pair in pairs:
            pairs[pair] = pairs[pair] + 1
        else:
            pairs[pair] = 1
        if a == b:
            if indx < L-2:
                if s[indx+2] == a:
                    indx = indx + 1
        indx = indx + 1
    freq_pair = max(pairs, key=pairs.get)
    return freq_pair


def substitute(s,freq_pair):
    "Pair substitution step of NSRPS"
    s_new = []
    L = len(s)
    alphabet = unique(s)
    M = len(alphabet)
    I = max(alphabet)
    rep_sym = I+1 # new symbol
    indx = 0
    while (indx<L-1):
        a = s[indx]
        b = s[indx+1]
        if a==freq_pair[0] and b==freq_pair[1]:
            s_new.append(rep_sym)
            indx = indx+1
        else:
            s_new.append(a)
        indx = indx+1
        if indx == L-1:
            s_new.append(s[indx])
    return s_new
