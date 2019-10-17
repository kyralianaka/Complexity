# -*- coding: utf-8 -*-
"""
Created on Tue Jul 30 13:47:27 2019

@author: kyra
"""

import gzip, bz2
import snappy
from complexity.ncd import ppm_compress_mod, ppmc
from complexity.lzc import lzc
import contextlib, sys, os
from complexity.ncd import arithmeticcoding, ppmmodel
import numpy as np
import matplotlib.pyplot as plt
import pylab

'''
This function will calculate the normalized compression distance matrix
of a spike array which can be plotted to analyze structure in a network.

    INPUT: spike array, numpy int array of arrays
           compressor: gzip, 'gzip_padded', bz2, 'ppm', or snappy

    OUTPUT: hmap, numpy float32 array of arrays

'''

def NCD(spike_array,compressor):
    '''
    Generates the normalized compression distance matrix of a spike array.

        INPUT: spike_array, numpy array of arrays; compressor
            Note: The best compressor to be used in this function with
            spike array data has been found to be gzip, and the compressor
            achieves most idempotency when each sequence is cast as np.uint8s
            and padded such that the length of the sequences are close to the
            window size (32768 bits).

        OUTPUT: hmap, numpy array of arrays
    '''
    L = len(spike_array)
    hmap = np.zeros([L,L],dtype=np.float32)

    for m in range(len(spike_array)):
        for n in range(len(spike_array)):
            hmap[n,m] = NCD_pairwise(spike_array[m],spike_array[n],compressor)

    return hmap

def NCD_pairwise(s1,s2, compressor):
    '''
    Calculates the normalized compression distance between two
    sequences using the compressor specified. Currently, the gzip compressor is
    considered the most accurate due to testing of idempotency. Padding the sequence
    using the gzip compressor may also help.

        INPUT:
            s1, s2: numpy arrays, any int dtype acceptable but uint8 recommended
            compressor: gzip, bz2, 'gzip_padded'(pads x2), snappy, 'ppm', 'ppmc', 'lz'
        OUTPUT:
            NCD: numpy float

    Note: all compressed lengths are a value in bytes (except ppmc which returns bits),
    but NCD is a ratio, so a conversion to bits would be pointless. Compressed lengths
    in bits can be found using the function clen below.
    '''
    if compressor in [gzip, bz2, 'gzip_padded']:
        if compressor == 'gzip_padded':
            s1 = np.tile(s1,3)
            s2 = np.tile(s2,3)
            compressor = gzip
        C_x = len(compressor.compress(s1))
        C_y = len(compressor.compress(s2))
        xy = np.concatenate((s1,s2))
        yx = np.concatenate((s2,s1))
        C_xy = len(compressor.compress(xy))
        C_yx = len(compressor.compress(yx))
        NCD = (np.amin([C_xy,C_yx])-np.amin([C_x, C_y]))/np.amax([C_x, C_y])
    if compressor == snappy:
        s1_str = ''.join(map(str,s1))
        s2_str = ''.join(map(str,s2))
        C_x = len(compressor.compress(s1_str))
        C_y = len(compressor.compress(s2_str))
        xy = ''.join(map(str, np.concatenate((s1,s2))))
        yx = ''.join(map(str,np.concatenate((s2,s1))))
        C_xy = len(compressor.compress(xy))
        C_yx = len(compressor.compress(yx))
        NCD = (np.amin([C_xy,C_yx])-np.amin([C_x, C_y]))/np.amax([C_x, C_y])
    if compressor == 'ppm':
        s12 = np.concatenate((s1,s2))
        s21 = np.concatenate((s2,s1))
        stream1 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s1))
        stream2 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s2))
        stream12 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s12))
        stream21 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s21))
        C_x = len(stream1.getbuffer())
        C_y = len(stream2.getbuffer())
        C_xy = len(stream12.getbuffer())
        C_yx = len(stream21.getbuffer())
        NCD = (np.amin([C_xy,C_yx])-np.amin([C_x, C_y]))/np.amax([C_x, C_y])
    if compressor == 'ppmc':
        x = ''.join(map(str,s1))
        y = ''.join(map(str,s2))
        xy = ''.join(map(str, np.concatenate((s1,s2))))
        yx = ''.join(map(str,np.concatenate((s2,s1))))
        xobj, yobj = ppmc.ppmc(x), ppmc.ppmc(y)
        xyobj, yxobj = ppmc.ppmc(xy), ppmc.ppmc(yx)
        C_x, C_y = xobj.compress(), yobj.compress()
        C_xy, C_yx = xyobj.compress(), yxobj.compress()
        NCD = (np.amin([C_xy,C_yx])-np.amin([C_x, C_y]))/np.amax([C_x, C_y])
    if compressor == 'lz':
        C_x = lzc.lz_complexity(s1)
        C_y = lzc.lz_complexity(s2)
        C_xy = lzc.lz_complexity(np.concatenate((s1,s2)))
        C_yx = lzc.lz_complexity(np.concatenate((s2,s1)))
        NCD = (np.amin([C_xy,C_yx])-np.amin([C_x, C_y]))/np.amax([C_x, C_y])
        '''
        In order to test the least modified form of pppm from the Nayuki project.

        seqs = dict.fromkeys(['x', 'y', 'xy', 'yx'])
        C = dict.fromkeys(seqs)
        seqs['x'] = s1.astype(np.uint8)
        seqs['y'] = s2.astype(np.uint8)
        seqs['xy'] = np.concatenate((seqs['x'],seqs['y']))
        seqs['yx'] = np.concatenate((seqs['y'],seqs['x']))
        for i in seqs:
            stringi = ''.join(map(str, seqs[i]))
            filein = open("filein.txt","w")
            filein.write(stringi)
            filein.close()
            ppm_compress_mod2.main("filein.txt","fileout.txt")
            statinfo1 = os.stat("fileout.txt")
            C[i] = statinfo1.st_size
        NCD = (np.amin([C['xy'],C['yx']])-np.amin([C['x'],C['y']]))/np.amax([C['x'],C['y']])
        '''
    return NCD

#x = np.random.binomial(1,0.5,200)
#score = NCD_pairwise(x,x,'ppmc')

def clen(x, compressor):
    '''
    Determines the compressed length of a sequence in bits using the compressor
    specified.
    '''
    if compressor in [gzip, bz2, snappy]:
        clen = len(compressor.compress(x))*8 #to get the clen in bits of the bytes object
    if compressor == snappy:
        x = ''.join(map(str,x))
        clen = len(compressor.compress(x))*8 #to get the clen in bits of the bytes object
    if compressor == 'ppm':
        x_stream = ppm_compress_mod.compress(ppm_compress_mod.makestream(x))
        clen = len(x_stream.getbuffer())*8 #to get the clen in bits from the byte stream
    if compressor == 'ppmc':
        x = ''.join(map(str,x))
        xobj = ppmc.ppmc(x)
        clen = xobj.compress()
    return clen


if __name__ == '__main__':
    '''
    pylab_pretty_plot()

    pickle_in = open('toynet','rb')
    toy_net = pickle.load(pickle_in)
    pickle_in.close()

    hm = NCD(toy_net)
    im = plt.imshow(hm, cmap='hot', interpolation='none')
    r = np.arange(8)
    plt.xticks(r, r.astype(str))
    plt.colorbar(im)

    x = np.random.binomial(1,0.5,1000)
    xstr = ''.join(map(str,x))
    C_x = len(snappy.compress(xstr))
    y = np.random.binomial(1,0.5,1000)
    ystr = ''.join(map(str,y))
    C_y = len(snappy.compress(ystr))

    xy = np.concatenate((x,y))
    xy = ''.join(map(str,xy))
    C_xy = len(snappy.compress(xy))
    yx = np.concatenate((y,x))
    yx = ''.join(map(str,yx))
    C_yx = len(snappy.compress(yx))

    xx = np.concatenate((x,x))
    xx = ''.join(map(str,xx))
    C_xx = len(snappy.compress(xx))

    yy = np.concatenate((y,y))
    yy = ''.join(map(str,yy))
    C_yy = len(snappy.compress(yy))

    NCD_xy = (np.amin([C_xy,C_yx])-np.amin([C_x,C_y]))/np.amax([C_x,C_y])
    NCD_xx = (C_xx-C_x)/C_x
    NCD_yy = (C_yy-C_y)/C_y
    print(NCD_xy, NCD_xx, NCD_yy)

    x = np.random.binomial(1,0.5,1000)
    lz_x = NCD_pairwise(x,x,'lz')
    y = np.tile([1,0,1,0],250)
    lz_y = NCD_pairwise(y,y,'lz')
    lz_xy = NCD_pairwise(x,y,'lz')
    lz_yx = NCD_pairwise(y,x,'lz')
    print(lz_x, lz_y, lz_xy, lz_yx)
    '''
