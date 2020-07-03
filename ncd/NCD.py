import gzip, bz2
import snappy
from complexity.ncd import ppm_compress_mod, ppmc
from complexity.ncd import arithmeticcoding, ppmmodel
from complexity.lzc import lzc
import contextlib, sys, os
import scipy.cluster.hierarchy as hier
import numpy as np
import matplotlib.pyplot as plt
import pylab

def NCD(spike_array,compressor):
    '''
    Generates the normalized compression distance matrix of a spike array. Assumes that C_xy =~ C_yx.

        INPUT: spike_array, numpy array of arrays; compressor

        OUTPUT: hmap, numpy array of arrays
    '''
    L = len(spike_array)
    hmap = np.zeros([L,L],dtype=np.float32)

    #Calculate each node's lz_complexity
    lzcs = np.zeros(L)
    for i in range(L):
        lzcs[i] = lzc.lz_complexity(spike_array[i,:])

    #Calculate the pairwise NCDs
    for m in range(L):
        for n in range(m+1):
            xy = np.concatenate((spike_array[m],spike_array[n]))
            hmap[n,m] = (lzc.lz_complexity(xy)-min(lzcs[m],lzcs[n]))/max(lzcs[m],lzcs[n])

    #Mirror array over the diagonal 
    i_lower = np.tril_indices(L)
    hmap[i_lower] = np.transpose(hmap)[i_lower]  # make the matrix symmetric

    return hmap

def NCD_pairwise(s1,s2, compressor):
    '''
    Calculates the normalized compression distance between two
    sequences using the compressor specified. 

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
    return NCD

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

def NCD_clusters(NCD_mat):
    #Get upper triangle as list
    dmat_condensed = NCD_mat[np.triu_indices(NCD_mat.shape[0],k=1)]
    #Linkage matrix
    link_mat = hier.linkage(dmat_condensed,method='average')
    #Compute the dendrogram
    dendro = hier.dendrogram(link_mat,no_plot=True)
    #Row ordering according to the dendrogram
    leaves = dendro['leaves']
    #Create ordered matrix
    new_mat = NCD_mat[leaves,:]
    new_mat = new_mat[:,leaves]   
    return new_mat


    
 

