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

def NCD(spike_array,compressor,triu_only=False):
    '''
    Generates the normalized compression distance matrix of a spike array. Assumes that C_xy =~ C_yx.

        INPUT: spike_array, numpy array of arrays; compressor; triu_only=True to get just the pairwise NCDs 

        OUTPUT: hmap, numpy array of arrays
    '''
    L = len(spike_array)
    hmap = np.zeros([L,L],dtype=np.float32)

    if compressor == 'lz':
        #Calculate each node's lz_complexity
        lzcs = np.zeros(L)
        for i in range(L):
            lzcs[i] = lzc.lz_complexity(spike_array[i,:])

        #Calculate the pairwise NCDs
        for m in range(L):
            for n in range(m+1):
                xy = np.concatenate((spike_array[m],spike_array[n]))
                hmap[n,m] = (lzc.lz_complexity(xy)-min(lzcs[m],lzcs[n]))/max(lzcs[m],lzcs[n])
    else:
        #Calculate the pairwise NCDs
        for m in range(L):
            for n in range(m+1):
                hmap[n,m] = NCD_pairwise(spike_array[m],spike_array[n],compressor)

    if triu_only==False:
        #Mirror array over the diagonal 
        i_lower = np.tril_indices(L)
        hmap[i_lower] = np.transpose(hmap)[i_lower]  
    else:
        #Get the upper triangle NCD values only
        idxs = np.triu_indices(L,k=1)
        hmap = hmap[idxs]

    return hmap

def NCD_pairwise(s1,s2, compressor):
    '''
    Calculates the normalized compression distance between two
    sequences using the compressor specified. 

        INPUT:
            s1, s2: numpy arrays, any int dtype acceptable but uint8 recommended
            compressor: gzip, bz2, 'gzip_padded'(pads x2), snappy, 'ppm', 'ppmc', 'lz','nlz'
        OUTPUT:
            NCD: numpy float

    Note: all compressed lengths are a value in bytes (except ppmc which returns bits),
    but NCD is a ratio, so a conversion to bits would be pointless. Compressed lengths
    in bits can be found using the function clen below.
    '''
    if compressor in [gzip, bz2, 'gzip_padded']:
        #compressor.compress returns a bytes object, and len(bytesObject) gives the number of bytes
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
        #snappy.compress() also seems to return a bytes object, so len() returns the number of bytes of the compressed sequence
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
        #sequences are converted to byte streams and then compressed into other byte streams
        stream1 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s1))
        stream2 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s2))
        stream12 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s12))
        stream21 = ppm_compress_mod.compress(ppm_compress_mod.makestream(s21))
        #the length of the byte stream/number of bytes in the compressed argument
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
        xobj = ppmc.ppmc(x)
        yobj = ppmc.ppmc(y)
        xyobj = ppmc.ppmc(xy)
        yxobj = ppmc.ppmc(yx)
        C_x = round(xobj.compress()) 
        C_y = round(yobj.compress())
        C_xy = round(xyobj.compress())
        C_yx = round(yxobj.compress())
        NCD = (np.amin([C_xy,C_yx])-np.amin([C_x, C_y]))/np.amax([C_x, C_y])
    if compressor == 'lz':
        C_x = lzc.lz_complexity(s1)
        C_y = lzc.lz_complexity(s2)
        C_xy = lzc.lz_complexity(np.concatenate((s1,s2)))
        C_yx = lzc.lz_complexity(np.concatenate((s2,s1)))
        NCD = (np.amin([C_xy,C_yx])-np.amin([C_x, C_y]))/np.amax([C_x, C_y])
    if compressor == 'nlz':
        xy = np.concatenate((s1,s2))
        yx = np.concatenate((s2,s1))
        C_len = np.zeros(4)
        for i,seq in enumerate([s1,s2,xy,yx]):
            t = len(seq)
            p = float(sum(seq))/t
            h = -p*np.log2(p) - (1-p)*np.log2(1-p)
            c_st = h*t/np.log2(t)
            C_len[i] = lzc.lz_complexity(seq)/c_st
        NCD = (np.amin(C_len[2:])-np.amin(C_len[:2]))/np.amax(C_len[:2])
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

if __name__ == "__main__":
    s1 = np.tile(np.random.binomial(1,0.5,size=10),200)
    s2 = np.tile(np.random.binomial(1,0.5,size=10),200)
    val = NCD_pairwise(s1,s2,'ppm')
    val2 = NCD_pairwise(s1,s2,gzip)
    print("NCD(x,y): " + str(val))

