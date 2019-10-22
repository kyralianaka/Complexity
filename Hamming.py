import numpy as np
import matplotlib.pyplot as plt
import pickle

'''
This function will calculate the minimum hamming distance between each of the nodes in a 
network and produce a heat map array. This array can be used to design a figure for 
illustrating the existence of 3 different network states.

    INPUT: s, spike array, numpy int array of arrays
    
    OUTPUT: hmap, numpy float32 array of arrays
    
'''

def Hamming(s):
    
    #create square array to append hds into
    L = len(s)
    hmap = np.zeros([L,L],dtype=np.float32) #Make sure making floats default
    
    for i in range(len(s)): #for every node in spike array s
        
        j = i
        while j < len(s):
        
            counters = np.zeros(len(s[0])) 
            
            "The main loop that compares 2 sequences"
            i_startidx = 0
            j_startidx = 0 #j is the next node in s (s[i+1])
            
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
                    
            hd = np.amin(counters)/len(s[i])
            hmap[i,j] = hd
            hmap[j,i] = hd
            j += 1
    
    return hmap 

if __name__ == '__main__':
    
    
    pickle_in = open('toynet','rb')
    toy_net = pickle.load(pickle_in)

    
    hm = Hamming(toy_net)
    im = plt.imshow(hm, cmap='hot', interpolation='none')
    plt.colorbar(im)


