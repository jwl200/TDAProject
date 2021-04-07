import msprime
from ripser import ripser
import matplotlib.pyplot as plt
import numpy as np
from math import log

'''computes camara et al's estimate of rho for a coalescent 
process with 40 seggregating sites'''

def Hamming_mat(tree):
    #computes Hamming distance matrix for a coalescent tree
    A = tree.genotype_matrix().transpose()
    return (2 * np.inner(A-0.5,0.5-A) + A.shape[1] / 2)

#Set parameters
N = 100             #sample size
theta = 7.75         #mutation rate
Ne = 10000      #diploid population size

f1 = -0.08225
f2 = 3697124
g1 = 0.11483
g2 = 8.4455

f = f1*((N/f2)*log(N)-N)
g = g1*N+g2

rholst = []
rhoPHlst = []
rhoHlst = [] 
a = 0           #store values of rho and predicted rho and rho/rhoPH

for rho in range(1,250):
    for i in range(40):
        tree_sequence = msprime.simulate( 
            sample_size=N, Ne=Ne, mutation_rate = theta/(4*Ne), recombination_rate=rho/(4*Ne))
        x = tree_sequence.get_num_mutations()
        #if 39 <= x <= 41:
        if x == 40:    
            #a = a+1
            D = Hamming_mat(tree_sequence)
            dgms = ripser(D, distance_matrix = True)['dgms']
            b1 = len(dgms[1])
            rhoPH = g * ( (1 + (1/f) ) ** b1 - 1) #compute predicted rho
            rholst.append(rho)
            rhoPHlst.append(rhoPH)
            rhoHlst.append(rhoPH/rho)

avg = np.mean(rhoHlst)
var = np.var(rhoHlst)

plt.plot(rhoPHlst,rholst,'o', color = 'darkmagenta', alpha = 0.7, ms = 4)
plt.xlabel('$\\rho_{PH}$', size = 14)
plt.ylabel('$\\rho$', size = 14)
plt.ylim([0,300])
plt.xlim([0,350])
plt.title('$s = 40$')
#plt.title('$39 \leq s \leq 41$')
plt.plot([0,300],[0,300],'k--')
plt.show()

# print(avg, var)