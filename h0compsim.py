import msprime
from ripser import ripser
import matplotlib.pyplot as plt
import numpy as np

'''Simulates coalescent process for a range of rho, plots
 rho vs avg barcode length and rho vs Betti1'''

def Hamming_mat(tree):
    #computes Hamming distance matrix for a coalescent tree
    A = tree.genotype_matrix().transpose()
    return (2 * np.inner(A-0.5,0.5-A) + A.shape[1] / 2)

def avg_h0length(Rbarc):
    #computes average H0 barcode length
    A = np.array(Rbarc[0])
    blst = A[:,0]
    dlst = A[:,1]
    lengths = np.subtract(dlst,blst)[0:-1] #exclue the last term as it is inf
    avglen = np.average(lengths)
    return avglen

#Set parameters
N = 100             #sample size
theta = 500         #mutation rate
Ne = 10000          #diploid population size

rholst = []         #used to store values
b1lst = []
avbarlst = []

for rho in range(1000):
    rholst.append(rho)
    
    #simulate a coalescent tree for given parameters and a range of rho
    tree_sequence = msprime.simulate( 
    sample_size=N, Ne=Ne, mutation_rate = theta/(4*Ne), recombination_rate=rho/(4*Ne))
    D = Hamming_mat(tree_sequence)
    dgms = ripser(D, distance_matrix = True)['dgms'] #computes PH with ripser
    
    b1lst.append(len(dgms[1]))
    avbarlst.append(avg_h0length(dgms))
    

#plots rho against num of H1 features.
plt.plot(b1lst,rholst,'o', color = 'darkmagenta', alpha = 0.7, ms = 4)
plt.xlabel('$\\beta_1$', size = 14)
plt.ylabel('$\\rho$', size = 14)
plt.show()
    
#plots rho against avg length of H0 barcodes
plt.plot(avbarlst,rholst,'o', color = 'darkmagenta', alpha = 0.7, ms = 4)   
plt.xlabel('$\\psi$', size = 14)
plt.ylabel('$\\rho$', size = 14)
plt.show()
