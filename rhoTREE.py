import msprime
from ripser import ripser
import matplotlib.pyplot as plt
import numpy as np
from math import exp

'''Computes TREE esimate of \rho from a distance matrix'''

def Hamming_mat(tree):
    #computes Hamming distance matrix for a coalescent tree
    A = tree.genotype_matrix().transpose()
    return (2 * np.inner(A-0.5,0.5-A) + A.shape[1] / 2)

def TREErho(avg0, var0, b1):
    #print("Psi = {0}\nVar = {1}\nB1 = {2}".format(avg0, var0, b1))
    A =  5.53046629e-02
    B =  -3.74380246e-04
    C = 1.81333434e-02
    D = -1.79713403e-04
    E =  -5.93368387e-05
    y_int = 2.24756003254
    logrho = A*avg0 + B*var0 + C*b1 + D*(avg0**2) + E*(b1**2) +  y_int
    print(logrho)
    return exp(logrho)

def barstats(Rbarc):
    #computes average H0 barcode length
    A = np.array(Rbarc[0])
    blst = A[:,0]
    dlst = A[:,1]
    lengths = np.subtract(dlst,blst)[0:-1]#exclue the last term as it is inf
    if lengths != []:
        lengths = np.append(lengths,max(lengths))
    avglen = np.mean(lengths)
    varlen = np.var(lengths)
    betti1 = len(Rbarc[1])
    return avglen,varlen,betti1

#Set parameters
N = 240             #sample size
theta = 150         #mutation rate
Ne = 1000      #diploid population size
rho = 100

for N in [30,100,170,240]:
    for theta in [150,300]:
        rholst = []
        rhoTREElst = []
        rhoHlst = []   
        for rho in range(1,800):
            tree = msprime.simulate(sample_size=N, 
                                          Ne=Ne, mutation_rate = theta/(4*Ne), 
                                          recombination_rate=rho/(4*Ne))
            D = Hamming_mat(tree)
            endpoint = np.max(D)
            bars = ripser(D, distance_matrix = True)['dgms']
            avg0,var0,b1 = barstats(bars)
            
            Trho = TREErho(avg0,var0,b1)
            rholst.append(rho)
            rhoTREElst.append(Trho)
            rhoHlst.append(Trho/rho)
        
        avg = np.mean(rhoHlst)
        var = np.var(rhoHlst)
        
        plt.plot(rhoTREElst,rholst,'o', color = 'darkmagenta', alpha = 0.7, ms = 4)   
        plt.xlabel('$\\rho_T}$', size = 14)
        plt.ylabel('$\\rho$', size = 14)
        plt.title('$N = {}, \\theta = {}$'.format(N,theta))
        plt.xlim([0,900])
        plt.ylim([0,900])
        plt.plot([0,900],[0,900],'k--')
        plt.show()
        
        print('Theta = {}, N = {}'.format(theta,N))
        print('avg = {}, var = {}'.format(avg,var))
        
