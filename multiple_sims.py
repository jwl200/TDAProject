import msprime
from ripser import ripser
import matplotlib.pyplot as plt
import numpy as np

'''simmulates multiple coalecent processes for same variables 
and plots distribution of various H1 barcode statistics'''

def Hamming_mat(tree):
    #computes Hamming distance matrix for a coalescent tree
    A = tree.genotype_matrix().transpose()
    return (2 * np.inner(A-0.5,0.5-A) + A.shape[1] / 2)

def Rips_to_GUDHI(barcode): 
    #converts H1 elements of ripser output to a GUDHI compatible output
    Gbarc = []
    for i in range(len(barcode[1])):        
        Gbarc.append((1,(barcode[1][i][0],barcode[1][i][1])))
    return Gbarc

def descriptives(Gbarc): 
    #finds the number of features, birth times, death times, 
    #and lifespans for a set of GUDHI barcodes
    num = len(Gbarc)
    birth,death,length = [],[],[]
    
    for i in Gbarc:
        b = i[1][0]
        d = i[1][1]
        birth.append(b)
        death.append(d)
        length.append(d-b)
        
    return num, birth, death, length

def M_sims(N,theta,rho,Ne,m):
    #simulate for given parameters m times and save descriptives
    Nlst,Blst,Dlst,Llst = [],[],[],[]
    
    for i in range(m): 
        tree_sequence = msprime.simulate( 
            sample_size=N, Ne=Ne, mutation_rate = theta/(4*Ne), recombination_rate=rho/(4*Ne))
        
        D = Hamming_mat(tree_sequence)
        dgms = ripser(D, distance_matrix = True)['dgms']
        barcode = Rips_to_GUDHI(dgms)
        Num,B,D,L = descriptives(barcode)
        Nlst.append(Num)
        Dlst.append(D)
        Blst.append(B)
        Llst.append(L)
        
    #reduce outputs to 1D lists to make them pyplot compatible        
    Blst = [item for sublist in Blst for item in sublist] 
    Dlst = [item for sublist in Dlst for item in sublist]
    Llst = [item for sublist in Llst for item in sublist]
    
    return Nlst,Blst,Dlst,Llst

#Set parameters
N = 100            #sample size
theta = 500        #mutation rate
rho1 = 20          #recombnation rate first simulation
rho2 = 40          #recombnation rate second simulation
rho3 = 80          #recombnation rate third simulation
rho4 = 160         #recombnation rate fourth simulation
Ne = 10000         #diploid population size
m = 1000           #number of simulations


#simulate for each rho m times
N1,B1,D1,L1 = M_sims(N, theta, rho1, Ne, m)
N2,B2,D2,L2 = M_sims(N, theta, rho2, Ne, m)
N3,B3,D3,L3 = M_sims(N, theta, rho3, Ne, m)
N4,B4,D4,L4 = M_sims(N, theta, rho4, Ne, m)


#plot histograms of descriptors for each simulation
plt.figure(figsize=(8,6))
plt.hist(B1, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 20', color='navy', density = True)
plt.hist(B2, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 40', color='violet', density = True)
plt.hist(B3, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 80', color='darkmagenta', density = True)
plt.hist(B4, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 160', color='deeppink', density = True)

plt.title("Birth Times")
plt.legend(loc='upper right', fontsize = 14)
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
plt.show()

plt.figure(figsize=(8,6))
plt.hist(D1, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 20', color='navy', density = True)
plt.hist(D2, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 40', color='violet', density = True)
plt.hist(D3, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 80', color='darkmagenta', density = True)
plt.hist(D4, bins = 75, range = (0,1000), histtype='step', label = '$\\rho$ = 160', color='deeppink', density = True)

plt.title("Death Times")
plt.legend(loc='upper right', fontsize = 14)
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
plt.show()

plt.figure(figsize=(8,6))
plt.hist(N1, bins = 75, range = (0,80), histtype='step', label = '$\\rho$ = 20', color='navy', density = True)
plt.hist(N2, bins = 75, range = (0,80), histtype='step', label = '$\\rho$ = 40', color='violet', density = True)
plt.hist(N3, bins = 75, range = (0,80), histtype='step', label = '$\\rho$ = 80', color='darkmagenta', density = True)
plt.hist(N4, bins = 75, range = (0,80), histtype='step', label = '$\\rho$ = 160', color='deeppink', density = True)

plt.title("Number of Features")
plt.legend(loc='upper right', fontsize = 14)
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
plt.show()

plt.figure(figsize=(8,6))
plt.hist(L1, bins = 75, histtype='step',range=(0,250), label = '$\\rho$ = 20', color='navy', density = True)
plt.hist(L2, bins = 75, histtype='step',range=(0,250), label = '$\\rho$ = 40', color='violet', density = True)
plt.hist(L3, bins = 75, histtype='step',range=(0,250), label = '$\\rho$ = 80', color='darkmagenta', density = True)
plt.hist(L4, bins = 75, histtype='step',range=(0,250), label = '$\\rho$ = 160', color='deeppink', density = True)

plt.title("Lifespans")
plt.legend(loc='upper right', fontsize = 14)
ax = plt.gca()
ax.axes.yaxis.set_visible(False)
plt.show()
