import msprime
from ripser import ripser
import numpy as np
import gudhi as gd

'''simulates a coalescent process, computes persistenh homology and plots'''

#Set parameters
N = 100             #sample size
theta = 500        #mutation rate
rho = 70          #recombnation rate
Ne = 10000          #diploid population size

#simulate a coalescent tree for given parameters
tree_sequence = msprime.simulate( 
    sample_size=N, Ne=Ne, mutation_rate = theta/(4*Ne), recombination_rate=rho/(4*Ne))


def Hamming_mat(tree):
    #computes Hamming distance matrix for a coalescent tree
    A = tree.genotype_matrix().transpose()
    return (2 * np.inner(A-0.5,0.5-A) + A.shape[1] / 2)

def Rips_to_GUDHI(barcode): 
    #converts H1 elements of ripser output to a GUDHI compatible output
    Gbarc = []
    for i in range(len(barcode[0])):        
        Gbarc.append((0,(barcode[0][i][0],barcode[0][i][1])))
        
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

# print(tree_sequence.get_num_mutations())

D = Hamming_mat(tree_sequence)
dgms = ripser(D, distance_matrix = True)['dgms']#computes PH with ripser
barcode = Rips_to_GUDHI(dgms)

#plot barcode and persistence plot with GUDHI
gd.plot_persistence_diagram(barcode, colormap = ('navy','darkmagenta','fuchsia'))    
gd.plot_persistence_barcode(barcode, colormap = ('navy','darkmagenta','fuchsia'))

