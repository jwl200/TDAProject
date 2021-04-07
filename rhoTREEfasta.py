''' Calculating \rho_T for FASTA data'''

from ripser import ripser
import gudhi as gd
import numpy as np
from math import exp
from fastaprocessing import format_data,single_line_fasta,empty_matrix,populate_matrix

def TREErho(avg0, var0, b1):
    #print("Psi = {0}\nVar = {1}\nB1 = {2}".format(avg0, var0, b1))
    A =  5.53046629e-02
    B =  -3.74380246e-04
    C = 1.81333434e-02
    D = -1.79713403e-04
    E =  -5.93368387e-05
    y_int = 2.24756003254
    logrho = A*avg0 + B*var0 + C*b1 + D*(avg0**2) + E*(b1**2) +  y_int
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

def Rips_to_GUDHI(barcode): 
    #converts H1 elements of ripser output to a GUDHI compatible output
    Gbarc = []
    for i in range(len(barcode[0])):        
        Gbarc.append((0,(barcode[0][i][0],barcode[0][i][1])))
    for i in range(len(barcode[1])):        
        Gbarc.append((1,(barcode[1][i][0],barcode[1][i][1])))
    return Gbarc

def Compute_rhoT(fasta_input,empty_fasta, single_line = False): #fasta file and an empty fasta file
    print("Reading FASTA file...")
    if not single_line:
        single_line_fasta(fasta_input, empty_fasta)
        fasta_file = empty_fasta
    else:
        fasta_file = fasta_input
        
    lines = format_data(fasta_file,'fasta')
    print("Computing Hamming distance matrix...")
    #print(lines[0:5])
    matrix = empty_matrix(lines)
    hamm_matrix = populate_matrix(matrix, lines)
    print('Computing barcodes ...')
    dgms = ripser(hamm_matrix,  distance_matrix=True, maxdim=1)['dgms']
    Rbarc = dgms 
    Gbarc = Rips_to_GUDHI(Rbarc)
    print('Plotting persistence diagrams')
    gd.plot_persistence_diagram(Gbarc, colormap = ('navy','darkmagenta','fuchsia'))    
    gd.plot_persistence_barcode(Gbarc, colormap = ('navy','darkmagenta','fuchsia'))
    print('Computing barcode statistics...')
    avg0, var0, b1 = barstats(Rbarc)
    print('Computing estimate of rho')
    rhoT = TREErho(avg0, var0, b1)
    print('$\psi = {}, \Phi = {}, \\beta_1 = {}, \\rho_T = {}$'.format(avg0, var0, b1, rhoT))
    return hamm_matrix,Rbarc

fasta_input = r'C:\Users\jake-\OneDrive\Documents\Uni Year 3\Year 3 Project\Python files\datasets\covid\covid500.fasta'
empty_fasta = r'C:\Users\jake-\OneDrive\Documents\Uni Year 3\Year 3 Project\Python files\datasets\covid\covid500oneline.fasta'
hamm_mat, Rbarc = Compute_rhoT(fasta_input, empty_fasta, single_line = False)
