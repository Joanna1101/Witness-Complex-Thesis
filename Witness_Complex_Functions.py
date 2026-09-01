## Imports ###################################################################################################################################
##############################################################################################################################################
import matplotlib.pyplot as plt
import numpy as np
import gudhi as gd  
from prettytable import PrettyTable
from make_torus import *
import itertools



    
## PLOTTING ##################################################################################################################################
##############################################################################################################################################
def plot_witness_edges(witnesses, landmarks, simplex_tree, name, max_dim=1):
    """
    Plots the witness complex
    
    TODO: Progression of epsilons: as epsilon increases more things get connected
    """
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    ax.scatter(witnesses[:,0], witnesses[:,1], witnesses[:,2], c = "blue", s=5, alpha=0.3, label = "Witnesses")
    ax.scatter(landmarks[:,0], landmarks[:,1], landmarks[:,2], c = "red", s=10, alpha= 0.5, label = "Landmarks")

    # Edges BETWEEN LANDMARKS
    for simplex in simplex_tree.get_skeleton(max_dim):
        if len(simplex[0]) == 2:  
            i, j = simplex[0]
            p, q = landmarks[i], landmarks[j] 
            ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], 'k-', linewidth=0.5, alpha=0.3)

    ax.set_box_aspect([1,1,1])
    ax.legend(loc="upper right")
    plt.savefig(f"deSilva_Witness_Complex_{name}")
    plt.close(fig)   




## DESILVA_CARLSSON  #########################################################################################################################
##############################################################################################################################################
def desilva_carlsson_simplex_tree(landmarks, witnesses, max_dim = 3, k=0):
    """
    Builds a simplex tree with connections between landmarks
    Inputs:  landmarks......(np array)
             witnesses......(np array)
             max_dim........(int) max simplex dimension
             k..............(int) optional parameter to relax complex
    Outputs: st.............(simplex tree)
    """
    L = len(landmarks)
    W = len(witnesses)
    
    # Distance Matrix: row = landmark, col = witness
    diffs = landmarks[:, None, :]- witnesses[None, :, :]
    D = np.linalg.norm(diffs, axis = 2)    # norms equivalent
    
    st = gd.SimplexTree()
    
    # Landmarks are 0-simplices, numbered 0...L
    for landmark in range (L):
        st.insert([landmark], filtration=0.0)
        
    # Add simplices according to smallest distances in D
    for i in range(W):
        ordered_dists = np.argsort(D[:, i])
        
        # p = num vertices in simplex
        for p in range(2, max_dim + 2):
            # Just in case
            if p + k > L:
                print(f"Epsilon too large: {p+k} > number of landmarks")
                return []
            
            # p + epsilon smallest distances
            allowed = set(ordered_dists[:p+k])
            
            # vertices in the simplex
            verts = list(ordered_dists[:p+k])
            
            if set(verts).issubset(allowed):                
                # Filtration value for simplex 'verts' witnessed by i 
                # is the max distance between a vertex in the simplex and i
                filt = float(np.max(D[verts, i]))
                st.insert(verts, filtration = filt)

    st.initialize_filtration()
    return st
    
    
def desilva_carlsson_witness_complex(x, y, z, name, n_landmarks = 50, max_dim = 3, k = 0):
    """
    Constructs the witness complex for x, y, z data
    Inputs:  x, y, z......(List, List, List) torus mapping
             alpha........(float) max distance parameter
             name.........(String) plot names
             n_landmarks..(Int) number of landmarks (equispaced in time)
             max_dim......(Int) max simplex dimension
             k............(Int) optional parameter to relax complex
    Outputs: dim..........(Int) Dimension of simplex tree
             vertices.....(Int) Vertices in simplex tree
             simplices....(Int) Simplices in simplex tree
             betti........(List) Betti numbers of simplex tree
             plot.........(Figures) persistence diagram and plot of complex
    """
    # Formatting data
    points = np.vstack([x,y,z]).T
    N = len(points)
    
    # Landmarks equispaced in time
    landmark_is = np.linspace(0, N-1, n_landmarks, dtype = int)
    landmarks = points[landmark_is]
    mask = np.ones(N, dtype=bool)
    mask[landmark_is] = False
    witnesses = points[mask]
    
    # Building complex and simplex tree
    print("Building deSilva - Carlsson Witness Complex")
    simplex_tree = desilva_carlsson_simplex_tree(landmarks, witnesses, max_dim, k)
    
    # Persistent Homology
    print("Plotting")
    bar_codes = simplex_tree.persistence()
    
    fig = plt.figure(figsize = (6,6))
    gd.plot_persistence_diagram(bar_codes)
    plt.savefig(f"Witness_Persistence_{name}")
    plt.close(fig) 
    
    plot_witness_edges(witnesses, landmarks, simplex_tree, name, max_dim = 1)
    
    return (simplex_tree,
            simplex_tree.dimension(),
            simplex_tree.num_vertices(),
            simplex_tree.num_simplices(),
            simplex_tree.betti_numbers(),
            landmarks)
    
    
    
    
    
## FUZZY #####################################################################################################################################
##############################################################################################################################################
