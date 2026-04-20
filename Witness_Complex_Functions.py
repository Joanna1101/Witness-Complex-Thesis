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
    plt.close()   
    
## ALPHA??? BUILDING #########################################################################################################################
##############################################################################################################################################    
def alpha_simplex_tree(landmarks, witnesses, max_alpha_square, max_dim = 3, k_nearest = None):
    """
    Builds a simplex tree with connections between landmarks not witnesses
    Inputs:  landmarks...........(np array)
             witnesses...........(np array)
             max_alpha_square....(int) distance threshold
             max_dim.............(int) max simplex dimension
             k_nearest...........(bool/int) use up to k nearest landmarks
    Outputs: st..................(simplex tree)
    
    
    Just use k nearest? Just use alpha? Epsilon?? How many landmarks necessary to get homology correct? 
    
    distance matrix of landmarks and witnesses (euclidean)
    landmark gets witnessed if [insert method here]
    linked list of landmarks and their associated witness
    if 2 landmarks share a witness they get connected in the simplex tree 
    
    draw triangle if three edges = clique, easiest 
    only draw triangle if three landmarks share a single witness btwn them
    tada witness complex 
    """
    L = len(landmarks)
    st = gd.SimplexTree()

    # Landmarks first
    for i in range(L):
        st.insert([i], filtration=0.0)

    for w in witnesses:
        # distances from this witness to all landmarks
        diffs = landmarks - w
        
        # euclidean distance
        # IOU? other metrics here? 
        d2 = np.linalg.norm(diffs, axis=1)**2

        # Sort by distance
        order = np.argsort(d2)
        if k_nearest is not None:
            order = order[:k_nearest]

        # Distance threshold
        close = [i for i in order if d2[i] <= max_alpha_square]
        
        # If less than 2 witnesses dont witness the landmark
        if len(close) < 2:
            continue  

        # Simplices for landmarks
        for dim in range(1, max_dim + 1):
            for comb in itertools.combinations(close, dim + 1):
                filt = max(d2[list(comb)])  # Filter based on max squared distance from witness -> verticies
                st.insert(list(comb), filtration=filt)

    st.initialize_filtration()
    return st
    
    
def alpha_witness_complex(x, y, z, alpha, name):
    """
    Constructs the witness complex for x, y, z data
    Inputs:  x, y, z......(List, List, List) torus mapping
             alpha........(float) max distance parameter
             name.........(String) plot names
    Outputs: dim..........(Int) Dimension of simplex tree
             vertices.....(Int) Vertices in simplex tree
             simplices....(Int) Simplices in simplex tree
             betti........(List) Betti numbers of simplex tree
             plot.........(Figures) persistence diagram and plot of complex
    """
    # Formatting data
    n_landmarks = 50
    points = np.vstack([x,y,z]).T
    
    # Equispaced landmarks (in time) (?)
    N = len(points)
    landmark_is = np.linspace(0, N-1, n_landmarks, dtype = int)
    landmarks = points[landmark_is]
    mask = np.ones(N, dtype=bool)
    mask[landmark_is] = False
    witnesses = points[mask]
    
    # Building complex and simplex tree
    print("Building Witness Complex")
    # WC = gd.EuclideanWitnessComplex(witnesses, landmarks)
    
    # APPARENTLY, the simplex tree builds edges between witnesses 
    # GUDHI does this because it makes nearest neighbor landmark queries easier
    # simplex_tree = WC.create_simplex_tree(max_alpha_square = 100.0, limit_dimension=3)
    
    # 5 nearest neighbors seems reasonable? idk what im doing here
    simplex_tree = alpha_simplex_tree(landmarks, witnesses, alpha, 2, 5)
    
    # Persistent Homology
    print("Plotting")
    bar_codes = simplex_tree.persistence()
    
    fig = plt.figure(figsize = (6,6))
    gd.plot_persistence_diagram(bar_codes)
    plt.savefig(f"Witness_Persistence_{name}")
    plt.close() 
    
    plot_witness_edges(witnesses, landmarks, simplex_tree, name, max_dim = 1)
    
    return simplex_tree.dimension(), simplex_tree.num_vertices(), simplex_tree.num_simplices(), simplex_tree.betti_numbers()





## DESILVA_CARLSSON BUILDING #################################################################################################################
##############################################################################################################################################
def desilva_carlsson_simplex_tree(landmarks, witnesses, max_dim = 3,):
    """
    Builds a simplex tree with connections between landmarks not witnesses
    Inputs:  landmarks...........(np array)
             witnesses...........(np array)
             max_dim.............(int) max simplex dimension
    Outputs: st..................(simplex tree)
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
            
            # vertices are p smallest distances
            verts = list(ordered_dists[:p])
            
            # Filtration value for simplex 'verts' witnessed by i 
            # is the max distance between a vertex in the simplex and i
            # different metric maybe 
            filt = float(np.max(D[verts, i]))
            st.insert(verts, filtration = filt)

    st.initialize_filtration()
    return st
    
    
def desilva_carlsson_witness_complex(x, y, z, name, n_landmarks = 50, max_dim = 3):
    """
    Constructs the witness complex for x, y, z data
    Inputs:  x, y, z......(List, List, List) torus mapping
             alpha........(float) max distance parameter
             name.........(String) plot names
             n_landmarks..(Int) number of landmarks (equispaced in time)
             max_dim......(Int) max simplex dimension
    Outputs: dim..........(Int) Dimension of simplex tree
             vertices.....(Int) Vertices in simplex tree
             simplices....(Int) Simplices in simplex tree
             betti........(List) Betti numbers of simplex tree
             plot.........(Figures) persistence diagram and plot of complex
    """
    # Formatting data
    n_landmarks = 50
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
    simplex_tree = desilva_carlsson_simplex_tree(landmarks, witnesses, max_dim)
    
    # Persistent Homology
    print("Plotting")
    bar_codes = simplex_tree.persistence()
    
    fig = plt.figure(figsize = (6,6))
    gd.plot_persistence_diagram(bar_codes)
    plt.savefig(f"Witness_Persistence_{name}")
    plt.close() 
    
    plot_witness_edges(witnesses, landmarks, simplex_tree, name, max_dim = 1)
    
    return (simplex_tree.dimension(),
            simplex_tree.num_vertices(),
            simplex_tree.num_simplices(),
            simplex_tree.betti_numbers())