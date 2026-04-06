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
    plt.savefig(f"CORRECT_Witness_Complex_{name}")
    plt.close()   
    
    
    
    
## COMPLEX BUILDING ##########################################################################################################################
##############################################################################################################################################
def build_simplex_tree(landmarks, witnesses, max_alpha_square, max_dim = 3, k_nearest = None):
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
    
    
def witness_complex(x, y, z, alpha, name):
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
    
    # 5 nearest neighbors seems reasonable?
    simplex_tree = build_simplex_tree(landmarks, witnesses, alpha, 2, 5)
    
    # Persistent Homology
    print("Plotting")
    bar_codes = simplex_tree.persistence()
    
    fig = plt.figure(figsize = (6,6))
    gd.plot_persistence_diagram(bar_codes)
    plt.savefig(f"Witness_Persistence_{name}")
    plt.close() 
    
    plot_witness_edges(witnesses, landmarks, simplex_tree, name, max_dim = 1)
    
    return simplex_tree.dimension(), simplex_tree.num_vertices(), simplex_tree.num_simplices(), simplex_tree.betti_numbers()
 
 
 
 
## TESTS #####################################################################################################################################
##############################################################################################################################################
def classical_vs_wiggly():
    # Keep these the same
    R = 15
    r = 2
    w1 = 1
    p = (1+np.sqrt(5))/2
    q = 1
    
    # 10,000 looks like a nice orbit but it makes computer explodes
    # TODO: how does number of simplices increase with number of landmarks? exp? !?
    n = 300
    
    # Wiggle parameters
    # Try alpha >= 1 and note betti numbers
    a = 3
    k = 5
    
    x_c, y_c, z_c, t_c = make_torus(R, r, w1, p, q, 10e-3, n, "classical", True)
    x_w, y_w, z_w, t_w = make_torus(R, r, w1, p, q, 10e-3, n, "wiggly", True, a, k)
    
    alpha = 100.0
    
    print("Constructing Classical Witness Complex")
    dim_c, vert_c, simplices_c, betti_c = witness_complex(x_c, y_c, z_c, alpha, "Classical_Torus.png")
    
    print("Constructing Wiggly Witness Complex")
    dim_w, vert_w, simplices_w, betti_w = witness_complex(x_w, y_w, z_w, alpha, "Wiggly_Torus.png")
    
    table = PrettyTable()
    table.field_names = ["Type", "Dimension", "Vertices", "Simplices", "Betti Numbers"]
    table.add_row(["Classical", dim_c, vert_c, simplices_c, str(betti_c)])
    table.add_row(["Wiggly", dim_w, vert_w, simplices_w, str(betti_w)])
    print(table)
    
def test_alphas():
    alphas = [5.0, 20.0, 100.0]
    
    # Torus parameters
    R = 15
    r = 2
    w1 = 1
    p = (1+np.sqrt(5))/2
    q = 1
    n = 300  # Number of points
    
    x_c, y_c, z_c, t_c = make_torus(R, r, w1, p, q, 10e-3, n, "classical", True)
    
    table = PrettyTable()
    table.field_names = ["alpha", "dim", "vertices", "simplices", "betti"]
    for alpha in alphas:
        print(f"Constructing Complex For alpha = {alpha}")
        dim, vert, simplices, betti = witness_complex(x_c, y_c, z_c, alpha, f"{alpha}_Classical_Torus.png")
        table.add_row([alpha, dim, vert, simplices, betti])

    print(table)

    
## MAIN ######################################################################################################################################
##############################################################################################################################################
if __name__ == "__main__":
    # classical_vs_wiggly()
    test_alphas()