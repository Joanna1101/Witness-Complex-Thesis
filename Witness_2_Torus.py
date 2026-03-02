## Imports
import matplotlib.pyplot as plt
import numpy as np
import gudhi as gd  
from prettytable import PrettyTable
from make_torus import *
    
def plot_witness_edges(points, landmarks, simplex_tree, name, max_dim=1):
    """
    Plots the witness complex
    
    TODO:     
    # try plotting different alphas with their persistence diagrams
    # try different shapes, frequency ratios
    # look at rips complex
    """
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    # Points... ? 
    ax.scatter(points[:,0], points[:,1], points[:,2], s=5, alpha=0.3)
    ax.scatter(landmarks[:,0], landmarks[:,1], landmarks[:,2], c="red", alpha= 0.5)

    # Edges
    for simplex in simplex_tree.get_skeleton(max_dim):
        if len(simplex[0]) == 2:  
            i, j = simplex[0]
            p, q = points[i], points[j]
            ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], 'k-', linewidth=0.5, alpha=0.3)

    ax.set_box_aspect([1,1,1])
    plt.savefig(f"Witness_Complex_{name}")
    plt.close()   
    
def witness_complex(x, y, z, name):
    """
    Constructs the witness complex for x, y, z data
    Inputs:  x, y, z......(List, List, List) torus mapping
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
    
    # Prevent clustering: use equispaced landmarks
    # use approximate period (build for periodic or semiperiodic orbit )
    landmarks = points[np.random.choice(points.shape[0], n_landmarks, replace = False)]
    witnesses = []
    for point in points:
        if point not in landmarks: witnesses.append(point)

    # plot witnesses and landmarks    
    
    # Building complex and simplex tree
    # Here, the landmarks are a subset of the witnesses
    print("Building Witness Complex")
    WC = gd.EuclideanWitnessComplex(witnesses, landmarks)
    simplex_tree = WC.create_simplex_tree(max_alpha_square = 100.0, limit_dimension=3)
    
    # Persistent Homology
    print("Plotting")
    bar_codes = simplex_tree.persistence()
    
    fig = plt.figure(figsize = (6,6))
    gd.plot_persistence_diagram(bar_codes)
    plt.savefig(f"Witness_Persistence_{name}")
    plt.close() 
    
    plot_witness_edges(points, landmarks, simplex_tree, name, max_dim = 1)
    
    return simplex_tree.dimension(), simplex_tree.num_vertices(), simplex_tree.num_simplices(), simplex_tree.betti_numbers()
 
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
    
    print("Constructing Classical Witness Complex")
    dim_c, vert_c, simplices_c, betti_c = witness_complex(x_c, y_c, z_c, "Classical_Torus.png")
    
    print("Constructing Wiggly Witness Complex")
    dim_w, vert_w, simplices_w, betti_w = witness_complex(x_w, y_w, z_w, "Wiggly_Torus.png")
    
    table = PrettyTable()
    table.field_names = ["Type", "Dimension", "Vertices", "Simplices", "Betti Numbers"]
    table.add_row(["Classical", dim_c, vert_c, simplices_c, str(betti_c)])
    table.add_row(["Wiggly", dim_w, vert_w, simplices_w, str(betti_w)])
    print(table)
    
if __name__ == "__main__":
    classical_vs_wiggly()