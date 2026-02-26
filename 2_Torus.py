## Imports
import matplotlib.pyplot as plt
import sympy
import numpy as np
import gudhi as gd  

def make_torus(R, r, w1, w2, t):
    """
    Makes a 2-torus, returns x,y,z embedding in R3. 
    
    :param R: Radius of large loop
    :param r: Radius of small loop
    :param w1: angular velocity 1
    # :param w2: angular velocity 2
    :param t: sequence of times
    """
    theta1 = w1*t # Big angle
    theta2 = w2*t # Little angle

    # Flow on torus mapped to R2
    x = (R+r*np.cos(theta2))*np.cos(theta1)
    y = (R+r*np.cos(theta2))*np.sin(theta1)
    z = r*np.sin(theta2)
    
    return x,y,z

def plot_2_torus():
    # Parameters
    R = 10
    r = 5
    w1 = 1
    
    # Testing two different angular velocities
    # w2_1 = np.pi
    w2_1 = 3
    # t_1 = np.linspace(0, 200*w2_1, 1000)
    t_1 = np.linspace(0, 3*np.pi*1, 1000)
    
    # Why is this one smoother?
    # w2_2 = (1 + np.sqrt(5))/2
    w2_2 = 2
    # t_2 = np.linspace(0, 200*w2_2, 1000)
    t_2 = np.linspace(0, 2*np.pi*1, 1000)
    
    x1,y1,z1 = make_torus(R, r, w1, w2_1, t_1)
    x2,y2,z2 = make_torus(R,r,w1, w2_2, t_2)
    
    
    # Visualization
    fig, axes = plt.subplots(1,2, figsize = (12,5), subplot_kw={'projection': '3d'})
    ax1, ax2 = axes
    
    ax1.set_xlabel('x')
    ax1.set_ylabel('y')
    ax1.set_zlabel('z')
    ax1.set_title(r'$\omega_2 = $ ' +f'{w2_1:.4f}')
    ax2.set_xlabel('x')
    ax2.set_ylabel('y')
    ax2.set_zlabel('z')
    ax2.set_title(r'$\omega_2 = $ ' +f'{w2_2:.4f}')
    
    # Rainbow Torus
    colors = (t_1 - t_1.min()) / (t_1.max() - t_1.min())
    ax1.scatter(x1, y1, z1, c=colors, cmap='cool', s=1)
    
    colors = (t_2 - t_2.min()) / (t_2.max() - t_2.min())
    ax2.scatter(x2, y2, z2, c=colors, cmap='cool', s=1)

    ax1.set_box_aspect([1,1,1])
    ax2.set_box_aspect([1,1,1])
    plt.savefig("Torus_Knot.png")
    plt.show()
    
def plot_witness_edges(points, simplex_tree, max_dim=1):
    # try plotting different alphas with their persistence diagrams
    # try different shapes, frequency ratios
    # look at rips complex
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    # Points
    ax.scatter(points[:,0], points[:,1], points[:,2], s=5, alpha=0.3)

    # Edges
    for simplex in simplex_tree.get_skeleton(max_dim):
        if len(simplex[0]) == 2:  
            i, j = simplex[0]
            p, q = points[i], points[j]
            ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], 'k-', linewidth=0.5, alpha=0.3)

    ax.set_box_aspect([1,1,1])
    plt.savefig("Euclidean_Witness_1.png")
    plt.show()   
    
def witness_complex():
    # Torus parameters
    R = 10
    r = 5
    w1 = 1
    w2 = (1 + np.sqrt(5))/2
    t = np.linspace(0, 200*w2, 1000)
    x,y,z = make_torus(R, r, w1, w2, t)
    
    # Formatting data
    n_landmarks = 200
    points = np.vstack([x,y,z]).T
    landmarks = points[np.random.choice(points.shape[0], n_landmarks, replace = False)]
    
    # Building complex and simplex tree
    # Here, the landmarks are a subset of the witnesses
    WC = gd.EuclideanWitnessComplex(points, landmarks)
    simplex_tree = WC.create_simplex_tree(max_alpha_square = 100.0, limit_dimension=3)
    
    # Computing Persistent Homology
    bar_codes = simplex_tree.persistence()
    gd.plot_persistence_diagram(bar_codes)
    
    print(f"Dimension: {simplex_tree.dimension()}")
    print(f"Number of Vertices: {simplex_tree.num_vertices()}")
    print(f"Number of Simplices: {simplex_tree.num_simplices()}")
    print(f"Betti Numbers: {simplex_tree.betti_numbers()}")
    plot_witness_edges(points, simplex_tree, max_dim = 1)
 
    
if __name__ == "__main__":
    # plot_2_torus()
    witness_complex()