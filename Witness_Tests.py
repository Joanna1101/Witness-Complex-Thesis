## Imports ###################################################################################################################################
##############################################################################################################################################
from Witness_Complex_Functions import *
from make_torus import *

## CIRCLE #####################################################################################################################################
##############################################################################################################################################
def circle(N, r):
    theta = np.linspace(0, 2*np.pi, N, endpoint = False)
    x = r*np.cos(theta)
    y = r*np.sin(theta)
    z = np.zeros_like(x)
    return x, y, z

def plot_circle_points(x, y, landmarks_idx, name, simplex_tree = None, plot_edges = False):
    plt.figure(figsize=(6,6))

    # Landmarks in red
    plt.scatter(x[landmarks_idx], y[landmarks_idx], c='gold', s=10, label='Landmarks')

    # Witnesses in blue
    mask = np.ones(len(x), dtype=bool)
    mask[landmarks_idx] = False
    plt.scatter(x[mask], y[mask], c='mediumpurple', s=5, label='Witnesses')
    plt.title(f"Landmarks and Witnesses")
    
    if plot_edges:
        # Find the edges
        edges = []
        for simplex, filt in simplex_tree.get_filtration():
            if len(simplex) == 2:
                edges.append((simplex, filt))

        # Landmarks array
        points = np.vstack([x, y]).T
        landmarks = points[landmarks_idx]
        
        # Sort by filtration for color mapping
        edges.sort(key=lambda x: x[1])

        # Plot the edges
        for k, ((a, b), filt) in enumerate(edges):
            t = k/(len(edges)-1+1e-12)
            color = cm.Greens(t)
            x = [landmarks[a,0], landmarks[b,0]]
            y = [landmarks[a,1], landmarks[b,1]]
            plt.plot(x, y, color=color, linewidth=1)
           
        # Better edge visualization 
        plt.xlim(-0.5, 0.5)
        plt.ylim(-2.5, -1.5)
        plt.title("Complex Edges")


    plt.legend()
    plt.savefig(f"{name}")
    plt.close()
    
import matplotlib.cm as cm
    
def desilva_circle():
    x, y, z = circle(100, 2)

    simplex_tree, dim, verts, simplices, betti = desilva_carlsson_witness_complex(x, y, z, name="circle", n_landmarks=50, max_dim=2)

    N = len(x)
    landmarks_idx = np.linspace(0, N-1, 50, dtype=int)
    plot_circle_points(x, y, landmarks_idx, "Circle_Points.png")
    plot_circle_points(x, y, landmarks_idx, "Circle_Edges.png", simplex_tree, True)

    
## TORUS #####################################################################################################################################
##############################################################################################################################################
def classical_vs_wiggly():
    # Keep these the same
    R = 15
    r = 2
    w1 = 1
    p = (1+np.sqrt(5))/2
    q = 1
    
    # testing different n's 
    n = 300
    
    # Wiggle parameters
    # Try alpha >= 1 and note betti numbers
    a = 0.3
    k = 5
    
    x_c, y_c, z_c, t_c = make_torus(R, r, w1, p, q, 10e-3, n, "classical", True)
    x_w, y_w, z_w, t_w = make_torus(R, r, w1, p, q, 10e-3, n, "wiggly", True, a, k)
    
    alpha = 100.0
    
    print("Constructing Classical Witness Complex")
    (_, dim_c, vert_c, simplices_c, betti_c, _) = desilva_carlsson_witness_complex(x_c, y_c, z_c, "Classical_Torus.png")
    
    print("Constructing Wiggly Witness Complex")
    (_, dim_w, vert_w, simplices_w, betti_w, _) = desilva_carlsson_witness_complex(x_w, y_w, z_w, "Wiggly_Torus.png")
    
    table = PrettyTable()
    table.field_names = ["Type", "Dimension", "Vertices", "Simplices", "Betti Numbers"]
    table.add_row(["Classical", dim_c, vert_c, simplices_c, str(betti_c)])
    table.add_row(["Wiggly", dim_w, vert_w, simplices_w, str(betti_w)])
    print(table)
    
def test_ns():
    R = 15
    r = 2
    w1 = 1
    p = (1 + np.sqrt(5)) / 2
    q = 1
    ns = [10, 100, 300, 500, 900]

    table = PrettyTable()
    table.field_names = ["n", "dim", "vertices", "simplices", "betti"]

    for n in ns:
        print(f"Constructing Complex for n = {n}")
        x, y, z, t = make_torus(R, r, w1, p, q, 10e-3, n, "classical", True)
        _, dim, vert, simplices, betti, _ = desilva_carlsson_witness_complex(x, y, z, name=f"Classical_n{n}.png")
        table.add_row([n, dim, vert, simplices, betti])

    print(table)   
    
       
def test_k():
    R = 15
    r = 2
    w1 = 1
    p = (1+np.sqrt(5))/2
    q = 1
    
    x_c, y_c, z_c, t_c = make_torus(R, r, w1, p, q, 10e-3, 300, "classical", True)
    
    ks = [0, 1, 2, 3, 4, 5]
        
    table = PrettyTable()
    table.field_names = ["k", "Dimension", "Vertices", "Simplices", "Betti Numbers"]
    
    for k in ks:
        print("Constructing Classical Witness Complex")
        (simplex_tree, dim_c, vert_c, simplices_c, betti_c, landmarks) = desilva_carlsson_witness_complex(x_c, y_c, z_c, "Classical_Torus.png", n_landmarks = 50, max_dim = 3, k=k)
        table.add_row([k, dim_c, vert_c, simplices_c, str(betti_c)])
        
        fig, ax = plt.subplots(figsize=(5,5))
        # ax = fig.add_subplot(projection = '3d')
        # Extract edges
        edges = []
        for simplex, filt in simplex_tree.get_filtration():
            if len(simplex) == 2:
                edges.append((simplex, filt))

        edges.sort(key=lambda x: x[1])
        num_edges = len(edges)

        # Plot edges
        for k, ((a, b), filt) in enumerate(edges):
            t = k / (num_edges - 1 + 1e-12)
            color = cm.viridis(t)
            xa, ya, za = landmarks[a]   
            xb, yb, zb = landmarks[b]
            ax.plot([xa, xb], [ya, yb], color=color, linewidth=0.8)

        # Plot landmarks
        ax.scatter(landmarks[:,0], landmarks[:,1], color='gold', s=10)

        ax.set_title(f"K = {k}")
        ax.set_aspect('equal')
        ax.axis('off')

        fig.savefig(f"Torus_k_{k}.png")
        plt.close(fig)
    print(table)
    
## MAIN ######################################################################################################################################
##############################################################################################################################################
if __name__ == "__main__":
    # desilva_circle()           # simple demo desilva circle
    # test_ns()                  # varying number of points in complex, desilva
    # classical_vs_wiggly()      # perturbing torus, desilva
    test_k()                     # k nearest neighbors, desilva  