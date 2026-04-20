## Imports ###################################################################################################################################
##############################################################################################################################################
from Witness_Complex_Functions import *
from make_torus import *

## TESTS #####################################################################################################################################
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
    (dim_c, vert_c, simplices_c, betti_c) = desilva_carlsson_witness_complex(x_c, y_c, z_c, "Classical_Torus.png")
    
    print("Constructing Wiggly Witness Complex")
    (dim_w, vert_w, simplices_w, betti_w) = desilva_carlsson_witness_complex(x_w, y_w, z_w, "Wiggly_Torus.png")
    
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
        dim, vert, simplices, betti = alpha_witness_complex(x_c, y_c, z_c, alpha, f"{alpha}_Classical_Torus.png")
        table.add_row([alpha, dim, vert, simplices, betti])

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
        dim, vert, simplices, betti = desilva_carlsson_witness_complex(x, y, z, name=f"Classical_n{n}.png")
        table.add_row([n, dim, vert, simplices, betti])

    print(table)    

def test_desilva_circle():
    print("blep")
    
## MAIN ######################################################################################################################################
##############################################################################################################################################
if __name__ == "__main__":
    test_ns()
    # classical_vs_wiggly()
    # test_alphas()