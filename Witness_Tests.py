## IMPORTS ###################################################################################################################################
##############################################################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import gudhi as gd
from Witness_Complex import *
from make_torus import *
from Sprott_Q import *
from prettytable import PrettyTable


## TORUS #####################################################################################################################################
##############################################################################################################################################
def torusTest():
    # Keep these the same
    R = 15
    r = 2
    w1 = 1
    p = (1+np.sqrt(5))/2
    q = 1
    
    x, y, z, t = make_torus(R, r, w1, p, q, 10e-3, 500, "classical", True)
    print("made torus")
        
    # Boy Do I need to vectorize this >___<
    print("Constructing Fuzzy Witness Complex")
    (st, dim, vert, simplices, betti, lands) = fuzzy_witness_complex(x, y, z, "torus_phi_n150", n_landmarks = 100, max_dim = 3, eps = 3.5, R = 5)
    print("made complex")
    
    table = PrettyTable()
    table.field_names = ["Type", "Dimension", "Vertices", "Simplices", "Betti Numbers"]
    table.add_row(["Torus", dim, vert, simplices, str(betti)])
    print(table)

## Q!!!!!! ###################################################################################################################################
##############################################################################################################################################
def QTest():
    f = 0.5
    d = 3.1
    tmax = 1000
    tstep = 0.3   
    x0 = np.array([1.0, 1.0, 1.0])
    
    x, y, z = makeQ(x0, d, f, tmax, tstep)
    
    print("Constructing Fuzzy Witness Complex")
    (st, dim, vert, simplices, betti, lands) = fuzzy_witness_complex(x, y, z, "Q", n_landmarks = 100, max_dim = 3, eps = 3.5, R = 5)
    print("made complex")
    
    table = PrettyTable()
    table.field_names = ["Type", "Dimension", "Vertices", "Simplices", "Betti Numbers"]
    table.add_row(["Torus", dim, vert, simplices, str(betti)])
    print(table)
    
    

## MAIN ######################################################################################################################################
##############################################################################################################################################
if __name__ == "__main__":
    torusTest()