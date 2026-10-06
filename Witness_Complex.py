## IMPORTS ###################################################################################################################################
##############################################################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import gudhi as gd


## HELPERS ###################################################################################################################################
##############################################################################################################################################
def plotComplex(witnesses, landmarks, st, name, max_dim):
    """
    Plots the simplicial complex    
    """
    fig = plt.figure(figsize=(8,8))
    ax = fig.add_subplot(111, projection='3d')

    # ax.scatter(witnesses[:,0], witnesses[:,1], witnesses[:,2], c = "blue", s=5, alpha=0.3, label = "Witnesses")
    ax.scatter(landmarks[:,0], landmarks[:,1], landmarks[:,2], c = "red", s=10, alpha= 0.5, label = "Landmarks")

    # Edges BETWEEN LANDMARKS
    for simplex in st.get_skeleton(max_dim):
        if len(simplex[0]) == 2:  # wtf is going on here
            i, j = simplex[0]
            p, q = landmarks[i], landmarks[j] 
            ax.plot([p[0], q[0]], [p[1], q[1]], [p[2], q[2]], 'k-', linewidth=0.5, alpha=0.3)

    ax.set_box_aspect([1,1,1])
    ax.legend(loc="upper right")
    plt.savefig(f"{name}.png")
    plt.close(fig)   

## SIMPLEX TREE ##############################################################################################################################
##############################################################################################################################################
def make_distance_matrix(landmarks, witnesses, eps):
    """
    Constructs distance matrix for a given simplex tree
    Inputs:  landmarks......(np array)
             witnesses......(np array)
             eps............(float) Radius for which a witness "witnesses"  a landmark
    Outputs: D..............(np array) distance if distance < R, otherwise -1.
    """
    distances = np.linalg.norm(landmarks[:, None, :]-witnesses[None, :, :], axis = 2)
    D = np.where(distances <= eps, distances, -1)
    
    # lets not use a double for loop lol
    # D = np.full((len(landmarks), len(witnesses)), -1.0)
    # for i in range(len(landmarks)):
    #     for j in range(len(witnesses)):
    #         d = np.linalg.norm(landmarks[i] - witnesses[j])

    #         if d <= eps:
    #             D[i, j] = d
    
    # For finding a good epsilon
    # distances = np.linalg.norm(landmarks[:, None, :] - witnesses[None, :, :], axis=2)
    # print("Min landmark-witness distance:", distances.min())
    # print("Avg:", distances.mean())
    
    # LL = np.linalg.norm(landmarks[:, None, :] - landmarks[None, :, :], axis=2)
    # print("Nearest landmark distances:")
    # for i in range(len(landmarks)):
    #     distances = LL[i].copy()
    #     distances[i] = np.inf
    #     print(i, distances.min())

    return D


def fuzzy_simplex_tree(landmarks, witnesses, max_dim = 3, eps = 0, R = 0):
    """
    Connects landmarks if they share a witness
    Inputs:  landmarks......(np array)
             witnesses......(np array)
             max_dim........(int) max simplex dimension
             eps............(float) radius for landmark-witness
             R..............(int) radius for landmark-landmark
    Outputs: st.............(simplex tree) [((landmark index), filtration)]
    """
    # TODO: more intelligent filtration
    dim = 1
    st = [] 
    
    st = [((i,), 0) for i in range(len(landmarks))] 
    
    # keep track of newest added batch of simplices
    prev_sims = []
    new_sims = []
    
    D = make_distance_matrix(landmarks, witnesses, eps)    
    filt = 1

    # Is it even working
    # print("D shape:", D.shape)
    # print("Total landmark-witness pairs:", D.size)
    # print("Valid pairs:", np.sum(D != -1))

    # Pairwise compare landmarks
    if max_dim >= 2:
        witnessed = D >= 0
        shared = witnessed@witnessed.T
        
        for i in range(len(landmarks)):
            for j in range(i+1, len(landmarks)):
                if shared[i, j] > 0:
                    new_sims.append((i, j))
                    st.append(((i, j), filt))
        dim += 1
        filt += 1
        
    while dim <= max_dim:
        Dhat = np.zeros((len(new_sims), len(witnesses)))
        
        for i in range(len(new_sims)):     # simplices
            for j in range(D.shape[1]):  # witnesses
                # find closest landmark
                closest_d = np.inf
                closest_p = None
                
                for l in new_sims[i]:
                    if D[l, j] < closest_d and D[l, j] != -1:
                        closest_d = D[l, j] 
                        closest_p = landmarks[l]
                
                # No landmark is witnessed by current witness
                if closest_p is None: continue
                
                # are all landmarks <= R of closest
                witnessed = True
                for l in new_sims[i]:
                    L_dist = np.linalg.norm(landmarks[l] - closest_p)
                    if L_dist > R:
                        witnessed = False
                        break
                
                if witnessed:
                    Dhat[i, j] = 1 # TODO: make this a filtration value

        # Update
        prev_sims = new_sims
        new_sims = []
        
        # Pairwise compare
        for i in range(Dhat.shape[0]):     # simplices
            for j in range(Dhat.shape[0]): # simplices
                if i != j:
                    ithRow = Dhat[i]
                    jthRow = Dhat[j]
                    
                    for w in range(Dhat.shape[1]):      # witnesses
                        if ithRow[w] == jthRow[w] == 1:
                            union = tuple(sorted(set(prev_sims[i]).union(prev_sims[j])))
                            if len(union) - 1 == dim:
                                st.append((union, filt))
                                new_sims.append(union)
                            break 
                            
        # Finally, increase dim and filtration                    
        dim += 1
        filt += 1
                            
    return st         



                
## FUZZY WITNESS #############################################################################################################################
##############################################################################################################################################  
def fuzzy_witness_complex(x, y, z, name, n_landmarks, max_dim, eps, R):
    """
    Inputs:  x, y, z......(List, List, List) 3d mapping
             name.........(String) title for figures
             n_landmarks..(int) number of landmarks (equispaced in time)
             max_dim......(int) max simplex dimension
             eps..........(float) landmark-witness radius
             R............(float) landmark-landmark radius
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
        
    # I think eps < R??? 
    st = fuzzy_simplex_tree(landmarks, witnesses, max_dim, eps, R)

    
    simplex_tree = gd.SimplexTree()
    for simplex, filt in st:
        simplex_tree.insert(simplex, filtration = filt)
    
    # Persistent Homology
    print("Plotting")
    bar_codes = simplex_tree.persistence()
    
    fig = plt.figure(figsize = (6,6))
    gd.plot_persistence_diagram(bar_codes)
    plt.savefig(f"Witness_Persistence_{name}")
    plt.close(fig) 
    
    plotComplex(witnesses, landmarks, simplex_tree, name, 2)
    
    return (simplex_tree,
            simplex_tree.dimension(),
            simplex_tree.num_vertices(),
            simplex_tree.num_simplices(),
            simplex_tree.betti_numbers(),
            landmarks)
    