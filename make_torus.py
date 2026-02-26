## IMPORTS ############################################################
#######################################################################
import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction




## TORUS MAPPINGS #####################################################
#######################################################################
def classical_torus(R, r, w1, w2, t):
    """
    Returns x,y,z embedding of classical round 2-torus in R3. 
    
    Inputs: R.........(Float) major radius
            r.........(Float) minor radius
            w1........(Float) angular velocity 1
            w2........(Float) angular velocity 2
            t.........(List) sequence of times
    Outputs: x,y,z....(List, List, List) mapping to R^3
    """
    theta1 = w1*t # Major rotation
    theta2 = w2*t # Minor rotation

    x = (R+r*np.cos(theta2))*np.cos(theta1)
    y = (R+r*np.cos(theta2))*np.sin(theta1)
    z = r*np.sin(theta2)
    
    return x,y,z

def mobius_torus(R, r, w1, w2, t):
    """
    Returns x,y,z embedding of twisted mobius-like 2-torus in R3. 
    
    Inputs: R.........(Float) major radius
            r.........(Float) minor radius
            w1........(Float) angular velocity 1
            w2........(Float) angular velocity 2
            t.........(List) sequence of times
    Outputs: x,y,z....(List, List, List) mapping to R^3
    """
    theta1 = w1*t # Major rotation
    theta2 = w2*t # Minor rotation

    #TODO: Implement here 
    x = (R+r*np.cos(theta2))*np.cos(theta1)
    y = (R+r*np.cos(theta2))*np.sin(theta1)
    z = r*np.sin(theta2)
    
    return x,y,z
    
    
    
    
## HELPERS ############################################################ 
#######################################################################
def plot_orbit(x, y, z, t, w1, w2, name):
    """
    Plots colored torus mapping
    Inputs:  x, y, z....(List, List, List) mapping
             t..........(List) time sequence for colors
             w1, w2.....(Float, Float) angular velocities for plot title
             name.......(String) desired plot name
    Outputs: png with desired name
    """
    fig = plt.figure(figsize = (8,6))
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(r'$\omega_1 = $ ' +f'{w1:.4f}'+', $\omega_2 = $ ' +f'{w2:.4f}')
    colors = (t - t.min()) / (t.max() - t.min())
    ax.scatter(x, y, z, c=colors, cmap='cool', s=1)
    ax.set_box_aspect([1,1,1])
    plt.savefig(f"{name}")
    plt.close()

def continued_frac(alpha, eps, max_its = 1000, step = 5):
    """
    Performs continued fraction approximation of angular velocity ratio
    Inputs:  alpha....(Float) w2/w1 floating point velocity ratio
             eps......(Float) precision, eg. 10e-3
             max_its....(Int) maximum iterations increasing denominator limit
             step.......(Int) increase to denominator limit at each step
    Outputs: p, q.....(Int, Int) p/q rational approximation
    """
    p_best = 0
    q_best = 1
    err = float("inf")
    limit = 1
    
    for i in range(max_its):
        pq_ratio = Fraction(alpha).limit_denominator(limit)
        p = pq_ratio.numerator
        q = pq_ratio.denominator
        err = np.abs(alpha - pq_ratio)
        
        if err <= eps:
            print(f"Reached precision of {err:.4f} after {i} iterations.")
            print(f"Rational approximation p/q: p = {p}, q = {q}")
            return p,q 
        
        limit += step
        
    print(f"Reached precision of {err:.4f} after {max_its} iterations.")
    print(f"Rational approximation p/q: p = {p}, q = {q}")
    return p, q




## MAKE TORUS #########################################################
#######################################################################
def make_torus(R, r, w1, p, q, eps, type, plot):
    """
    Computes second angular velocity up to desired precision
    Returns coordinates in R^3 and time series
    
    Inputs: R.........(Float) major radius
            r.........(Float) minor radius
            w1........(Float) angular velocity 1
            p.........(Float or Int) numerator of velocity ratio
            q.........(Float or Int) denominator of velocity ratio
            eps.......(Float or String) float -> rational approximation precision
                                        "exact" -> double precision float
            type......(String) mapping type: classical, mobius, etc
            plot......(Boolean) if true, plots system
    Outputs: x,y,z....(List, List, List) mapping to R^3
    """
    if eps == "exact":
        w2 = w1*(p/q)
        t = np.linspace(0, 2*np.pi*10, 10000) # 10 loops, 10,000 samples
    else:    
        w2 = w1 * (p/q)
        
        if 2*np.pi*p/w1 == 2*np.pi*q/w2:
            print(r"$\omega_1/\omega_2$ is rational, period is exact.")
            T = 2*np.pi*p/w1
        else: 
            print(r"$\omega_1/\omega_2$ is rational, period is approximate.")
            alpha = w2/w1
            p,q = continued_frac(alpha, eps)
            T = 2*np.pi*p/w1
          
        t = np.linspace(0, T, 1000*p)
            
    if type == "classical":
        x,y,z = classical_torus(R, r, w1, w2, t)
    elif type == "mobius":
        x,y,z = mobius_torus(R, r, w1, w2, t)
    else:
        print("Unrecognized Mapping Type")
        return []
    
    if plot:
        name = f"{type}_eps={eps}_w1={w1:.4f}_w2={w2:.4f}.png"
        plot_orbit(x, y, z, t, w1, w2, name)
        
    return x, y, z, t



    
## EXAMPLE USAGE ######################################################
#######################################################################
def rational_classical_torus():
    R = 10
    r = 5
    w1 = 1
    p = 2
    q = 1
    x, y, z, t = make_torus(R, r, w1, p, q, "exact", "classical", True)
    
def pi_phi_classical_exact():
    R = 10
    r = 5
    w1 = 1
    q = 1
    x1, y1, z1, t1 = make_torus(R, r, w1, np.pi, q, "exact", "classical", True)
    x2, y2, z2, t2 = make_torus(R, r, w1, (1+np.sqrt(5))/2, q, "exact", "classical", True)
        
def pi_phi_classical_10e6():
    R = 10
    r = 5
    w1 = 1
    q = 1
    x1, y1, z1, t1 = make_torus(R, r, w1, np.pi, q, 10e-3, "classical", True)
    x2, y2, z2, t2 = make_torus(R, r, w1, (1+np.sqrt(5))/2, q, 10e-3, "classical", True)
    
    
    
    
## TEMPORARY MAIN #####################################################
#######################################################################
if __name__ == "__main__":
    # rational_classical_torus()
    # pi_phi_classical_exact()
    pi_phi_classical_10e6()
    
    
    
    