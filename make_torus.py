## IMPORTS ############################################################
#######################################################################
import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction




## TORUS MAPPINGS #####################################################
#######################################################################
def classical_torus(R, r, w1, w2, t, plot = False):
    """
    Returns x,y,z embedding of classical round 2-torus in R3. 
    
    Inputs: R.........(Float) major radius
            r.........(Float) minor radius
            w1........(Float) angular velocity 1
            w2........(Float) angular velocity 2
            t.........(List) sequence of times
            plot......(Boolean) plots surface shape if true
    Outputs: x,y,z....(List, List, List) mapping to R^3
    """
    theta1 = w1*t # Major rotation
    theta2 = w2*t # Minor rotation

    x = (R+r*np.cos(theta2))*np.cos(theta1)
    y = (R+r*np.cos(theta2))*np.sin(theta1)
    z = r*np.sin(theta2)
    
    if plot:
        u = np.linspace(0 ,2*np.pi, 200)
        v = np.linspace(0, 2*np.pi, 100)
        U, V = np.meshgrid(u,v)
        X = (R + r*np.cos(V)) * np.cos(U)
        Y = (R + r*np.cos(V)) * np.sin(U)
        Z = r * np.sin(V)
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, rstride=2, cstride=2, color='mediumpurple', alpha=0.35, linewidth=0)
        plt.savefig("Classical_Torus_Surface.png")
        plt.close()
        
    return x,y,z
    
def wiggle_torus(R, r, w1, w2, t, a, k, plot = False):
    """
    Returns x,y,z embedding of twisted mobius-like 2-torus in R3. 
    
    Inputs: R.........(Float) major radius
            r.........(Float) minor radius
            w1........(Float) angular velocity 1
            w2........(Float) angular velocity 2
            t.........(List) sequence of times
            a.........(Int) amplitude of wiggle
            k.........(Int) number of wiggles
            plot......(Boolean) plots surface shape if true
    Outputs: x,y,z....(List, List, List) mapping to R^3
    """
    theta1 = w1*t # Major rotation
    theta2 = w2*t # Minor rotation
    r_wiggle = r*(1+a*np.sin(k*theta1)) # Radius modulation

    x = (R+r_wiggle*np.cos(theta2))*np.cos(theta1)
    y = (R+r_wiggle*np.cos(theta2))*np.sin(theta1)
    z = r_wiggle*np.sin(theta2)
    
    if plot:
        u = np.linspace(0 ,2*np.pi, 200)
        v = np.linspace(0, 2*np.pi, 100)
        U, V = np.meshgrid(u,v)
        r_wiggle_s = r*(1+a*np.sin(k*U))
        X = (R + r_wiggle_s*np.cos(V)) * np.cos(U)
        Y = (R + r_wiggle_s*np.cos(V)) * np.sin(U)
        Z = r_wiggle_s * np.sin(V)
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111, projection='3d')
        ax.plot_surface(X, Y, Z, rstride=2, cstride=2, color='mediumpurple', alpha=0.35, linewidth=0)         
        plt.savefig("Wiggly_Torus_Surface.png")
        plt.close()
    
    return x,y,z
    
    
    
## HELPERS ############################################################ 
#######################################################################
def plot_orbit(x, y, z, t, w1, w2, name, col):
    """
    Plots colored torus mapping
    Inputs:  x, y, z....(List, List, List) mapping
             t..........(List) time sequence for colors
             w1, w2.....(Float, Float) angular velocities for plot title
             name.......(String) desired plot name
             col........(String) color map or black
    Outputs: png with desired name
    """
    fig = plt.figure(figsize = (8,6))
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(r'$\omega_1 = $ ' +f'{w1:.4f}'+', $\omega_2 = $ ' +f'{w2:.4f}')
    if col == "black":
        ax.scatter(x, y, z, c = "black", s=0.5)
    else: 
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
def make_torus(R, r, w1, p, q, eps, samples, type, plot, a=0, k=0):
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
            samples...(Int) number of samples to use
            type......(String) mapping type: classical, wiggly
            plot......(Boolean) if true, plots system
            a, k......(Ints) for wiggly torus
    Outputs: x,y,z....(List, List, List) mapping to R^3
    """
    if eps == "exact":
        w2 = w1*(p/q)
        t = np.linspace(0, 2*np.pi*10, samples) # 10 loops, 10,000 samples
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
          
        t = np.linspace(0, T, samples)
            
    if type == "classical":
        x,y,z = classical_torus(R, r, w1, w2, t, True)
    elif type == "wiggly":
        x,y,z = wiggle_torus(R, r, w1, w2, t, a, k, True)
    else:
        print("Unrecognized Mapping Type")
        return []
    
    if plot:
        name = f"{type}_eps={eps}_w1={w1:.4f}_w2={w2:.4f}.png"
        plot_orbit(x, y, z, t, w1, w2, name, "color")
        
    return x, y, z, t



    
## EXAMPLE USAGE ######################################################
#######################################################################
def rational_classical_torus():
    R = 10
    r = 5
    w1 = 1
    p = 2
    q = 1
    n = 1000
    x, y, z, t = make_torus(R, r, w1, p, q, n, "exact", "classical", True)
    
def pi_phi_classical_exact():
    R = 10
    r = 5
    w1 = 1
    q = 1
    n = 1000
    x1, y1, z1, t1 = make_torus(R, r, w1, np.pi, q, n, "exact", "classical", True)
    x2, y2, z2, t2 = make_torus(R, r, w1, (1+np.sqrt(5))/2, q, n, "exact", "classical", True)
        
def pi_phi_classical_10e6():
    R = 10
    r = 5
    w1 = 1
    q = 1
    n = 1000
    x1, y1, z1, t1 = make_torus(R, r, w1, np.pi, q, 10e-3, "classical", True)
    x2, y2, z2, t2 = make_torus(R, r, w1, (1+np.sqrt(5))/2, q, 10e-3, "classical", True)
    
def pi_wiggly_10e6():
    R = 15
    r = 2
    w1 = 1
    q = 1
    n = 1000
    a = 0.3 # amplitude (> 1 causes weird flips)
    k = 5   # Number of wiggles
    x1, y1, z1, t1 = make_torus(R, r, w1, np.pi, q, 10e-3, n, "wiggly", True, a, k)
    
    
    
    
## TEMPORARY MAIN #####################################################
#######################################################################
if __name__ == "__main__":
    # rational_classical_torus()
    # pi_phi_classical_exact()
    # pi_phi_classical_10e6(
    pi_wiggly_10e6()    
    
    
    