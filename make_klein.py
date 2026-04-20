## IMPORTS ############################################################
#######################################################################
import numpy as np
import matplotlib.pyplot as plt
from make_torus import *

## MAPPING TO R4 ######################################################
#######################################################################
# Following Dłotko1, arXiv:2301.06753v1  See https://doi.org/10.1137/23M1594728
def klein_bottle(theta, phi):
    x = np.cos(0.5 * theta) * np.cos(phi) - np.sin(0.5 * theta) * np.sin(2 * phi)
    y = np.sin(0.5 * theta) * np.cos(phi) + np.cos(0.5 * theta) * np.sin(2 * phi)
    z = 8 * np.cos(theta) * (1 + 0.5 * np.sin(0.5 * phi))
    w = 8 * np.sin(theta) * (1 + 0.5 * np.sin(0.5 * phi))
    return x, y, z, w


def make_klein_bottle(w1, p, q, eps, samples, plot):
    """
    Makes dense orbit on klein bottle
    Inputs:  w1.........(float) angular velocity 1
             p,q........(float or int) ratio parameters for angular velocity 2
             eps........(float or "exact") rational approximation precision
             samples....(int) number of samples
             plot.......(bool) whether to plot two 3D projections
    Outputs: x,y,z,w....(np arrays) coordinates in R^4
             t..........(array) time samples
    """
    if eps == "exact":
        w2 = w1 * (p / q)
        t = np.linspace(0, 2 * np.pi * 10, samples)
    else:
        w2 = w1 * (p / q)

        # Rational ratio = exact period like torus
        if 2*np.pi*p/w1 == 2*np.pi*q/w2:
            print(r"$\omega_1/\omega_2$ is rational, period is exact.")
            T = 2 * np.pi * p / w1
        else:
            print(r"$\omega_1/\omega_2$ is rational, period is approximate.")
            alpha = w2 / w1
            p, q = continued_frac(alpha, eps)
            T = 2 * np.pi * p / w1

        t = np.linspace(0, T, samples)

    theta = w1 * t
    
    
    phis  = [0, 0.1, 0.5, 0.9, 1, 5]
    for phi in phis: 
        x, y, z, w = klein_bottle(theta, phi)

        # --- Plotting ---
        if plot: 
            fig1 = plt.figure()
            ax1 = fig1.add_subplot(111, projection='3d')
            ax1.plot(x, y, z, color='black', linewidth=0.8)
            ax1.set_title("(x, y, z) Projection")
            ax1.set_xlabel("x")
            ax1.set_ylabel("y")
            ax1.set_zlabel("z")
            plt.savefig(f"Klein_xyz_{phi}.png")

            fig2 = plt.figure()
            ax2 = fig2.add_subplot(111, projection='3d')
            ax2.plot(x, z, w, color='purple', linewidth=0.8)
            ax2.set_title("(x, z, w) Projection")
            ax2.set_xlabel("x")
            ax2.set_ylabel("z")
            ax2.set_zlabel("w")
            plt.savefig(f"Klein_xzw_{phi}.png")

            plt.show()

    return x, y, z, w, t

if __name__ == "__main__":
    # Same parameters as Klein.m:
    w1 = 1
    p=1
    q=1
    eps = "exact"
    samples = 100
    plot =True
    x, y, z, w, t = make_klein_bottle(w1, p, q, eps, samples, plot)