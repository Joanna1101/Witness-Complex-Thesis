import matplotlib.pyplot as plt
import numpy as np
from fractions import Fraction

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

    # Flow on torus mapped to R2 (wiggly torus?) 
    # cursed macaronis lol
    x = (R*(1+0*np.sin(theta1))+r*(1+0.8*np.sin(3*theta1))*np.cos(theta2))*np.cos(theta1)
    y = (R*(1+0*np.sin(theta1))+r*(1+0.8*np.sin(3*theta1))*np.cos(theta2))*np.sin(theta1)
    z = r*(1+0.8*np.sin(theta1))*np.sin(theta2)
    
    return x,y,z

def rational_torus():
    # Parameters
    R = 10
    r = 5
    w1 = 1
    
    # Rational angular velocities
    w2_1 = 3
    t_1 = np.linspace(0, (3/2)*np.pi*1, 1000)
    w2_2 = 2
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
    plt.savefig("Rational_Torus.png")
    
def rational_winding_torus():
    # Parameters
    R = 10
    r = 5
    w1 = 1
    
    # Rational angular velocities
    w2_1 = 3
    t_1 = np.linspace(0, 200*w2_1, 1000)
    w2_2 = 2
    t_2 = np.linspace(0, 200*w2_2, 1000)
    
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
    plt.savefig("Rational_Winding_Torus.png")
    
def irrational_torus():
    # Parameters
    R = 10
    r = 5
    w1 = 1
    
    # Irrational angular velocities
    w2_1 = np.pi
    t_1 = np.linspace(0, w2_1*np.pi*1, 1000)
    w2_2 = (1 + np.sqrt(5))/2
    t_2 = np.linspace(0, w2_2*np.pi*1, 1000)
    
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
    plt.savefig("Irrational_Torus.png")
    plt.close()

def irrational_winding_torus():
    # Parameters
    R = 10
    r = 5
    w1 = 1
    
    # Irrational angular velocities
    w2_1 = np.pi
    t_1 = np.linspace(0, 200*w2_1, 1000)
    w2_2 = (1 + np.sqrt(5))/2
    t_2 = np.linspace(0, 200*w2_2, 1000)
    
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
    plt.savefig("Irrational_Winding_Torus.png")
    
   
##  WORKING VERSION ################################################### 
#######################################################################
def plot_semiperiodic_orbit(R, r, w1, p, q):
    w2 = w1 * (p/q)
    
    # Rational case
    if 2*np.pi*p/w1 == 2*np.pi*q/w2:
        print("Rational")
        T = 2*np.pi*p/w1
    # Irrational case: continued fractions
    else: 
        print("Irrational")
        alpha = w2/w1
        
        # would be cool to do this within a precision: loop through different denominator limits
        # We want a bounded denominator.
        pq_ratio = Fraction(alpha).limit_denominator(10)
        p = pq_ratio.numerator
        q = pq_ratio.denominator
        print(f"p: {p}")
        print(f"q: {q}")
        T = 2*np.pi*p/w1
      
    # Use "period" to define time series  
    t = np.linspace(0, T, 1000*p)
        
    x,y,z = make_torus(R, r, w1, w2, t)
    
    # Visualization
    fig = plt.figure(figsize = (8,6))
    ax = fig.add_subplot(111, projection = '3d')
    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')
    ax.set_title(r'$\omega_1 = $ ' +f'{w1:.4f}'+', $\omega_2 = $ ' +f'{w2:.4f}')
    colors = (t - t.min()) / (t.max() - t.min())
    ax.scatter(x, y, z, c=colors, cmap='cool', s=1)
    ax.set_box_aspect([1,1,1])
    plt.savefig(f"w1 = {w1:.4f}, w2 = {w2:.4f}.png")
    plt.close()
    
    
if __name__ == "__main__":
    R = 10
    r = 5
    w1 = 1
    
    # If p, q integers the orbit will be periodic.
    # p = np.pi
    p = (1 + np.sqrt(5))/2
    q = 1
        
    plot_semiperiodic_orbit(R, r, w1, p, q)
    
    
    # Individual test functions 
    # rational_torus()
    # rational_winding_torus()
    # irrational_torus()
    # irrational_winding_torus()