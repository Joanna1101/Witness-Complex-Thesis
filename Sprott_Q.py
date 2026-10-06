## Imports ###################################################################################################################################
##############################################################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from scipy.integrate import solve_ivp
    
    
    
## System Q My Beloved #######################################################################################################################
##############################################################################################################################################
def Q(w, d, f):
    """
    System Q: 
    dx/dt = -z
    dy/dt = x-y
    dz/dt = dx + y^2 + fz
    
    Inputs: w.......(list) samples ??? 
            d, f....(floats) parameters
    Outputs: ????
    """
    # WHERE DOES T EVEN GO IN HERE??
    x, y, z = w
    wdot = np.zeros(3)
    
    wdot[0] = -z
    wdot[1] = x-y
    wdot[2] = d*x + y**2 + f*z
    
    return wdot
    
def makeQ(x0, d, f, tmax, tstep):
    # Integrate from 1 to tmax/10
    # Scipy wants ode as a function of t and w, but Q is autonomous
    soln = solve_ivp(lambda t, w: Q(w, d, f), (1.0, tmax/10), x0, method = "RK45")
    
    # shape (3, number of timesteps) [0,:] is x, [1,:] is y, [2, :] is z
    x0 = soln.y[:, -1] 
    
    # Build trajectory
    W = []
    T = []
    time = 1
    while time <= tmax:
        soln = solve_ivp(lambda t, w: Q(w, d, f), (time, time + tstep), x0, method = "RK45")
        W.append(soln.y.T)
        T.append(soln.t)
        x0 = soln.y[:, -1]
        time += tstep
    
    W = np.vstack(W)
    T = np.concatenate(T)
    x = W[:, 0]
    y = W[:, 1]
    z = W[:, 2]
    
    # fig = go.Figure(data = [go.Scatter3d(x=x, y=y, z=z, mode="lines", line = dict(width = 3, color = T, colorscale = "viridis"))])
    fig = go.Figure(data = [go.Scatter3d(x=x, y=y, z=z, mode="lines", line = dict(width = 3, color = "deepskyblue"))])
    fig.update_layout(title = f"Q: f = {f} and d = {d}", scene = dict(xaxis_title = "x", yaxis_title = "y", zaxis_title = "z"))
    fig.write_html(f"Q_f_{f}_d_{d}.html")

    return x, y, z
    
    
    
## MAIN ######################################################################################################################################
##############################################################################################################################################
if __name__ == "__main__":
    f = 0.5
    d = 3.1
    tmax = 1000
    tstep = 0.3   
    x0 = np.array([1.0, 1.0, 1.0])
    
    x, y, z = makeQ(x0, d, f, tmax, tstep)
