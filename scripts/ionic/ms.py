# Mitchell Schaeffer Model

from pylab import *
from scipy.integrate import ode

def mitchell_schaeffer(time, states, params):

    # parameters
    v_gate     = 0.13
    tau_in     = 0.3   # ms
    tau_out    = 6.0   # ms
    tau_open   = 120.0 # ms
    tau_closed = 150.0 # ms   
    
    # stimulation parameters
    idurat   = 1.0
    istimamp = 0.1
    istart   = 0.1
    isend    = istart + idurat
    if(time >= istart and time <= isend):
        Jstim = istimamp
    else:
        Jstim = 0.0

    # state variables
    v, h = states
    values = zeros(2)
    
    Jin   = (h*(v*v*(1.0-v)))/tau_in
    Jout  = -v/tau_out

    values[0] = Jin + Jout + Jstim
    values[1] = (1.0-h)/tau_open if v < v_gate else -h/tau_closed
  
    return values    
   
if __name__ == "__main__":

    # Simulation setup
    T       = 500.0
    dt      = 0.1
    states0 = [0,1]
    params  = {}

    t = arange(0, T+dt, dt)
  
    # Using ode
    r = ode(mitchell_schaeffer).set_integrator('dopri5')
    r.set_initial_value(states0, t[0])
    r.set_f_params(params)
    sol = []
    
    while r.successful() and r.t < T:
        r.integrate(r.t + dt)
        sol.append( [ r.t, r.y[0], r.y[1] ] )

    # Plot v
    sol  = array(sol)
    time = sol[:,0]
    v, h = sol[:,1], sol[:,2]
    plot(time,v,label='v')
    xlabel('time')
    ylabel('V')
    grid(True)
    show()

