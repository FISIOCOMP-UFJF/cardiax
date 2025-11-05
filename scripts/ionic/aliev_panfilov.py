# Aliev & Panfilov cell model based on the paper by Goektepe & Kuhl (2010)

import sys
from pylab import *
from scipy.integrate import odeint
from scipy.integrate import ode
from calcAPD import calculate_apd

def aliev_panfilov(time, states, params):

    # parameters
    a, b, k = 0.01, 0.15, 8.0
    eps, mu1, mu2 = 0.002, 0.2, 0.3
    kv  = 0.005

    # twitch function parameters
    e0, einf, psi, vb, vr = 0.1, 1, 1, 0, -80
    
    # stimulation parameters
    istart, isend  = 0.1, 10000000000.0
    idurat   = 1.0
    period   = params['bcl']
    istimamp = 0.5
    
    cond1 = ( time >= istart )
    cond2 = ( time <= isend )
    cond3 = ( ((time-istart) - (floor((time-istart)/period)*period))  <= idurat )

    if(cond1 and cond2 and cond3):
        istim = istimamp
    else:
        istim = 0.0

    # state variables
    v,r,ta = states[0], states[1], states[2]
    
    term1  = eps + (mu1*r)/(mu2+v)    
    twitch = lambda v: e0 + (einf-e0)*exp(-exp(-psi*(v-vb)))  
    
    values = zeros(3)
    
    values[0] = k*v*(v-a)*(1-v) - r*v + istim
    values[1] = term1 * (-r -k*v*(v-b-1))
    values[2] = twitch(v) * (kv*v - ta)
    
    return values    

def solve_with_bcl(bcl):
    
    dt      = 0.1
    size    = bcl/dt
    stims   = 6
    states0 = [0,0,0]

    T = bcl*stims
    n = T/dt
    t = linspace(0, T, n)
    
    params  = {}
    params['bcl'] = bcl

    # Using ode
    r = ode(aliev_panfilov).set_integrator('dopri5')
    r.set_initial_value(states0, t[0])
    r.set_f_params(params)

    sol = []

    while r.successful() and r.t < T:
        r.integrate(r.t + dt)
        sol.append( [ r.t, r.y[0], r.y[1], r.y[2] ] )

    sol = array(sol)
    
    # handle exception when len(sol[:,1]) is n+1
    if(shape(sol[:,1])[0] > stims*size):
        sol = array(sol[:-1])

    i0, i1 = (stims-2)*size, (stims-1)*size
    apd = calculate_apd(t[1]-t[0], sol[i0:i1,1])
    print(" APD90 = %f" % (apd))
    return apd

def convert_vm_physiological(t,v):
    time = 12.9 * t
    vm   = 100.0 * v - 80.0
    return (time, vm)
    
if __name__ == "__main__":

    if (len(sys.argv) > 1):
        bcl = float(sys.argv[1])
    else:
        bcl = 80.0

    dt      = 0.1
    size    = bcl/dt
    stims   = 6
    states0 = [0,0,0]
    T       = bcl*stims
    
    params  = {}
    params['bcl'] = bcl

    t = arange(0, (stims*bcl)+dt, dt)

    print(" Start       = %f" % 0.0)
    print(" Stop        = %f" % T)
    print(" BCL         = %f" % bcl)
    print(" Total size  = %d" % len(t))
    print(" Slice size  = %d" % int(size))
   
    # Using ode
    r = ode(aliev_panfilov).set_integrator('dopri5')
    r.set_initial_value(states0, t[0])
    r.set_f_params(params)
    sol = []
    
    while r.successful() and r.t < T:
        r.integrate(r.t + dt)
        sol.append( [ r.t, r.y[0], r.y[1], r.y[2] ] )

    sol = array(sol)
    # handle exception when len(sol[:,1]) is n+1
    if(shape(sol[:,1])[0] > stims*size):
        sol = array(sol[:-1])

    #time,vm = convert_vm_physiological(sol[:,0],sol[:,1])
    #plot(time,vm,label='V')
    #xlabel('time [ms]')
    #ylabel('Vm [mV]')
    #grid(True)
    #show()

    plot(sol[:,0],sol[:,1],label='V')
    plot(sol[:,0],sol[:,2]*0.4,label='r')
    plot(sol[:,0],sol[:,3]*100,label='Ta')
    grid(True)
    show()       
    
    # reference for restitution curve
    apd0 = solve_with_bcl(500)
    print(" APD_ref = %f\n" % (apd0))        

    # compute and plot restitution curve
    apdlst = []
    bcllst = [25,30,40,50,80,90,100,120,150,200,220,250,300,350,400,450,500]
    for b in bcllst:
        apd = solve_with_bcl(b)
        apdlst.append(apd)

    apds = array(apdlst)
    bcls = array(bcllst)
    plot(bcls/apd0, apds/apd0, 'k-o')
    xlabel('$T/APD_0$')
    ylabel('$APD/APD_0$')
    title('Restitution curve')
    show()
