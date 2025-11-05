#!/usr/bin/env python

import sys
import numpy as np
from pylab import *

def calculate_apd (dt, vm, perc=0.9, verbose=False):
    """
    APD90, ou seja, duracao do potencial de acao a 90% da fase de repolarizacao,
    calculado usando a diferenca entre o tempo de ativacao e o de repolarizacao.
    O tempo de ativacao e' o tempo de derivada maxima de um PA.
    O tempo de repolarizacao e' o instante em que o potencial da membrana atinge
    o nivel correspondente a 10% da amplitude do PA apos a ativacao.

    >>> filename='testdata/test_apd.dat'
    >>> calculate_apd (filename)
    Activation time:           001.0000 ms
    Repolarization time:       163.4500 ms
    Action potential duration: 162.4500 ms
    162.44999999999999
    """
    tm = np.arange(0,len(vm)*dt,dt)
    
    ### find the activation time
    idmax = 0
    dvmax = 0.0
    dt = tm[1] - tm[0]
    size = len(vm)-1
    
    for i in range(size):
        dv = vm[i+1] - vm[i]
        if (dv > dvmax):
            dvmax = dv
            idmax = i
            
    act_time = tm[idmax]
    
    amp = np.max(vm) - np.min(vm)           # find the AP amplitude
    
    #### find the repolarization time
    vchk = amp*(1.0-perc) + np.min(vm)
    imax = vm.argmax()
    imin = vm.argmin()
    for i in range(imax,size):
        if ((vm[i] <= vchk)):
            rep_time = tm[i]
            break
    
    apd = np.abs(rep_time - act_time)       # calculate APD
    
    if (verbose):
        print('Activation time:           %08.4f ms' % act_time)
        print('Repolarization time:       %08.4f ms' % rep_time)
        print('Action potential duration: %08.4f ms' % apd)
   
    return apd

def calculate_apd_file (data_file, perc=0.9):
    """
    APD90, ou seja, duracao do potencial de acao a 90% da fase de repolarizacao,
    calculado usando a diferenca entre o tempo de ativacao e o de repolarizacao.
    O tempo de ativacao e' o tempo de derivada maxima de um PA.
    O tempo de repolarizacao e' o instante em que o potencial da membrana atinge
    o nivel correspondente a 10% da amplitude do PA apos a ativacao.

    >>> filename='testdata/test_apd.dat'
    >>> calculate_apd (filename)
    Activation time:           001.0000 ms
    Repolarization time:       163.4500 ms
    Action potential duration: 162.4500 ms
    162.44999999999999
    """
    data = np.loadtxt(data_file)            # format: %f %f \n (time voltage)
    tm = data[:,0]                          # time
    vm = data[:,1]                          # transmembrane potential

    #print(tm, vm)
    
    ### find the activation time
    idmax = 0
    dvmax = 0.0
    dt = tm[1] - tm[0]
    size = len(vm)-1
    
    for i in range(size):
        dv = vm[i+1] - vm[i]
        if (dv > dvmax):
            dvmax = dv/dt
            idmax = i
            
    act_time = tm[idmax]
    
    amp = np.max(vm) - np.min(vm)           # find the AP amplitude
    
    #### find the repolarization time
    vchk = amp*(1.0-perc) + np.min(vm)
    stop = False
    for i in range(size):
        dv = vm[i+1] - vm[i]
        if ((dv<0.) and (vm[i] < vchk) and (not stop)):
            rep_time = tm[i]
            stop = True
    
    apd = np.abs(rep_time - act_time)       # calculate APD
    
    print('Activation time:           %08.4f ms' % act_time)
    print('Repolarization time:       %08.4f ms' % rep_time)
    print('Action potential duration: %08.4f ms' % apd)
   
    return apd
    
if __name__ == "__main__":
    
    #import doctest
    #doctest.testmod()

    if len(sys.argv) < 3:
        print("\n Usage: calcAPD.py <file.dat> <timestep>\n")
        exit(1)

    #filename = sys.argv[1]
    #apd = calculate_apd (filename)

    filename = sys.argv[1]
    data = np.loadtxt(filename)            # format: %f %f \n (time voltage)
    tm = data[:,0]                          # time
    vm = data[:,1]   

    dt = float(sys.argv[2])
    step = tm[1]-tm[0]
    print('timestep %f' % step)

    apd = calculate_apd(dt, vm, verbose=True)
    print("APD90 = %f" % apd)
    

    

