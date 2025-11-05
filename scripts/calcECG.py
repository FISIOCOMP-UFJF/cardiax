#!/usr/bin/env python
# coding: utf-8
#
# This script is used to calculate the transmural ECG
# It requires two files with Ve from two nodes extracted
# using Paraview in .csv format.
#
# Bernardo M. Rocha
#
import os, sys
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

  if (len(sys.argv) < 3):
    print("\n Usage: calcECG [endo] [epi] <endo2> <epi2>\n")
    sys.exit(1)

  endfile1 = sys.argv[1]
  epifile1 = sys.argv[2]

  # normal
  end1 = np.loadtxt(endfile1, delimiter=',', skiprows=1)
  epi1 = np.loadtxt(epifile1, delimiter=',', skiprows=1)

  end1 = end1[:,0]
  epi1 = epi1[:,0]

  ve1 = epi1 - end1

  if (len(sys.argv) == 5):
    endfile2 = sys.argv[3]
    epifile2 = sys.argv[4]
    end2 = np.loadtxt(endfile2, delimiter=',', skiprows=1)
    epi2 = np.loadtxt(epifile2, delimiter=',', skiprows=1)
    end2 = end2[:,0]
    epi2 = epi2[:,0]
    ve2 = epi2 - end2

  plt.plot(ve1, lw=1.5, label='ve1')  
  if (len(sys.argv) == 5):
    plt.plot(ve2, lw=1.5, label='ve2')

  plt.xlabel('Tempo')
  plt.ylabel('ECG (mV)')
  plt.legend(loc='best')
  plt.show()

  # save data 1
  savepath1 = os.path.dirname(endfile1)
  savefile1 = savepath1 + '/ecg.txt'  
  np.savetxt(savefile1,ve1)

  # save data 2
  if (len(sys.argv) == 5):
    savepath2 = os.path.dirname(endfile2)
    savefile2 = savepath2 + '/ecg.txt'
    np.savetxt(savefile2,ve2)