#!/usr/bin/env python
# coding: utf-8
#
# Script to plot {ACT,APD,REP} along a 1D transmural
# specified by the (y,z) coordinates of the cube/slab
#
# Bernardo M. Rocha, Fev, 2014
#
import sys
from pylab import *

import matplotlib as mpl

w,h = 3, 3 / ((1 + math.sqrt(5)) / 2)

mpl.rcParams['axes.labelsize' ] = 11
mpl.rcParams['xtick.labelsize'] = 10
mpl.rcParams['ytick.labelsize'] = 10
mpl.rcParams['legend.fontsize'] = 10
mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['font.serif'] = 'Times'
mpl.rcParams['text.usetex'] = True
mpl.rcParams['figure.figsize'] = w, h

if(len(sys.argv)<5):
	print("\n Usage: plot_1D_transmural_line.py [mesh.pts] [apdmapcolor.dat] [y,z]\n")
	sys.exit(1)

ptsname = sys.argv[1]
apdname = sys.argv[2]
yc = float(sys.argv[3])
zc = float(sys.argv[4])

ylabstr = apdname.split('.')
ylabstr = ylabstr[0][-3:] 

f = open(ptsname)
apd = loadtxt(apdname)

npts = int(f.readline())
vpts = zeros((npts,3))
for i in range(npts):
	l = f.readline()
	p = l.split()
	vpts[i,0] = float(p[0])
	vpts[i,1] = float(p[1])
	vpts[i,2] = float(p[2])
f.close()

xmin = vpts[:,0].min()
xmax = vpts[:,0].max()*2

ry = array( where(vpts[:,1]==yc) )
rz = array( where(vpts[:,2]==zc) )
ry = ry[0,:]
rz = rz[0,:]

if (len(ry)==0 or len(rz)==0):
	print(" Erro: ponto da linha transmural nao encontrado.")
	sys.exit(1)

r=[]
for i in range(len(ry)):
	for j in range(len(rz)):
		if ry[i] == rz[j]: r.append(ry[i])

print(r)
                
x, apdline = zeros(len(r)), zeros(len(r))
for i in range(len(r)):
	idx = r[i]
	x[i] = vpts[idx,0]
	apdline[i] = apd[idx]
x = 100*(x-xmin)/xmax

print("Dispersao APD")
print("Min: %f" % apdline.min())
print("Max: %f" % apdline.max())
print("Dif: %f" % (apdline.max()-apdline.min())) 
print("\nX\t Dados")
for i in range(len(x)):
	print("%f %f" % (x[i],apdline[i]))

# plota grafico
plot(x,apdline,zorder=5,lw=1.25)
#xlabel('Endo $\hspace{0.25cm}$ Dist\^ancia transmural (\%) $\hspace{0.25cm}$ Epi')
xlabel('Endo $\hspace{0.4cm}$ Dist\^ancia transmural (\%) $\hspace{0.4cm}$ Epi')
#title('DPA (ms)', fontsize=11)
title('%s (ms)' % ylabstr, fontsize=11)

#grid(True, lw=0.5, color='0.75', zorder=0)
savefig('apd.pdf', bbox_inches='tight')
show()

