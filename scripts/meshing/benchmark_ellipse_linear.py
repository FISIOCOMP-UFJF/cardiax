
import numpy as np 
import math



def benchmark_ellipse(h, params):


	top = params[0]
	outer_long_axis = params[1]
	outer_short_axis = params[2]
	wall_thickness = params[3]
	#fiber angles
	fiber_angle_epi  = params[4]
	fiber_angle_endo = params[5]

	if len(h)==3:
		nab=h[0]
		ncirc=h[1]
		nr=h[2]
	else:
		outer=10
		ncirc = round(2*np.pi*(outer_long_axis-wall_thickness/2)/h)
		nab =  round(  ((outer_short_axis/2+top)/(2*outer_short_axis) * np.pi*math.sqrt((outer_short_axis-wall_thickness)**2+outer_long_axis**2) ) / h);
		nr = round(wall_thickness/h);



	#mesh generation: points to base nodes off
	points_circ  = ncirc
	points_ab    = nab + 1
	points_trans = nr  + 1

	iz = 1

	coords = np.zeros([iz,points_ab,5,points_circ+1])


	for long_r in np.arange(outer_long_axis-wall_thickness, outer_long_axis + wall_thickness/(points_trans-1), wall_thickness/(points_trans-1) ):

		short_r = outer_short_axis - (outer_long_axis-long_r)
		#print 'long axis: ', long_r
		#print 'short axis: ',short_r
		#print

		uu = -math.acos(top/long_r)
		u = np.array(np.arange(-np.pi, uu + (np.pi+uu)/(points_ab-1)/2, (np.pi+uu)/(points_ab-1) ) )
		v = np.array(np.arange(-np.pi, np.pi + (2*np.pi/points_circ), (2*np.pi/points_circ)))


		coords = np.resize(coords, [iz,points_ab,5,points_circ+1])
		#print '\nCoords\n'
		#print coords
		
		coords[iz-1,:,0,:] = short_r * np.outer(np.sin(u),np.cos(v))
		coords[iz-1,:,1,:] = short_r * np.outer(np.sin(u), np.sin(v))
		coords[iz-1,:,2,:] = long_r * np.outer(np.cos(u), np.ones([1,len(v)]))
		coords[iz-1,:,3,:] = np.transpose(np.tile(u,(len(v),1)))
		coords[iz-1,:,4,:] = np.tile(v,(len(u),1))




		iz = iz + 1

	#print np.shape(coords)
	#print coords

	dof_per_element=8
	el_nodes = np.zeros([dof_per_element,5])
	naive_elements = np.zeros([nab*ncirc*nr,dof_per_element],dtype=int)
	naive_nodes    = np.zeros([dof_per_element*nab*ncirc*nr,5])
	naive_epi1endo0 = np.zeros([dof_per_element*nab*ncirc*nr,1])

	ei = 0

	for el_ab in range(0, nab):
		for el_circ in range(0, ncirc):
			for el_r in range(0,nr):
				ix_ab = (el_ab-1) + np.array(range(1,3))
				ix_c = (el_circ-1) + np.array(range(1,3))
				ix_r = (el_r-1) + np.array(range(1,3))
				#print 'ix_ab: ', ix_ab
				#print 'ix_c: ', ix_c
				#print 'ix_r: ', ix_r

				for xi3 in range(0,2):
					for xi2 in range(0,2):
						for xi1 in range(0,2):			
							#print ix_ab[xi2], ix_c[xi1], ix_r[xi3]
							#print coords[ ix_r[xi3] , ix_ab[xi2] , :, ix_c[xi1] ]    
							#print 4*(xi3) + 2*(xi2) + xi1
							aux=coords[ ix_r[xi3] , ix_ab[xi2] , :, ix_c[xi1] ] 
							#print aux
							el_nodes[4*(xi3) + 2*(xi2) + xi1,:] = coords[ ix_r[xi3] , ix_ab[xi2] , :, ix_c[xi1] ]    


				ei = ei + 1
				nn = (ei-1)*8 + np.array(range(0,8))
				naive_elements[ei-1][:] = nn
				naive_nodes[nn,:] = el_nodes
				#print float(el_r)/nr
				#print float(el_r)/nr+(np.floor(np.array(range(0,8))/4.))/(nr)
				naive_epi1endo0[nn,0] = float(el_r)/nr + (np.floor(np.array(range(0,8))/4.))/nr

	#print naive_nodes
	#print el_nodes
	#merge nodes to create mesh
	
	Celems,Cnodes,nodemap,revmap = merge_nodes(naive_elements,naive_nodes[:,0:3])
	UV = naive_nodes[revmap,3:5]
	#print Celems
	
	
	epi1endo0 = np.zeros([np.size(Cnodes,0),1])

	for i in range(0,np.size(Cnodes,0)):
		f = naive_epi1endo0[nodemap==i]
		epi1endo0[i] = f[0]

	#print epi1endo0

	fib = np.zeros([np.size(Cnodes,0), 3])
	angle = np.zeros([np.size(Cnodes,0),1])

	for ni in range(0,np.size(Cnodes, 0)):
		u = UV[ni,0]
		v = UV[ni,1]
		fiber_angle = (fiber_angle_endo + epi1endo0[ni] * (fiber_angle_epi - fiber_angle_endo)) * (np.pi/180)

		long_r  = outer_long_axis - (1-epi1endo0[ni])*wall_thickness
		short_r = outer_short_axis - (outer_long_axis-long_r)


		deriv_dir = np.array([math.sin(fiber_angle), math.cos(fiber_angle)])[:, np.newaxis]

		#these are simply the d/du and d/dv in a matrix
		M = np.array( [ [short_r*math.cos(u)*math.cos(v), -short_r*math.sin(u)*math.sin(v)], [short_r*math.cos(u)*math.sin(v), short_r*math.sin(u)*math.cos(v)], [-long_r*math.sin(u), 0] ] )

		M[:,0] = M[:,0]/np.linalg.norm(M[:,0])
		M[:,1] = M[:,1]/np.linalg.norm(M[:,1])


		fib[ni, :] = (M.dot(deriv_dir)).T
		fib[ni, :] = fib[ni, :] / np.linalg.norm(fib[ni, :])

		angle[ni] = fiber_angle*(180/np.pi)
        
		#apex nodes:
		if np.abs(u+np.pi) < 1e-6:
			fib[ni,:] = 0
			angle[ni] = 0
			#print ni

	
	for i in range(len(Celems[:,2])):
		aux=0
		aux = Celems[i,2]
		Celems[i,2] = Celems[i,3]	
		Celems[i,3]	= aux

		aux=0
		aux = Celems[i,6]
		Celems[i,6] = Celems[i,7]	
		Celems[i,7]	= aux

	base = np.where(Cnodes[:,2]==5)[0]

	boundary_elem = []
	for i in range(np.size(Celems,0)):
		el = Celems[i,:]
		endocardio = Celems[i, np.where(epi1endo0[el]==0)[0] ]
		if(endocardio!=[]):
			aux, freq = np.unique(endocardio,return_counts=True)
			if(len(freq)<4):
				endocardio = np.array([endocardio[2], endocardio[1], endocardio[0],endocardio[3]])
			else:
				endocardio = np.array([endocardio[0], endocardio[3], endocardio[2],endocardio[1]])
			boundary_elem.append(endocardio)

	boundary_elem = np.array(boundary_elem)
	
	
	return Cnodes, Celems, fib, base, boundary_elem


def merge_nodes(elements, nodes):
	os = np.size(nodes,0)

	#print nodes

	nodes[:,[0,1,2]] = nodes[:,[2,1,0]]

	nodes= np.round(nodes,10)

	nodes+=0.


	uniques, reversenodemap, nodemap = np.unique([str(i) for i in nodes], return_index=True, return_inverse=True)
	nodes = nodes[reversenodemap]

	#print uniques
	ids = np.lexsort(np.transpose(nodes)[::-1])	
	nodes=nodes[ids]
	reversenodemap=reversenodemap[ids]

	l=[]
	for i in range(0,len(nodemap)):
		l.append( int(np.where(ids==nodemap[i])[0]) )
	nodemap=np.array(l)
	
	nodes[:,[0,1,2]] = nodes[:,[2,1,0]]


	elements = nodemap[elements]	

	print 'Removed duplicates %d -> %d\n' % (os,np.size(nodes,0))

	return elements, nodes, nodemap, reversenodemap




def carp2xml (Cnodes, Celems, fib, base, param, ninc, pressure, boundary_elem, outputMesh):

    
    
    pts  = Cnodes
    elem = Celems

    elem_xml = None
    num_dim  = None
    num_pts  = np.shape(pts)[0]
    num_elem = np.shape(elem)[0]

    num_dim, elem_xml = 3, "hexahedron"

    # print info
    print("Number of nodes: %d" % num_pts)
    print("Number of elements: %d" % num_elem)
    print("Element type: %s" % elem_xml)

    sparam = str(param)
    sparam=sparam.replace('[','')
    sparam=sparam.replace(']','')
    sparam=sparam.strip().split()
    print sparam
    #sparam=sparam.replace('   ',', ')

    # basic header
    outputFile = file(outputMesh, 'w')
    outputFile.write('<?xml version=\"1.0\"?>\n')    
    outputFile.write('<problem>\n')
    outputFile.write('<elasticity type=\"THREE_DIM\">\n')
    outputFile.write('  <parameters>\n')
    outputFile.write('    <material>Guccione</material>\n')
    outputFile.write('    <coefficients>')
    for s in sparam[0:-1]:
    	outputFile.write('%s, ' % (s))
    outputFile.write('%s' % (sparam[-1]))
    outputFile.write('</coefficients>\n')

    outputFile.write('    <ninc>%d</ninc>\n' % (ninc))
    outputFile.write('  </parameters>\n')
    outputFile.write('  <pressure>\n')
    outputFile.write('    <node id=\"1\" marker=\"11\" value=\"%f\"/>\n' % (pressure))
    outputFile.write('  </pressure>\n')
    outputFile.write('  <prescribed_displacement>\n')
    for i in range(0,len(base)):
        outputFile.write('    <node id=\"%d\" direction=\"0\" value=\"0.000000\" />\n' % (base[i]))
        outputFile.write('    <node id=\"%d\" direction=\"1\" value=\"0.000000\" />\n' % (base[i]))
        outputFile.write('    <node id=\"%d\" direction=\"2\" value=\"0.000000\" />\n' % (base[i]))

    outputFile.write('  </prescribed_displacement>\n')
    outputFile.write('</elasticity>\n')



    outputFile.write('<mesh celltype="%s" dim="%d">\n' % (elem_xml, num_dim))

    # nodes section
    outputFile.write('  <nodes size="%d">\n' % (num_pts))
    for i in xrange(num_pts):
        outputFile.write('    <node id="%d" ' % (i) ) 
        outputFile.write('x="%f" y="%f" z="%f" />\n' % (pts[i,0],pts[i,1],pts[i,2]))
    outputFile.write('  </nodes>\n')

    # elements section
    outputFile.write('  <elements size="%d">\n' % (num_elem))
    for i in xrange(num_elem):
        etype = "Hx"
        outputFile.write('    <element id="%d" ' % (i))
        write_element(outputFile, np.array(elem[i][:],dtype=int) )
    outputFile.write('  </elements>\n')



    fibtype = 1
    fibname = None

    vecs = fib
    fibsize = np.shape(vecs)[0]
    if(fibsize != num_elem):
        print("Error: number of fiber vectors must match the number of elements")
        sys.exit(1)
    if(fibtype==1):
        fibname = "fiber_transversely_isotropic"
    else:
        fibname = "fiber_orthotropic"
    print("Fiber model: %s" % fibname)
    print("Fiber vectors: %d" % fibsize)

    # element_data section
    outputFile.write('  <element_data type="%s">\n' % (fibname))
    for i in xrange(fibsize):
        outputFile.write('    <element id="%d">\n' % (i))
        write_vecs(outputFile, vecs[i,:])
        outputFile.write('    </element>\n')
    outputFile.write('  </element_data>\n')       

    outputFile.write('  <boundary celltype=\"quadrilateral\" dim=\"2\">\n')

    for i in range(np.size(boundary_elem,0)):
        outputFile.write('    <element id=\"%d\" marker=\"11\"  v0=\"%d\" v1=\"%d\" v2=\"%d\" v3=\"%d\" />\n' % (i,boundary_elem[i,0], boundary_elem[i,1], boundary_elem[i,2], boundary_elem[i,3]))

    outputFile.write('  </boundary>\n')       
    

    outputFile.write('</mesh>\n')
    outputFile.write('</problem>\n')


    ## electrophysiology parameters
    #outputFile.write('<electrophysiology>\n')
    #if( os.path.isfile(stimFile) ):
    #    sfile = open(stimFile, "r")
    #    nstim = int(sfile.readline())
    #    vstim = []
    #    for i in range(nstim):
    #        l = sfile.readline()
    #        vstim.append(l)
    #    sfile.close()
    #    outputFile.write('  <stimuli number=\"%d\">\n' % nstim)
    #    for i in range(nstim):
    #        write_stimuli(outputFile, vstim[i])
    #    outputFile.write('  </stimuli>\n')
    #outputFile.write('</electrophysiology>\n')
    #
    #outputFile.close()
    print("Done")

# ------------------------------------------------------------------------------    

def write_element(out, conec):
    num_nodes = len(conec)
    for i in range(num_nodes):
        out.write('v%d="%d" ' % (i,conec[i]))
    out.write(' />\n')

# ------------------------------------------------------------------------------

def write_vecs(out, vec):
    size = np.shape(vec)[0]
    # transversely isotropic - fiber only
    if(size == 3):
        out.write('        <fiber>%f,%f,%f</fiber>\n' % (vec[0],vec[1],vec[2]))
    elif(size == 9):
        out.write('        <fiber>%f,%f,%f</fiber>\n' % (vec[0],vec[1],vec[2]))
        out.write('        <sheet>%f,%f,%f</sheet>\n' % (vec[3],vec[4],vec[5]))
        out.write('        <normal>%f,%f,%f</normal>\n' % (vec[6],vec[7],vec[8]))
        # TERMINAR DE IMPLEMENTAR
        pass

# ------------------------------------------------------------------------------

def write_stimuli(out, s):
    stim = s.split()
    s  = float(stim[0])
    d  = float(stim[1])
    v  = float(stim[2])
    x0 = float(stim[3])
    x1 = float(stim[4])
    y0 = float(stim[5])
    y1 = float(stim[6])
    z0 = float(stim[7])
    z1 = float(stim[8])
    out.write('      ')
    out.write('<stim')
    out.write(' start=\"%.2f\" duration=\"%.2f\" value=\"%.4f\" ' % (s,d,v))
    out.write('x0=\"%.2f\" x1=\"%.2f\" ' % (x0,x1))
    out.write('y0=\"%.2f\" y1=\"%.2f\" ' % (y0,y1))
    out.write('z0=\"%.2f\" z1=\"%.2f\" ' % (z0,z1))
    out.write('/>\n')




def teste():
	top = 5.
	outer_long_axis = 20.
	outer_short_axis = 10.
	wall_thickness = 3.
	#fiber angles
	fiber_angle_epi  = -60.
	fiber_angle_endo = 90.

	mesh_params = [top, outer_long_axis, outer_short_axis, wall_thickness, fiber_angle_epi, fiber_angle_endo]

	#Number of elements: [nab, ncirc, ntrans]
	h=[12,27,6]
	#h=0.5

	Cnodes, Celems, fib, base, boundary_elem=benchmark_ellipse(h, mesh_params)


	model_params = np.array([2, 8, 2, 4, 300])
	#param = np.array([10, 1, 1, 1, 200])
	ninc=100
	pressure=15

	#Fiber nodes to fiber elements:
	fib_elem = []
	for e in Celems:
		f=np.sum(fib[e,:], axis=0)/8.
		fib_elem.append(f/np.linalg.norm(f))
	fib_elem=np.array(fib_elem)


	carp2xml (np.around(Cnodes,4), Celems, fib_elem, base, model_params, ninc, pressure, boundary_elem, "teste_prob3_fibra.xml")



def createXml(filename, mesh_params, model_params, h, ninc, pressure):

	Cnodes, Celems, fib, base, boundary_elem=benchmark_ellipse(h, mesh_params)

	#Fiber nodes to fiber elements:
	fib_elem = []
	for e in Celems:
		f=np.sum(fib[e,:], axis=0)/8.
		fib_elem.append(f/np.linalg.norm(f))
	fib_elem=np.array(fib_elem)

        out = file('teste_malha.pts', 'w')

        for p in np.around(Cnodes,4):
            out.write("%f %f %f\n" % (p[0], p[1], p[2]))

        out.close()

        out = file('teste_malha.elem', 'w')
        for e in Celems:
            out.write("Hx %d %d %d %d %d %d %d %d 1\n" % (e[0], e[1], e[2], e[3], e[4], e[5], e[6], e[7]))

        out.close()


        print np.size(np.around(Cnodes,4),0)
        print np.size(Celems,0)
#        print fib_elem


	carp2xml (np.around(Cnodes,4), Celems, fib_elem, base, model_params, ninc, pressure, boundary_elem, filename)

