#include "fem.h"
#include "malloc.h"

void mesh_read (get_mesh *msh, const char *filename, element_list *el)
{
  char nf[256]; // node file
  char ef[256]; // element file
  char af[256]; // cond. axes file
  FILE *fp;
  int i,retval;
  float a1, a2, a3, b1, b2, b3, c1, c2, c3;

  sprintf(nf, "%s.pts", filename);
  msh->nodefile = (char *) calloc(strlen(nf)+1, sizeof(char));
  sprintf(ef, "%s.elem", filename);
  msh->elemfile = (char *) calloc(strlen(ef)+1, sizeof(char));
  sprintf(af, "%s.fib", filename);
  msh->axesfile = (char *) calloc(strlen(af)+1, sizeof(char));

  strcpy(msh->nodefile, nf);
  strcpy(msh->elemfile, ef);
  strcpy(msh->axesfile, af);

  mesh_read_elements (msh);
  mesh_read_nodes    (msh);

  // allocate element list
  el->elems = (element *) malloc(msh->nel * sizeof(element));

  // set mesh properties
  msh->nx = (int) sqrt(msh->ndof);
  msh->ny = (int) sqrt(msh->ndof);

  // read conductivity vectors from file
  element *elem;
  fp = fopen(af, "r");
  if(!fp){
    fprintf(stderr,"Error opening fibers file %s \n", af);
    exit(1);
  }

  for(i=0; i<msh->nel; i++)
  {
    retval = fscanf(fp, "%f", &a1);
    retval = fscanf(fp, "%f", &a2);
    retval = fscanf(fp, "%f", &a3);
    retval = fscanf(fp, "%f", &b1);
    retval = fscanf(fp, "%f", &b2);
    retval = fscanf(fp, "%f", &b3);
    retval = fscanf(fp, "%f", &c1);
    retval = fscanf(fp, "%f", &c2);
    retval = fscanf(fp, "%f", &c3);
    
    // build element info
    elem = &el->elems[i];
    elem->elid   = i;
    elem->eltype = msh->eltype;
    elem->ndim   = msh->ndim;
    elem->nint   = msh->nint;
    elem->ndofs  = msh->nen;

    elem->al[0] = a1;
    elem->al[1] = a2;
    elem->al[2] = a3;

    elem->at[0] = b1;
    elem->at[1] = b2;
    elem->at[2] = b3;

    elem->an[0] = c1;
    elem->an[1] = c2;
    elem->an[2] = c3;

  }

  fclose(fp);

  cout << endl;
  cout << "Mesh information" << endl << endl;
  cout << " Node file:                " <<  msh->nodefile << endl;
  cout << " Element file:             " << msh->elemfile << endl;
  cout << " Number of nodes:          " << msh->ndof << endl;
  cout << " Number of elements:       " << msh->nel << endl;
  cout << flush;
  
}

void mesh_read_elements (get_mesh *msh)
{
  int i,**conec, tmp, ret, nel;
  int n1,n2,n3,n4;
  char strElemType[2];
  FILE *fp;

  // open elemfile
  fp = fopen(msh->elemfile,"r");
  if(!fp){
    error("Error opening file in read_elements\n");
  }

  // read element conectivities
  ret = fscanf(fp,"%d\n",&nel);
  msh->nel = nel;

  ret = fscanf(fp,"%s ",strElemType);

  if( strcmp( strElemType, "Tr" ) == 0 ) {
    msh->eltype = 1;
    msh->nen    = 3;
    msh->nint   = 3;
    msh->ndim   = 2;
  }
  if( strcmp( strElemType, "Qd" ) == 0 ) {
    msh->eltype = 0;
    msh->nen    = 4;
    msh->nint   = 4;
    msh->ndim   = 2;
  }
  if( strcmp( strElemType, "Hx" ) == 0 ) {
    msh->eltype = 2;
    msh->nen    = 8;
    msh->nint   = 8;
    msh->ndim   = 3;
  }

  conec = new_2D_array<int> (msh->nen, msh->nel);

  if(msh->nen==3){

    ret = fscanf(fp,"%d %d %d %d\n",&n1,&n2,&n3,&tmp);
    conec[0][0] = n1; 
    conec[1][0] = n2; 
    conec[2][0] = n3;
    for(i=1; i<msh->nel; i++){
      ret = fscanf(fp,"%d %d %d\n",&n1,&n2,&n3);
      conec[0][i] = n1; 
      conec[1][i] = n2; 
      conec[2][i] = n3;
    }

  } else if (msh->nen==4) {

    ret = fscanf(fp,"%d %d %d %d %d\n",
		    &n1,&n2,&n3,&n4,&tmp);
    conec[0][0] = n1; 
    conec[1][0] = n2; 
    conec[2][0] = n3;  
    conec[3][0] = n4;
    for(i=1; i<msh->nel; i++){
      ret = fscanf(fp,"%s %d %d %d %d %d\n",
		       strElemType,&n1,&n2,&n3,&n4,&tmp);
      conec[0][i] = n1; 
      conec[1][i] = n2; 
      conec[2][i] = n3; 
      conec[3][i] = n4;
    }
  }

  // close elemfile
  fclose(fp);

  // set mesh localToGlobal
  msh->localToGlobal = conec;
}

void mesh_read_nodes (get_mesh *msh)
{
  int i, j, node, ret, ndof;
  double xtemp[3];
  Real **x;
  FILE *fp;

  // open nodefile
  fp = fopen(msh->nodefile,"r"); 
  if(!fp){
    fprintf(stderr,"Error opening nodes file %s in read_nodes()\n",
	    msh->nodefile);
    exit(1);
  }

  // read nodes
  ret = fscanf(fp,"%d\n",&ndof); 
  msh->ndof = ndof;

  x = new_2D_array<Real> (msh->ndim, msh->ndof);
  node = 0;
  for(i=0; i<msh->ndof; i++){
    ret = fscanf(fp, "%lf %lf %lf\n",&xtemp[0],&xtemp[1],&xtemp[2]);
    for(j=0; j<msh->ndim; j++){
      x[j][node] = (Real) xtemp[j];
    }
    node++;
  }
  
  // close file
  fclose(fp);

  // set the mesh properties
  msh->nodes = x;
}

/*!
 Read boundary nodes and values
 @param[in]  bconfile   The input file in to read boundary values
 @param[out] nbc        The number of boundary conditions to be applied
 */

Real ** mesh_read_bc (char bconfile[], int *nbc)
{
  int i, bcnode, ret;
  FILE *fp;
  Real **bc;
  double bcvalue;
  fp = fopen(bconfile,"r");
  if(!fp){
    fprintf(stderr," Error opening boundary file %s in read_bcon()\n",bconfile);
    exit(1);
  }

  // read boundary nodes
  ret = fscanf(fp,"%d\n",nbc);
  bc = new_2D_array<Real> (*nbc,2);
  for(i=0; i<*nbc; i++){
    ret = fscanf(fp,"%d %lf\n",&bcnode,&bcvalue);
    bc[i][0] = bcnode;
    bc[i][1] = (Real) bcvalue;
  }
  fclose(fp);
  return bc;
}

void delete_mesh (get_mesh *msh)
{
  // free local to global dof mapping
  if(msh->localToGlobal != NULL)
  {
    delete_2D_array<int> (msh->localToGlobal, msh->nen, msh->nel);
  }

  // free nodes coordinates
  if(msh->nodes != NULL)
  {
    delete_2D_array<Real> (msh->nodes, msh->ndim, msh->ndof);
  }
 
  // free file strings
  if(msh->nodefile != NULL) free(msh->nodefile);
  if(msh->elemfile != NULL) free(msh->elemfile);
  if(msh->axesfile != NULL) free(msh->axesfile);

}

void delete_element_list (element_list *el)
{
  if(el->elems != NULL) free(el->elems);
}

