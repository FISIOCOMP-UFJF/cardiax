#include <stdio.h>
#include <stdlib.h>
#include <boost/program_options/option.hpp>
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>
	 
using namespace std;
namespace po = boost::program_options;

int main(int argc, char *argv[])
{
  // ARGUMENTS CHECK
  if(argc < 4){
    fprintf(stderr, 
     "Usage: %s <width-pts> <height-pts> <element-side-size> <flag-triangle>\n", 
	    argv[0]);
    return 1;		
  }

  // <flag-triangle>
  // = 0 - False
  // = 1 - True - Use Triangular mesh
	
  int width, height, flag=0;
  float elem_size;

  FILE *felem, *fpts, *ffib, *fstim, *fbdn;

  width = atoi(argv[1]);
  height = atoi(argv[2]);
  elem_size = atof(argv[3]);

  if (argc == 5) flag = atoi(argv[4]);

  char elem_filename[256];
  char fib_filename[256];
  char pts_filename[256];
  char bdn_filename[256];
  char stim_filename[256];

  if (flag) {
    sprintf(elem_filename,"tissue_tri_%dx%d_%.2f_um.elem", width, height, elem_size);
    sprintf(fib_filename,"tissue_tri_%dx%d_%.2f_um.fib", width, height, elem_size);
    sprintf(pts_filename,"tissue_tri_%dx%d_%.2f_um.pts", width, height, elem_size);
    sprintf(stim_filename,"tissue_tri_%dx%d_%.2f_um.stim", width, height, elem_size);
  } else {
    sprintf(elem_filename,"tissue_%dx%d_%.2f_um.elem", width, height, elem_size);
    sprintf(fib_filename,"tissue_%dx%d_%.2f_um.fib", width, height, elem_size);
    sprintf(pts_filename,"tissue_%dx%d_%.2f_um.pts", width, height, elem_size);
    sprintf(bdn_filename,"tissue_%dx%d_%.2f_um.bdn", width, height, elem_size);
    sprintf(stim_filename,"tissue_%dx%d_%.2f_um.stim", width, height, elem_size);
  }

  felem = fopen(elem_filename, "w");
  ffib  = fopen(fib_filename, "w");
  fpts  = fopen(pts_filename, "w");
  fbdn  = fopen(bdn_filename, "w");
  fstim = fopen(stim_filename, "w");
	
  // writing elements file
  int num_elem;
  if (!flag){
    // use quadrilaterals
    num_elem = (width-1)*(height-1);
    fprintf(felem, "%d\n", num_elem);
    
    for(int i = 0; i < (height -1); i++){
      for(int j = 0; j< width -1; j++){
	fprintf(felem, "Qd %d %d %d %d 1\n", 
		i*width + j, 
		i*width + j + 1, 
		(i+1)*width + j + 1, 
		(i+1)*width + j);
      }
    }
  } else {
    // use triangles
    num_elem = 2*(width-1)*(height-1);
    fprintf(felem, "%d\n", num_elem);
    
    for(int i = 0; i<(height-1); i++){
      for(int j = 0; j<(width-1); j++){
	fprintf(felem, "Tr %d %d %d 1\n", 
		i*width + j, 
		i*width + j + 1, 
		(i+1)*width + j + 1);
	fprintf(felem, "Tr %d %d %d 1\n", 
		i*width + j, 
		(i+1)*width + j + 1, 
		(i+1)*width + j);
      }
    }
  }

  float a1, a2, a3, b1, b2, b3, c1, c2, c3;

  // writing fib file
  fprintf(ffib, "%d\n", 0);
  for(int i=0; i<num_elem; i++){
    /*a1 = 0.0;
    a2 = 0.0;
    a3 = -1.0;
    b1 = 0.707107;
    b2 = 0.707106;
    b3 = 0.0;
    c1 = 0.707107;
    c2 = -0.707106;
    c3 = 0.0;*/
    a1 = 1.0;
    a2 = 0.0;
    a3 = 0.0;
    b1 = 0.0;
    b2 = 1.0;
    b3 = 0.0;
    c1 = 0.0;
    c2 = 0.0;
    c3 = 1.0;
    fprintf(ffib, "%f %f %f %f %f %f %f %f %f\n", 
	    a1, a2, a3, b1, b2, b3, c1, c2, c3);
  }
	
  // writing pts file
  fprintf(fpts, "%d\n", width*height);
  for(int i=0; i<height; i++){
    for(int j=0; j<width; j++){
      fprintf(fpts, "%f %f %f\n", j*elem_size, i*elem_size, 0.0);
    }
  }

  // writing bdn file
  fprintf(fbdn, "%d\n", width+height+(width-2)+(height-2));

  for(int i=0; i<width; i++){ // bottom
    fprintf(fbdn, "%d %f\n", i, 0.0);
  }  

  for(int j=1; j<height-1; j++){ // right
    fprintf(fbdn, "%d %f\n", (j*width)+(width-1), 0.0);
  }  

  for(int j=1; j<height-1; j++){ // left
    fprintf(fbdn, "%d %f\n", j*width, 0.0);
  }  

  for(int i=0; i<width; i++){ // top
    fprintf(fbdn, "%d %f\n", (height-1)*width+i, 0.0);
  }  

	
  // writing stim file
  int stim_dimx = 0.1*width*elem_size;
  int stim_dimy = height*elem_size;

  // plane wave starting at left
  int stim_startx = 0;
  int stim_starty = 0;

  // central stimulus
  //int stim_startx = elem_size*(width/2 - stim_dimx/2);
  //int stim_starty = elem_size*(height/2 - stim_dimy/2);

  int stim_itstart = 1;
  float stim_amplitude = -8.0;

  fprintf(fstim, "1\n");
  fprintf(fstim, "%d %g %d %d %d %d\n", 
	  stim_itstart, stim_amplitude, 
	  stim_startx, stim_startx + stim_dimx, 
	  stim_starty, stim_starty + stim_dimy);

  fclose(felem);
  fclose(ffib);
  fclose(fpts);
  fclose(fstim);
  fclose(fbdn);

  return 0;
}
